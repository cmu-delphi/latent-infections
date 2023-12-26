# standard
import json
from ctypes import *
from datetime import date, timedelta
from functools import partial
from typing import Callable, Dict

# third party
import numpy as np
from scipy.linalg import toeplitz
from scipy.sparse import diags as band

# load DP algorithm
so_file = "dp_1d.so"
c_dp_1d = CDLL(so_file)


def construct_day_specific_convolution_matrix(y, run_date, delay_kernels):
    """Construct convolution matrix with time-varying rows."""
    n = y.shape[0]
    if n > len(delay_kernels):
        assert ValueError('Not enough dates in delay_kernel dictionary')

    kernel_length = len(delay_kernels[run_date])
    C = np.zeros((n, n))
    first_kernel_date = min(delay_kernels.keys())
    for i in range(C.shape[0]):
        kernel_day = max(run_date - timedelta(i), first_kernel_date)
        end_index = max(0, C.shape[1] - i)
        start_index = max(0, end_index - kernel_length)
        day_specific_kernel = np.array(delay_kernels[kernel_day][::-1])
        row = C.shape[0] - i - 1
        if end_index > 0:
            C[row, start_index:end_index] = day_specific_kernel[
                                            -(end_index - start_index):]
    return C, np.array(delay_kernels[run_date])



def _impute_with_neighbors(x: np.ndarray) -> np.ndarray:
    """
    Impute missing values with the average of the elements immediately
    before and after.
    Parameters
    ----------
    x
        signal with missing values
    Returns
    -------
        imputed signal
    """
    # handle edges
    if np.isnan(x[0]):
        x[0] = x[1]

    if np.isnan(x[-1]):
        x[-1] = x[-2]

    imputed_x = np.copy(x)
    for i, (a, b, c) in enumerate(zip(x, x[1:], x[2:])):
        if np.isnan(b):
            imputed_x[i + 1] = (a + c) / 2

    assert np.isnan(imputed_x).sum() == 0

    return imputed_x



def deconvolve_tf(y: np.ndarray,
                  C: np.ndarray,
                  lam: float,
                  n_iters: int = 10000,
                  k: int = 3,
                  clip: bool = False,
                  alpha_k = None,
                  u_k = None,
                  tolerance = 1e-3) -> np.ndarray:
    """
        Perform deconvolution with trend filtering smoothing
        Parameters
        ----------
            y
               array of values to convolve
            C
               convolution matrix operator
            lam
               trend filtering parameter
            n_iters
               number of ADMM interations to perform
            k
               order of the trend filtering penalty
            clip
               Boolean to clip count values to [0, infty)
            alpha_k
               starting alpha_k value if using a warm start
            u_k
               starting u_k value if using a warm start
            tolerance
               termination threshold for checking against the primal and residual duals
        Returns
        -------
           array of the deconvolved signal values
           final alpha_k value
           final u_k value
    """

    n = y.shape[0]
    C = C[:n]
    rho = lam  # set equal
    D = band([-1, 1], [0, 1], shape=(n - 1, n)).toarray()
    D = np.diff(D, n=k, axis=0)

    # pre-calculations
    Cty = C.T @ y
    first_x_update = np.linalg.inv(C.T @ C + rho * D.T @ D)

    if alpha_k is None:
        alpha_k = np.zeros(n - k - 1)
    if u_k is None:
        u_k = np.zeros(n - k - 1)

    x_k = first_x_update @ (Cty + rho * D.T @ (alpha_k + u_k))

    niter = 0
    for i in range(n_iters):
        niter += 1
        alpha_old = alpha_k # Store previous value of alpha for convergence checking

        x_k = first_x_update @ (Cty + rho * D.T @ (alpha_k + u_k))

        c_dp_1d.read.argtypes = c_int, POINTER(c_double), c_double, POINTER(
            c_double)
        c_dp_1d.read.restype = None
        x = D @ x_k - u_k
        alpha_k = (c_double * len(x))()
        x_c = (c_double * len(x))(*x)
        c_dp_1d.tf_dp(c_int(len(x)), x_c, c_double(lam/rho), alpha_k) #%%
        u_k = u_k + alpha_k - D @ x_k

        # Compute primal and dual residual norms
        rr = np.linalg.norm(D @ x_k - alpha_k)
        ss = np.linalg.norm(rho * D.T @ ([alpha_k[i]-alpha_old[i] for i in range(len(alpha_k))]))

        # Check convergence criterion
        if rr/np.sqrt(n) < tolerance and ss/np.sqrt(n) < tolerance:
            break

    if clip:
        x_k = np.clip(x_k, 0, np.infty)
    return x_k, alpha_k, u_k



def deconvolve_tf_cv(y: np.ndarray,
                     kernel_dict: Dict,
                     as_of_date: date,
                     fit_func: Callable = deconvolve_tf,
                     lam_cv_grid: np.ndarray = np.logspace(1, 3.5, 10),
                     n_iters: int = 10000,
                     k: int = 3,
                     clip: bool = False) -> np.ndarray:
    """
       Run cross-validation to tune smoothness over deconvolve_tf.
       First, leave-every-third-out CV is performed over lambda
       and the lambda with the smallest squared error is chosen.

       Parameters
       ----------
       y
           array of values to convolve
       kernel_dict
           dictionary of convolution kernel values, indexed by as_of dates
       fit_func
           deconvolution function to use
       lam_cv_grid
           grid of trend filtering penalty values to search over
       n_iters
           number of ADMM interations to perform.
       k
           order of the trend filtering penalty.
       clip
           Boolean to clip count values to [0, infty)
       output_tuning
           Boolean whether to output a file storing metadata on tuning
           parameters
       output_tuning_file
           name of output file, not used if None or output_tuning is False
       Returns
       -------
           array of the deconvolved signal values
    """
    fit_func = partial(fit_func, n_iters=n_iters, k=k, clip=clip)
    n = y.shape[0]
    lam_cv_loss = np.zeros((lam_cv_grid.shape[0],))

    # use leave-every-third-out cv for finding lambda, this controls smoothness
    # of entire curve
    for i in range(3):
        test_split = np.zeros((n,), dtype=bool)
        test_split[i::3] = True

        # Handle the largest lambda in lam_cv_grid first (outside of a loop)
        x_hat = np.full((n,), np.nan)
        C, kernel = construct_day_specific_convolution_matrix(
            y[~test_split],
            as_of_date,
            kernel_dict)
        x_hat[~test_split], alpha_w, u_w = fit_func(y=y[~test_split],
                                      lam=max(lam_cv_grid),
                                      C=C)
        x_hat = _impute_with_neighbors(x_hat)
        C1, _ = construct_day_specific_convolution_matrix(x_hat, as_of_date,
                                                          kernel_dict)
        y_hat = (C1 @ x_hat)[:len(x_hat)]
        lam_cv_loss[0] += np.sum((y[test_split] - y_hat[test_split]) ** 2)

        # Handle the rest of the lambdas
        for j, reg_par in enumerate(reversed(lam_cv_grid[:-1]), start = 1):
            x_hat = np.full((n,), np.nan)
            x_hat[~test_split], alpha_w, u_w = fit_func(y=y[~test_split],
                                          lam=reg_par,
                                          C=C, alpha_k = alpha_w, u_k = u_w) # use the previous alpha_w and u_w and update them
            x_hat = _impute_with_neighbors(x_hat)
            C1, _ = construct_day_specific_convolution_matrix(x_hat, as_of_date,
                                                              kernel_dict)
            y_hat = (C1 @ x_hat)[:len(x_hat)]
            lam_cv_loss[j] += np.sum((y[test_split] - y_hat[test_split]) ** 2)
    lam = lam_cv_grid[np.argmin(lam_cv_loss[::-1])] # Work on reversed lam_cv_loss because the input for that is the lambda from smallest to largest
    C, kernel = construct_day_specific_convolution_matrix(y, as_of_date,
                                                           kernel_dict)
    x_hat, _, _ = fit_func(y=y, lam=lam, C=C)

    return x_hat




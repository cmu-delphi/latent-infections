#' Leave-kth-out cross validation for choosing a optimal parameter lambda
#'
#' @param nfold Integer. This number of folds to conduct the leave-kth-out
#' cross validation. For leave-kth-out cross validation, every kth
#' observed_counts and their corresponding position (evenly or unevenly
#' spaced) are placed into the same fold. The first and last observed_counts are
#' not assigned to any folds. Smallest allowable value is `nfold = 2`.
#' @param error_measure Metric used to calculate cross validation scores.
#' Must be choose from `mse` and `mae`
#' `mse` calculates the mean square error; `mae` calculates the mean absolute error
#' @param regular_splits Logical.
#'   If `TRUE`, the folds for k-fold cross-validation are chosen by placing
#'   every kth point into the same fold. The first and last points are not
#'   included in any fold and are always included in building the predictive
#'   model. As an example, with 15 data points and `kfold = 4`, the points are
#'   assigned to folds in the following way:
#'   \deqn{
#'   0 \; 1 \; 2 \; 3 \; 4 \; 1 \; 2 \; 3 \;  4 \; 1 \; 2 \; 3 \; 4 \; 1 \; 0
#'   }{0 1 2 3 4 1 2 3 4 1 2 3 4 1 0} where 0 indicates no assignment.
#'   Therefore, the folds are not random and running `estim_path_single()` twice
#'    will give the same result.
#' @param invert_splits Logical.
#'   Typical K-fold CV would use K-1 folds for the training
#'   set while reserving 1 fold for evaluation (repeating the split K times).
#'   Setting this to true inverts this process, using a much smaller training
#'   set with a larger evaluation set. This tends to result in larger values
#'   of `lambda` that minimize CV.
#' @param ... additional parameters passed to `estim_path_single()` function

#' @return An object with S3 class `"cv_tf"`. Among the list components:
#' * `full_fit` A list with the `thetas`, `lambda`, `korder`, fitted with all
#' `observed_counts` and `lambda`
#' * `cv_scores` leave-kth-out cross validation scores
#' * `cv_se` leave-kth-out cross validation standard error
#' * `lambda.min` lambda which achieved the optimal cross validation score
#' * `lambda.1se` lambda that gives the optimal cross validation score
#' within one standard error.
#' * `lambda` the value of `lambda` used in the algorithm.
#' @export
#'
#' @examples
#' y <- c(1, rnorm(100, dnorm(1:100, 50, 15) * 500 + 1))
#' cv <- cv_estimate_tf(y, korder = 3, nfold = 3, nsol = 30)
#' cv
#' 
#' params <- read_rds("data/ca-convolution-mat.rds")
#' y <- read_rds("data/deconvolved_ca_ga.rds") |>
#' filter(time_value <= ymd("2023-03-01"), geo_value == "ca") |>
#'   pull(cases)
#' cmats <- params$Cmat
#' cmat <- reduce(cmats, `+`)
#'
#' res <- cv_estimate_tf(y, 
#'                cmat = cmat,
#'                error_measure = "mse")
#' # Plot of the thetas for all lambdas
#' matplot(res$full_fit$thetas, type = "l")
#' # Plot line for lambda_min
#' matplot(res$full_fit$thetas[, which(res$lambda == res$lambda.min)], type = "l")
#' 
cv_estimate_tf <- function(
    observed_counts,
    x = 1:n,
    cmat,
    korder = 3L,
    nfold = 3L,
    error_measure = c("mse", "mae"),
    lambda = double(50),
    lambda_max = -1,
    maxiter = 1e3L,
    regular_splits = TRUE, 
    invert_splits = FALSE,
    ...) {
  arg_is_pos_int(nfold)
  n <- length(observed_counts) 
  xin <- x
  if (inherits(xin, "Date")) x <- as.numeric(x)
  arg_is_numeric(x)
  
  if (nfold == 1) cli::cli_abort("nfold must be greater than 1")
  
  ## Run program one time to create lambda 
  Rcpp::sourceCpp("src/estim_path.cpp") 
  
  full_fit = estim_path_single(observed_counts, x, cmat, korder = 3, lambda = lambda, nsol = length(lambda), lambdamax = lambda_max, maxadmm_iter = maxiter) 
  
  if (is.null(lambda)) lambda <- full_fit$lambda 
  
  if (regular_splits) { 
    middle_fold <- rep_len(1:nfold, n - 2)
  } else {
    middle_fold <- sample.int(nfold, n - 2, replace = TRUE)
  }
  foldid <- c(0, middle_fold, 0)
  cvall <- matrix(NA, nfold, length(lambda))
  
  error_measure <- match.arg(error_measure)
  err_fun <- switch(error_measure,
                    mse = function(y, m) (y - m)^2,
                    mae = function(y, m) abs(y - m)
  )
  
  for (f in 1:nfold) {
    train_idx <- if (invert_splits) foldid %in% c(0, f) else foldid != f
    test_idx <- !train_idx

    mod <- estim_path_single(observed_counts[train_idx], x, cmat[train_idx,], korder = 3, lambda = lambda, nsol = length(lambda), maxadmm_iter = maxiter)
    
    score <- colMeans(err_fun(observed_counts[test_idx], mod$thetas[test_idx,])) 
    cvall[f, seq_along(score)] <- score
    if (length(lambda) != length(score)) {
      cli::cli_warn(c(
        "Estimation for the full `lambda` sequence did not occur for fold {.val {f}}",
        "because the maximum number of iterations was exhausted.",
        "i" = "You may wish to increase `maxiter` from the current {.val {maxiter}}."
      ))
    }
  }
  index <- apply(cvall, 2, function(x) any(is.na(x)))
  cvall <- cvall[, !index]
  lambda <- lambda[!index]
  
  
  ### Calculate CV summary
  cv_scores <- colMeans(cvall)
  cv_se <- apply(cvall, FUN = stats::sd, MARGIN = 2) / sqrt(nfold)
  i0 <- which.min(cv_scores)
  
  structure(
    enlist(
      full_fit,
      cv_scores,
      cv_se,
      lambda,
      lambda.min = lambda[i0],
      lambda.1se = max(
        lambda[cv_scores <= cv_scores[i0] + cv_se[i0]],
        na.rm = TRUE
      ),
      call = match.call()
    ),
    class = "cv_tf"
  )
}

############################################################################################################

# Test using CA data:

# params <- read_rds("data/ca-convolution-mat.rds")
# y <- read_rds("data/deconvolved_ca_ga.rds") |>
#  filter(time_value <= ymd("2023-03-01"), geo_value == "ca") |>
#  pull(cases)
# x <- seq_along(y)
# cmats <- params$Cmat
# cmat <- reduce(cmats, `+`)

# res <- cv_estimate_tf(y, 
#               cmat = cmat,
#               error_measure = "mse")
# Plot of the thetas for all lambdas
# matplot(res$full_fit$thetas, type = "l")
# Plot line for lambda_min
# matplot(res$full_fit$thetas[,which(res$lambda == res$lambda.min)], type = "l")

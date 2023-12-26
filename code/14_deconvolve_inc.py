from config_spec_to_ablation import *

#######################################################################################################################

for state in states:
    # Load inc pmfs    # This assumes that the report date = symptom onset date.
    inc_dir = '/Users/admin/Downloads/inc_dist/data/incubation_period_distributions/' + str(state)
    inc_pmfs = pickle.load(open(f'{inc_dir}/{state}_inc_distribution_nov5.p', 'rb'))

    # Get JHU 7 day averages of the number of new confirmed COVID-19 cases per 100,000 population for plotting
    jhu_csse_7dav_data_state = jhu_csse_7dav_data[jhu_csse_7dav_data['geo_value'] == state]

    # Get y data for state
    jhu_csse_data_state = jhu_csse_data[jhu_csse_data['geo_value'] == state]

    y = jhu_csse_data_state.value.values

    n = y.shape[0]
    D = band([-1, 1], [0, 1], shape=(n - 1, n)).toarray()
    D = np.diff(D, n = 3, axis=0)

    # Get Phat matrix to calculate lambda_max
    C, _ = construct_day_specific_convolution_matrix(y, run_date, inc_pmfs)

    # Calculate lambda_max
    lambda_max = max(np.linalg.pinv(D @ D.T) @ D @ C.T @ y)

    # Create lambda grid, assuming the values are ordered from smallest to largest
    lam_grid = np.logspace(np.log10(10e-15 *lambda_max), np.log10(lambda_max), 20)

    # Deconvolution tf w/three fold cv
    res = deconvolve_tf_cv(y = y,
                     kernel_dict = inc_pmfs,
                     as_of_date = run_date,
                     lam_cv_grid = lam_grid,
                           clip = True)

    # Make folder for state if it doesn't exist
    if not os.path.exists(storage_dir + "/" + state):
        os.makedirs(storage_dir + "/" + state)

    # Store deconvolution results in folder for state
    np.save(f'{storage_dir}/{state}/{state}_inc_deconvolution_res_d60c_nov5.npy', res)

    # Save a plot comparing the tf results to the jhu incident case counts
    plt.figure(figsize=(12, 5))
    plt.scatter(jhu_csse_7dav_data_state.time_value,
                jhu_csse_7dav_data_state.value,
                color='gray', s=0.3, label='JHU 7 day avg cases')
    plt.plot(inc_pmfs.keys(), res, ls='--', label='Est. incident infections (assuming report = symptom onset)')
    plt.xlabel('Time')
    plt.ylabel('New infections in ' + state)
    plt.legend(loc='upper left')
    plt.show()
    plt.savefig(storage_dir + "/" + str(state) + "/" + str(state) + "_inc_deconvolution_plot_d60c_nov5_w_7dav.png")
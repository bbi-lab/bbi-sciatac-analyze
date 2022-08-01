#! /usr/bin/env python3

from __future__ import print_function
from __future__ import division
from scipy import stats
import numpy as np
import collections
import pandas as pd
import argparse
import json
from scipy.cluster.vq import whiten, kmeans2

MixedModelParams = collections.namedtuple("MixedModelParams",
                              ["mu_noise", "alpha_noise", "mu_signal", "alpha_signal", "frac_noise"])

def convert_params(mu, theta):
    """
    Convert mean/dispersion parameterization of a negative binomial to the ones scipy supports

    See https://en.wikipedia.org/wiki/Negative_binomial_distribution#Alternative_formulations
    """
    r = theta
    var = mu + 1 / r * mu ** 2
    p = (var - mu) / var
    return r, 1 - p


def pmf(counts, mu, theta):
    return stats.nbinom.pmf(counts, *convert_params(mu, theta))

def logpmf(counts, mu, theta):
    return stats.nbinom.logpmf(counts, *convert_params(mu, theta))


def cdf(counts, mu, theta):
    return stats.nbinom.cdf(counts, *convert_params(mu, theta))


def sf(counts, mu, theta):
    return stats.nbinom.sf(counts, *convert_params(mu, theta))

def weighted_mean(counts, weights=None):
    """The mean of the input counts.  If weights are provided, returns the weighted mean."""
    if weights is None:
        return counts.mean()

    return (counts * weights).sum() / weights.sum()


def weighted_variance(counts, weights=None):
    """The variance of the input counts.  If weights are provided, returns the weighted variance."""
    if weights is None:
        return counts.var()
    sum_of_squares = (weighted_mean(counts, weights) - counts) ** 2
    return (sum_of_squares * weights).sum() / weights.sum()


def MOM_dispersion(counts, weights=None):
    """Simple --- and biased --- estimator of the dispersion parameter of a negative binomial distribution
    using the method of moments.
    """
    mu = weighted_mean(counts, weights)
    var = weighted_variance(counts, weights)
    alpha = (var - mu) / np.power(mu, 2)
    return max(1e-6, alpha)


def MLE_dispersion(counts, weights=None, maxiter=10000, epsilon=1e-6):
    """Borrowed from RiboPy (https://bitbucket.org/caethan/ribopy/src):
    Fixed point estimation of the dispersion (alpha) parameter for a negative binomial distribution.
    Provides the maximum likelihood estimate for alpha.
    """

    v = np.arange(max(counts))
    def G_mle(alpha, mu, counts, weights=None):
        """This is the fixed point for the maximum likelihood estimate of alpha (the NB dispersion parameter).

        alpha - float in (0, inf), the dispersion parameter
        mu - float, the fitted mean of the sample data
        counts - array of ints of length N, where N is the number of samples

        returns a single float that is a new estimate of alpha
        """
        subfunc_cumsum = np.cumsum((alpha * v) / (1 + alpha * v))

        indices = counts - 1
        indices[indices < 0] = 0

        return np.log(1 + alpha * mu) / (mu - weighted_mean(subfunc_cumsum[indices], weights))

    if (counts <= 0).all():
        return 1e-6

    alpha = MOM_dispersion(counts, weights)

    # If max(counts) is too big, this will be slow; fall back on method-of-moments
    if max(counts) > 1e8:
        return alpha

    mu = weighted_mean(counts, weights)
    var = weighted_variance(counts, weights)
    if var <= mu:
        return alpha
    
    for i in range(maxiter):
        new = G_mle(alpha, mu, counts, weights)
    
        if 2 * abs(new - alpha) / (new + alpha) < epsilon:
            break
        alpha = new

    return new

def robust_divide(a,b):
    """Handles 0 division and conversion to floats automatically
    """
    a = float(a)
    b = float(b)
    if b == 0:
        return float('NaN')
    else:
        return a/b

def goodness_of_fit(count_data, params):
    """Goodness of fit metric using KS distance

    Count data is a simple array or list of counts
    """
    def mixed_nbinom_cdf(k):
        noise_cdf = stats.nbinom.cdf(k, 1 / params.alpha_noise, 1 / (1 + params.alpha_noise * params.mu_noise))
        signal_cdf = stats.nbinom.cdf(k, 1 / params.alpha_signal, 1 / (1 + params.alpha_signal * params.mu_signal))
        return params.frac_noise * noise_cdf + (1 - params.frac_noise) * signal_cdf

    gof, _ = stats.kstest(count_data, mixed_nbinom_cdf)
    return gof


def calculate_probability_distribution(count_dict, params):
    """Returns the normalized probability that each count comes from the signal or the noise distributions,
    given some parameters for those.  Counts where both PMFs are zero are returned with a (non-normalized)
    zero probability of coming from either."""
    count_values = np.array(sorted([k for k in count_dict]))

    noise_pmf = pmf(count_values, params.mu_noise, params.alpha_noise)
    signal_pmf = pmf(count_values, params.mu_signal, params.alpha_signal)


    output = np.vstack((params.frac_noise * noise_pmf, (1 - params.frac_noise) * signal_pmf))
    mask = output.sum(axis=0) == 0
    output[:, mask] = 0.5
    output /= np.sum(output, axis=0)
    output[:, mask] = 0.0
    return output

def weighted_parameter_estimates(count_dict, prob_matrix):
    """Given N counts and a 2 by N probability matrix with the likelihood that each count is from signal or noise
    distributions, estimate the maximum likelihood parameters for the distributions"""
    prob_noise = prob_matrix[0, :]
    prob_signal = prob_matrix[1, :]

    noise_mask = prob_noise > 0
    signal_mask = prob_signal > 0

    count_values = np.array(sorted([k for k in count_dict]))
    count_weights = np.array([count_dict[k] for k in count_values])

    mu_noise = weighted_mean(count_values[noise_mask], count_weights[noise_mask] * prob_noise[noise_mask])
    mu_signal = weighted_mean(count_values[signal_mask], count_weights[signal_mask] * prob_signal[signal_mask])
    alpha_noise = MLE_dispersion(count_values[noise_mask], count_weights[noise_mask] * prob_noise[noise_mask])
    alpha_signal = MLE_dispersion(count_values[signal_mask], count_weights[signal_mask] * prob_signal[signal_mask])
    frac_noise = (prob_noise * count_weights).sum() / (prob_matrix * count_weights).sum()
    return normalize_parameters(MixedModelParams(mu_noise, alpha_noise, mu_signal, alpha_signal, frac_noise))

def estimate_threshold(params, odds_ratio=20, signal_cdf_threshold = 0.995):
    """Estimate the separation threshold between signal and noise using the fits - look for the lowest count that
    gives an appropriate odds ratio in favor of signal
    """

    def get_ratio(k):
        p_signal = pmf(k, params.mu_signal, params.alpha_signal)
        signal = (1 - params.frac_noise) * p_signal
        p_noise = pmf(k, params.mu_noise, params.alpha_noise)
        noise = params.frac_noise * p_noise
        
        signal_captured = 1 - cdf(k, params.mu_signal, params.alpha_signal)
        return (signal / noise, signal_captured)

    for test_count in np.arange(int(params.mu_noise), int(params.mu_signal) + 1):
        ratio, signal_captured = get_ratio(test_count)
        if ratio >= odds_ratio and signal_captured <= signal_cdf_threshold:
            return test_count
    return test_count


def k_means_estimate(count_data):
    """
    Performs a very rough estimate of threshold separating noise and signal with 1-D kmeans.
    """
    centroids, labels = kmeans2(k=[np.log10(1), np.log10(count_data.max())], data=np.log10(count_data + 1), missing='raise', minit='matrix')
    noise_guess = np.argmin(centroids)
    max_noise_guess = count_data[labels == noise_guess].max()
    return max_noise_guess


def normalize_parameters(params):
    """Ensure that the larger of the two distributions is always the signal distribution."""
    if params.mu_noise > params.mu_signal:
        return MixedModelParams(params.mu_signal, params.alpha_signal,
                                params.mu_noise, params.alpha_noise,
                                1 - params.frac_noise)
    return params


def estimate_parameters(count_data, epsilon=1e-6, maxiter=250, threshold=None):
    """Uses an expectation-maximization method to estimate the mixed model parameters."""
    if threshold is None:
        threshold = k_means_estimate(count_data)
    
    count_dict = collections.Counter(count_data)
    
    count_values = np.array(sorted([k for k in count_dict]))

    prob_matrix = np.empty((2, len(count_values)), dtype=float)
    prob_matrix[0, :] = count_values <= threshold
    prob_matrix[1, :] = count_values > threshold
    params = weighted_parameter_estimates(count_dict, prob_matrix)
    
    last_params = None
    i = 0
    while i < maxiter and (last_params is None or
                           any([abs(last_params[j] - params[j]) / params[j] > epsilon for j in range(len(params))])):

        i += 1
        prob_matrix = calculate_probability_distribution(count_dict, params)
        last_params = params
        params = weighted_parameter_estimates(count_dict, prob_matrix)

    return params


def compute_cell_index(species_list, cell_barcodes):
    """Produces a cell index dictionary that is keyed on barcode and links cell barcodes to per-species indices
    """
    cell_index = {}
    for species in species_list:
        bc_list = cell_barcodes.get(species, {}).keys()
        bc_list.sort()
        for index, barcode in enumerate(bc_list):
            label = "{species}_cell_{index}".format(**locals())
            if barcode not in cell_index:
                cell_index[barcode] = label
            else:
                cell_index[barcode] = "_".join([cell_index[barcode], label])
    return cell_index

def choose_shift(count_data, min_shift=2, max_shift=12):
    """
    Helper function to choose the best shift to apply to the data before applying mixture model.
    This helps automatically pick a shift vs. just having a fixed
    """
    best_fit = None
    best_shift = None

    for count_shift in range(min_shift, max_shift, 1):
        print('Testing shift of %s...' % count_shift)
        count_data_shifted = count_data[count_data >= count_shift] - count_shift

        fitted_params = estimate_parameters(count_data_shifted, maxiter=50)
        current_fit = goodness_of_fit(count_data_shifted, fitted_params)
        
        if best_fit is None or current_fit > best_fit:
            best_fit = current_fit
            best_shift = count_shift

    return best_shift

if __name__ == '__main__':
    """
    The barcode count distribution is likely bimodal so modeled as two negative binomial
    distributions -- one is (mostly) cells and the other is background. This script finds
    good parameters for the two distributions. Then it finds a location between the two
    distributions: "where the minimum count that both yields an odds ratio (in favor of
    signal) of 20 or higher and removes at least 0.5% of the signal distribution as
    estimated from the signal distribution's CDF (we found that this second criterion
    prevents fits with thresholds that appeared to be too permissive otherwise)."
    This description is from
      A human cell atlas of fetal chromatin accessibility
      Domcke et al.
      Science 370 (2020)
      DOI: 10.1126/science.aba7612
    """
    parser = argparse.ArgumentParser('Script to call cells given the barcode count distribution.')
    parser.add_argument('--count_report', required=True, help='Count report as output by earlier steps in pipline.')
    parser.add_argument('called_cells_counts', help='Output of counts for all called cells.')
    parser.add_argument('cell_whitelist', help='Output of counts for all called cells.')
    parser.add_argument('--reads_in_peaks', action='store_true', help='Set this flag to call cells using only reads in peaks rather than all reads.')
    parser.add_argument('--reads_threshold', type=int, help='Count report as output by earlier steps in pipline.')
    parser.add_argument('--fit_metadata', help='JSON file with stats on fit and chosen shift.')
    args = parser.parse_args()

    counts = pd.read_csv(args.count_report, sep='\t')

    # Take MOLECULES (divide by two) either in or out of peaks
    if args.reads_in_peaks:
        counts['input_data'] = (counts.total_deduplicated_peaks / 2).astype(int)
    else:
        counts['input_data'] = (counts.total_deduplicated / 2).astype(int)

    # Do full EM
    parameters = {}
    parameters['noise_mean'] = None
    parameters['noise_dispersion'] = None
    parameters['signal_mean'] = None
    parameters['signal_dispersion'] = None
    parameters['fraction_noise'] = None
    parameters['cell_threshold'] = None
    parameters['cells_detected'] = None
    parameters['estimated_cells_present'] = None
    parameters['chosen_shift'] = None

    if args.reads_threshold:
        parameters['cell_threshold'] = args.reads_threshold
        parameters['cells_detected'] = counts.input_data[counts.input_data > args.reads_threshold].size
    else:
        # Search for the best shift in the data over reasonable set of values
        count_data = np.array(counts.input_data)
        
        print('Searching for best shift:')
        count_shift = choose_shift(count_data)
        print('Shift of %s chosen.' % count_shift)

        count_data = count_data[count_data >= count_shift] - count_shift

        print('Performing EM to estimate mixture model parameters...')
        fitted_params = estimate_parameters(count_data)

        print('Parameter estimation done.')

        # Threshold reported should be for reads not molecules (so multiply by 2)
        signal_molecules_threshold = estimate_threshold(fitted_params, 100) + count_shift
        signal_reads_threshold = 2 * signal_molecules_threshold
    
        print('Threshold based on mixture of NB (reads): {}'.format(signal_reads_threshold))
        parameters['noise_mean'] = float(fitted_params.mu_noise)
        parameters['signal_mean'] = float(fitted_params.mu_signal)
        parameters['noise_dispersion'] = float(fitted_params.alpha_noise)
        parameters['signal_dispersion'] = float(fitted_params.alpha_signal)
        parameters['fraction_noise'] = float(fitted_params.frac_noise)
        parameters['cell_threshold'] = int(signal_reads_threshold)
        parameters['goodness_of_fit'] = float(goodness_of_fit(count_data, fitted_params))
        called_cell_count = int(np.sum(count_data >= signal_reads_threshold))
        parameters['cells_detected'] = int(called_cell_count)
        parameters['estimated_cells_present'] = int((1 - fitted_params.frac_noise) * len(count_data))
        parameters['chosen_shift'] = int(count_shift)

    # Save model results/parameters to file if requested
    if args.fit_metadata:
        with open(args.fit_metadata, 'w') as f:
            f.write(json.dumps(parameters, indent=4))

    # Generate other output files
    cell_molecules_threshold = parameters['cell_threshold'] / 2
    under_depth_threshold = counts.input_data < cell_molecules_threshold
    counts = counts.loc[~under_depth_threshold]

    counts.drop(['input_data'], axis=1).to_csv(args.called_cells_counts, index=False, sep='\t')
    counts[['cell']].to_csv(args.cell_whitelist, index=False, header=False)


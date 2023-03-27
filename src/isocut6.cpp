/*
 * Copyright 2016-2017 Flatiron Institute, Simons Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "isocut5.h"
#include "isocut6.h"
#include "jisotonic5.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void isocut6(double* dipscore_out, double* cutpoint_out, bigint N, double* samples, isocut6_opts opts)
{
    double* samples_sorted = (double*)malloc(sizeof(double) * N);

    // sort the samples if needed
    if (opts.already_sorted)
        ns_isocut5::copy_samples(N, samples_sorted, samples);
    else
        jisotonic5_sort(N, samples_sorted, samples);
    double* X = samples_sorted;

    double* spacings = (double*)malloc(sizeof(double) * (N - 1));
    double* multiplicities = (double*)malloc(sizeof(double) * (N - 1));
    double* log_densities = (double*)malloc(sizeof(double) * (N - 1));
    for (bigint i = 0; i < N - 1; i++) {
        spacings[i] = X[i + 1] - X[i];
        multiplicities[i] = 1;
        if (spacings[i]) {
            log_densities[i] = log(multiplicities[i] / spacings[i]);
        }
        else {
            log_densities[i] = log(0.000000001); // I hope this is small enough
        }
    }

    double* log_densities_unimodal_fit = (double*)malloc(sizeof(double) * (N - 1));
    jisotonic5_updown(N - 1, log_densities_unimodal_fit, log_densities, multiplicities);

    double* densities_unimodal_fit_times_spacings = (double*)malloc(sizeof(double) * (N - 1));
    for (bigint i = 0; i < N - 1; i++)
        densities_unimodal_fit_times_spacings[i] = exp(log_densities_unimodal_fit[i]) * spacings[i];

    bigint critical_range_min, critical_range_max;
    bigint peak_index = ns_isocut5::find_max_index(N - 1, log_densities_unimodal_fit);

    *dipscore_out = ns_isocut5::compute_ks5(&critical_range_min, &critical_range_max, N - 1, multiplicities, densities_unimodal_fit_times_spacings, peak_index);

    bigint critical_range_length = critical_range_max - critical_range_min + 1;

    double* log_densities_resid = (double*)malloc(sizeof(double) * (N - 1));
    double* log_densities_resid_on_critical_range = (double*)malloc(sizeof(double) * (critical_range_length));
    double* log_densities_resid_fit_on_critical_range = (double*)malloc(sizeof(double) * (critical_range_length));
    double* weights_for_downup = (double*)malloc(sizeof(double) * (critical_range_length));

    for (bigint i = 0; i < N - 1; i++)
        log_densities_resid[i] = log_densities[i] - log_densities_unimodal_fit[i];

    for (bigint i = 0; i < critical_range_length; i++) {
        log_densities_resid_on_critical_range[i] = log_densities_resid[critical_range_min + i];
        weights_for_downup[i] = 1; // should we use spacings for the weights here?
    }
    jisotonic5_downup(critical_range_length, log_densities_resid_fit_on_critical_range, log_densities_resid_on_critical_range, weights_for_downup);

    bigint cutpoint_index = ns_isocut5::find_min_index(critical_range_length, log_densities_resid_fit_on_critical_range);
    *cutpoint_out = (X[critical_range_min + cutpoint_index] + X[critical_range_min + cutpoint_index + 1]) / 2;

    free(samples_sorted);
    free(log_densities_unimodal_fit);
    free(log_densities_resid);
    free(weights_for_downup);
    free(multiplicities);
    free(spacings);
    free(log_densities);
    free(densities_unimodal_fit_times_spacings);
    free(log_densities_resid_on_critical_range);
    free(log_densities_resid_fit_on_critical_range);
}
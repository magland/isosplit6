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
#ifndef ISOCUT5_H
#define ISOCUT5_H

#include <stdlib.h>
#include <cstdint>

typedef int64_t bigint;

struct isocut5_opts {
    bool already_sorted = false;
};

namespace ns_isocut5 {
void copy_samples(bigint N, double* out, double* in);
double sum(bigint N, double* X);
bigint find_min_index(bigint N, double* X);
bigint find_max_index(bigint N, double* X);
double compute_ks4(bigint N, double* counts1, double* counts2);
double compute_ks5(bigint* critical_range_min, bigint* critical_range_max, bigint N, double* counts1, double* counts2, bigint peak_index);
void debug_print_array(bigint N, double* X);
}

void isocut5(double* dipscore_out, double* cutpoint_out, bigint N, double* samples, isocut5_opts opts);

/*
 * MCWRAP [ dipscore[1,1], cutpoint[1,1] ] = isocut5_mex(samples[1,N])
 * SET_INPUT N = size(samples,2)
 * SOURCES isocut5.cpp jisotonic5.cpp
 * HEADERS isocut5.h jisotonic5.h
 */
void isocut5_mex(double* dipscore, double* cutpoint, int N, double* samples);

#endif // ISOCUT5_H

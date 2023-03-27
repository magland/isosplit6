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
#ifndef ISOSPLIT6_H
#define ISOSPLIT6_H

//#include "mlcommon.h"
#include "isocut6.h"

struct isosplit6_opts {
    double isocut_threshold = 2.0;
    int min_cluster_size = 10;
    int K_init = 200;
    bool refine_clusters = false;
    int max_iterations_per_pass = 500;
};

bool isosplit6(int* labels_out, bigint M, bigint N, double* X, isosplit6_opts opts);
 
#endif // ISOSPLIT6_H

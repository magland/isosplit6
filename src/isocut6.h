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
#ifndef ISOCUT6_H
#define ISOCUT6_H

#include <stdlib.h>
#include <cstdint>

typedef int64_t bigint;

struct isocut6_opts {
    bool already_sorted = false;
};

void isocut6(double* dipscore_out, double* cutpoint_out, bigint N, double* samples, isocut6_opts opts);

#endif // ISOCUT6_H

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

#ifndef jisotonic5_h
#define jisotonic5_h

#include <stdlib.h>
#include <cstdint>

typedef int64_t bigint;

void jisotonic5(bigint N, double* BB, double* MSE, double* AA, double* WW);
void jisotonic5_updown(bigint N, double* out, double* in, double* weights);
void jisotonic5_downup(bigint N, double* out, double* in, double* weights);
void jisotonic5_sort(bigint N, double* out, const double* in);

#endif

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
#ifndef ISOSPLIT5_H
#define ISOSPLIT5_H

//#include "mlcommon.h"
#include "isocut5.h"

#include <vector>

struct isosplit5_opts {
    double isocut_threshold = 1.0;
    int min_cluster_size = 10;
    int K_init = 200;
    bool refine_clusters = false;
    int max_iterations_per_pass = 500;
};

typedef std::vector<std::vector<bigint> > intarray2d;

namespace ns_isosplit5 {
struct kmeans_opts {
    bigint num_iterations = 0;
};

bigint compute_max(bigint N, int* labels);
bigint compute_max(bigint N, bigint* inds);
void kmeans_multistep(int* labels, bigint M, bigint N, double* X, bigint K1, bigint K2, bigint K3, kmeans_opts opts);
void kmeans_maxsize(int* labels, bigint M, bigint N, double* X, bigint maxsize, kmeans_opts opts);
void compare_clusters(double* dip_score, std::vector<bigint>* new_labels1, std::vector<bigint>* new_labels2, bigint M, bigint N1, bigint N2, double* X1, double* X2, double* centroid1, double* centroid2);
void compute_centroids(double* centroids, bigint M, bigint N, bigint Kmax, double* X, int* labels, std::vector<bigint>& clusters_to_compute_vec);
void compute_covmats(double* covmats, bigint M, bigint N, bigint Kmax, double* X, int* labels, double* centroids, std::vector<bigint>& clusters_to_compute_vec);
void get_pairs_to_compare(std::vector<bigint>* inds1, std::vector<bigint>* inds2, bigint M, bigint K, double* active_centroids, const intarray2d& active_comparisons_made);
void compare_pairs(std::vector<bigint>* clusters_changed, bigint* total_num_label_changes, bigint M, bigint N, double* X, int* labels, const std::vector<bigint>& inds1, const std::vector<bigint>& inds2, const isosplit5_opts& opts, double* centroids, double* covmats); //the labels are updated
bool matinv(bigint M, double* out, double* in);
void matvec(bigint M, bigint N, double* out, double* mat, double* vec);
void extract_subarray(double* X_sub, bigint M, double* X, const std::vector<bigint>& inds);
}

namespace smi {
bool get_inverse_via_lu_decomposition(int M, double* out, double* in);
};

struct parcelate2_opts {
    bool final_reassign = false; //not yet implemented
};
bool parcelate2(int* labels, bigint M, bigint N, double* X, bigint target_parcel_size, bigint target_num_parcels, const parcelate2_opts& p2opts);
void alloc(intarray2d& X, bigint N1, bigint N2);

bool isosplit5(int* labels_out, bigint M, bigint N, double* X, isosplit5_opts opts);

/*
 * MCWRAP [ labels_out[1,N] ] = isosplit5_mex(X[M,N])
 * SET_INPUT M = size(X,1)
 * SET_INPUT N = size(X,2)
 * SOURCES isosplit5.cpp isocut5.cpp jisotonic5.cpp
 * HEADERS isosplit5.h isocut5.h jisotonic5.h
 */
void isosplit5_mex(double* labels_out, int M, int N, double* X);

#endif // ISOSPLIT5_H

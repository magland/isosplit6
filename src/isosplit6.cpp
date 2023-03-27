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
#include "isosplit5.h"
#include "isosplit6.h"
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <math.h>
#include "isocut6.h"

namespace ns_isosplit6 {
void compare_pairs(std::vector<bigint>* clusters_changed, bigint* total_num_label_changes, bigint M, bigint N, double* X, int* labels, const std::vector<bigint>& inds1, const std::vector<bigint>& inds2, const isosplit6_opts& opts, double* centroids, double* covmats); //the labels are updated
}

bool isosplit6(int* labels, bigint M, bigint N, double* X, isosplit6_opts opts)
{
    // compute the initial clusters
    bigint target_parcel_size = opts.min_cluster_size;
    bigint target_num_parcels = opts.K_init;
    // !! important not to do a final reassign because then the shapes will not be conducive to isosplit iterations -- hexagons are not good for isosplit!
    parcelate2_opts p2opts;
    p2opts.final_reassign = false;
    if (!parcelate2(labels, M, N, X, target_parcel_size, target_num_parcels, p2opts)) {
        for (bigint i=0; i<N; i++) {
            labels[i]=-1;
        }
        printf("Failure in parcelate2.\n");
        return false;
    }
    int Kmax = ns_isosplit5::compute_max(N, labels);

    double* centroids = (double*)malloc(sizeof(double) * M * Kmax);
    double* covmats = (double*)malloc(sizeof(double) * M * M * Kmax);
    std::vector<bigint> clusters_to_compute_vec;
    for (bigint k = 0; k < Kmax; k++)
        clusters_to_compute_vec.push_back(1);
    ns_isosplit5::compute_centroids(centroids, M, N, Kmax, X, labels, clusters_to_compute_vec);
    ns_isosplit5::compute_covmats(covmats, M, N, Kmax, X, labels, centroids, clusters_to_compute_vec);

    // The active labels are those that are still being used -- for now, everything is active
    std::vector<int> active_labels_vec(Kmax, 1);
    std::vector<int> active_labels;
    for (bigint i = 0; i < Kmax; i++)
        active_labels.push_back(i + 1);

    // Repeat while something has been merged in the pass
    bool final_pass = false; // plus we do one final pass at the end
    intarray2d comparisons_made; // Keep a matrix of comparisons that have been made in this pass
    alloc(comparisons_made, Kmax, Kmax);
    for (bigint i1 = 0; i1 < Kmax; i1++)
        for (bigint i2 = 0; i2 < Kmax; i2++)
            comparisons_made[i1][i2] = 0;
    while (true) { //passes
        bool something_merged = false; //Keep track of whether something has merged in this pass. If not, do a final pass.
        std::vector<bigint> clusters_changed_vec_in_pass(Kmax); //Keep track of the clusters that have changed in this pass so that we can update the comparisons_made matrix at the end
        for (bigint i = 0; i < Kmax; i++)
            clusters_changed_vec_in_pass[i] = 0;
        bigint iteration_number = 0;
        while (true) { //iterations

            std::vector<bigint> clusters_changed_vec_in_iteration(Kmax); //Keep track of the clusters that have changed in this iteration so that we can update centroids and covmats
            for (bigint i = 0; i < Kmax; i++)
                clusters_changed_vec_in_iteration[i] = 0;

            iteration_number++;
            if (iteration_number > opts.max_iterations_per_pass) {
                printf("Warning: max iterations per pass exceeded.\n");
                break;
            }

            if (active_labels.size() > 0) {
                // Create an array of active centroids and comparisons made, for determining the pairs to compare
                double* active_centroids = (double*)malloc(sizeof(double) * M * active_labels.size());
                for (bigint i = 0; i < (bigint)active_labels.size(); i++) {
                    for (bigint m = 0; m < M; m++) {
                        active_centroids[m + M * i] = centroids[m + M * (active_labels[i] - 1)];
                    }
                }
                intarray2d active_comparisons_made;
                alloc(active_comparisons_made, active_labels.size(), active_labels.size());
                for (bigint i1 = 0; i1 < (bigint)active_labels.size(); i1++) {
                    for (bigint i2 = 0; i2 < (bigint)active_labels.size(); i2++) {
                        active_comparisons_made[i1][i2] = comparisons_made[active_labels[i1] - 1][active_labels[i2] - 1];
                    }
                }

                // Find the pairs to compare on this iteration
                // These will be closest pairs of active clusters that have not yet
                // been compared in this pass
                std::vector<bigint> inds1, inds2;
                ns_isosplit5::get_pairs_to_compare(&inds1, &inds2, M, active_labels.size(), active_centroids, active_comparisons_made);
                std::vector<bigint> inds1b, inds2b; //remap the clusters to the original labeling
                for (bigint i = 0; i < (bigint)inds1.size(); i++) {
                    inds1b.push_back(active_labels[inds1[i] - 1]);
                    inds2b.push_back(active_labels[inds2[i] - 1]);
                }

                // If we didn't find any, break from this iteration
                if (inds1b.size() == 0) {
                    break;
                }

                // Actually compare the pairs -- in principle this operation could be parallelized
                std::vector<bigint> clusters_changed;
                bigint total_num_label_changes = 0;
                ns_isosplit6::compare_pairs(&clusters_changed, &total_num_label_changes, M, N, X, labels, inds1b, inds2b, opts, centroids, covmats); //the labels are updated
                for (bigint i = 0; i < (bigint)clusters_changed.size(); i++) {
                    clusters_changed_vec_in_pass[clusters_changed[i] - 1] = 1;
                    clusters_changed_vec_in_iteration[clusters_changed[i] - 1] = 1;
                }

                // Update which comparisons have been made
                for (bigint j = 0; j < (bigint)inds1b.size(); j++) {
                    comparisons_made[inds1b[j] - 1][inds2b[j] - 1] = 1;
                    comparisons_made[inds2b[j] - 1][inds1b[j] - 1] = 1;
                }

                // Recompute the centers for those that have changed in this iteration
                ns_isosplit5::compute_centroids(centroids, M, N, Kmax, X, labels, clusters_changed_vec_in_iteration);
                ns_isosplit5::compute_covmats(covmats, M, N, Kmax, X, labels, centroids, clusters_changed_vec_in_iteration);

                // For diagnostics
                //printf ("total num label changes = %d\n",total_num_label_changes);

                // Determine whether something has merged and update the active labels
                for (bigint i = 0; i < Kmax; i++)
                    active_labels_vec[i] = 0;
                for (bigint i = 0; i < N; i++)
                    active_labels_vec[labels[i] - 1] = 1;
                std::vector<int> new_active_labels;
                for (bigint i = 0; i < Kmax; i++)
                    if (active_labels_vec[i])
                        new_active_labels.push_back(i + 1);
                if (new_active_labels.size() < active_labels.size())
                    something_merged = true;
                active_labels = new_active_labels;

                free(active_centroids);
            }
        }

        // zero out the comparisons made matrix only for those that have changed in this pass
        for (bigint i = 0; i < Kmax; i++) {
            if (clusters_changed_vec_in_pass[i]) {
                for (bigint j = 0; j < Kmax; j++) {
                    comparisons_made[i][j] = 0;
                    comparisons_made[j][i] = 0;
                }
            }
        }

        if (something_merged)
            final_pass = false;
        if (final_pass)
            break; // This was the final pass and nothing has merged
        if (!something_merged)
            final_pass = true; // If we are done, do one last pass for final redistributes
    }

    // We should remap the labels to occupy the first natural numbers
    std::vector<bigint> labels_map(Kmax, 0);
    for (bigint i = 0; i < (bigint)active_labels.size(); i++) {
        labels_map[active_labels[i] - 1] = i + 1;
    }
    for (bigint i = 0; i < N; i++) {
        labels[i] = labels_map[labels[i] - 1];
    }

    // If the user wants to refine the clusters, then we repeat isosplit on each
    // of the new clusters, recursively. Unless we only found only one cluster.
    bigint K = ns_isosplit5::compute_max(N, labels);

    if ((opts.refine_clusters) && (K > 1)) {
        int* labels_split = (int*)malloc(sizeof(int) * N);
        isosplit6_opts opts2 = opts;
        opts2.refine_clusters = true; // Maybe we should provide an option on whether to do recursive refinement
        bigint k_offset = 0;
        for (bigint k = 1; k <= K; k++) {
            std::vector<bigint> inds_k;
            for (bigint i = 0; i < N; i++)
                if (labels[i] == k)
                    inds_k.push_back(i);
            if (inds_k.size() > 0) {
                double* X_k = (double*)malloc(sizeof(double) * M * inds_k.size()); //Warning: this may cause memory problems -- especially for recursive case
                int* labels_k = (int*)malloc(sizeof(int) * inds_k.size());
                for (bigint i = 0; i < (bigint)inds_k.size(); i++) {
                    for (bigint m = 0; m < M; m++) {
                        X_k[m + M * i] = X[m + M * inds_k[i]];
                    }
                }
                isosplit6(labels_k, M, inds_k.size(), X_k, opts2);
                for (bigint i = 0; i < (bigint)inds_k.size(); i++) {
                    labels_split[inds_k[i]] = k_offset + labels_k[i];
                }
                k_offset += ns_isosplit5::compute_max(inds_k.size(), labels_k);
                free(labels_k);
                free(X_k);
            }
        }
        for (bigint i = 0; i < N; i++)
            labels[i] = labels_split[i];
        free(labels_split);
    }

    free(centroids);
    free(covmats);

    return true;
}

namespace ns_isosplit6 {
bool merge_test(std::vector<bigint>* L12, bigint M, bigint N1, bigint N2, double* X1, double* X2, const isosplit6_opts& opts, double* centroid1, double* centroid2, double* covmat1, double* covmat2)
{
    L12->resize(N1 + N2);
    for (bigint i = 0; i < N1 + N2; i++)
        (*L12)[i] = 1;
    if ((N1 == 0) || (N2 == 0)) {
        printf("Error in merge test: N1 or N2 is zero.\n");
        return true;
    }

    //std::vector<double> centroid1 = compute_centroid(M, N1, X1);
    //std::vector<double> centroid2 = compute_centroid(M, N2, X2);

    std::vector<double> V(M);
    for (bigint m = 0; m < M; m++) {
        V[m] = centroid2[m] - centroid1[m];
    }

    std::vector<double> avg_covmat;
    avg_covmat.resize(M * M);
    for (bigint rr = 0; rr < M * M; rr++) {
        avg_covmat[rr] = (covmat1[rr] + covmat2[rr]) / 2;
    }
    std::vector<double> inv_avg_covmat;
    inv_avg_covmat.resize(M * M);
    if (!ns_isosplit5::matinv(M, inv_avg_covmat.data(), avg_covmat.data())) {
        fprintf(stderr, "Unable to invert matrix. This may be due to the fact that you have duplicate events. Contact Jeremy if this is not the case, or if you would prefer the program to continue in this case. Aborting.\n");
        abort();
        return false;
    }

    std::vector<double> V2(M);
    ns_isosplit5::matvec(M, M, V2.data(), inv_avg_covmat.data(), V.data());
    //ns_isosplit5::matvec(M,M,V.data(),inv_avg_covmat.data(),V2.data());
    for (bigint i = 0; i < M; i++)
        V[i] = V2[i];

    /*
    printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");
    print_matrix(M,M,avg_covmat.data());
    printf("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB\n");
    print_matrix(M,M,inv_avg_covmat.data());
    printf("\n\n");
    */

    double sumsqr = 0;
    for (bigint m = 0; m < M; m++) {
        sumsqr += V[m] * V[m];
    }
    if (sumsqr) {
        for (bigint m = 0; m < M; m++)
            V[m] /= sqrt(sumsqr);
    }

    std::vector<double> projection1(N1), projection2(N2), projection12(N1 + N2);
    for (bigint i = 0; i < N1; i++) {
        double tmp = 0;
        for (bigint m = 0; m < M; m++)
            tmp += V[m] * X1[m + i * M];
        projection1[i] = tmp;
        projection12[i] = tmp;
    }
    for (bigint i = 0; i < N2; i++) {
        double tmp = 0;
        for (bigint m = 0; m < M; m++)
            tmp += V[m] * X2[m + i * M];
        projection2[i] = tmp;
        projection12[N1 + i] = tmp;
    }

    bool do_merge;
    isocut6_opts oo;
    oo.already_sorted = false;
    double dipscore, cutpoint;
    isocut6(&dipscore, &cutpoint, N1 + N2, projection12.data(), oo);

    if (dipscore < opts.isocut_threshold) {
        do_merge = true;
    }
    else {
        do_merge = false;
    }
    for (bigint i = 0; i < N1 + N2; i++) {
        if (projection12[i] < cutpoint)
            (*L12)[i] = 1;
        else
            (*L12)[i] = 2;
    }

    return do_merge;
}

void compare_pairs(std::vector<bigint>* clusters_changed, bigint* total_num_label_changes, bigint M, bigint N, double* X, int* labels, const std::vector<bigint>& k1s, const std::vector<bigint>& k2s, const isosplit6_opts& opts, double* centroids, double* covmats)
{
    bigint Kmax = ns_isosplit5::compute_max(N, labels);
    std::vector<bigint> clusters_changed_vec(Kmax);
    for (bigint i = 0; i < Kmax; i++)
        clusters_changed_vec[i] = 0;
    int* new_labels = (int*)malloc(sizeof(bigint) * N);
    *total_num_label_changes = 0;
    for (bigint i = 0; i < N; i++)
        new_labels[i] = labels[i];
    for (bigint i1 = 0; i1 < (bigint)k1s.size(); i1++) {
        int k1 = k1s[i1];
        int k2 = k2s[i1];
        std::vector<bigint> inds1, inds2;
        for (bigint i = 0; i < N; i++) {
            if (labels[i] == k1)
                inds1.push_back(i);
            if (labels[i] == k2)
                inds2.push_back(i);
        }
        if ((inds1.size() > 0) && (inds2.size() > 0)) {
            std::vector<bigint> inds12;
            inds12.insert(inds12.end(), inds1.begin(), inds1.end());
            inds12.insert(inds12.end(), inds2.begin(), inds2.end());
            std::vector<bigint> L12_old(inds12.size());
            for (bigint i = 0; i < (bigint)inds1.size(); i++)
                L12_old[i] = 1;
            for (bigint i = 0; i < (bigint)inds2.size(); i++)
                L12_old[inds1.size() + i] = 2;
            std::vector<bigint> L12(inds12.size());

            bool do_merge;
            if (((bigint)inds1.size() < opts.min_cluster_size) || ((bigint)inds2.size() < opts.min_cluster_size)) {
                do_merge = true;
            }
            else {
                double* X1 = (double*)malloc(sizeof(double) * M * inds1.size());
                double* X2 = (double*)malloc(sizeof(double) * M * inds2.size());
                ns_isosplit5::extract_subarray(X1, M, X, inds1);
                ns_isosplit5::extract_subarray(X2, M, X, inds2);
                do_merge = merge_test(&L12, M, inds1.size(), inds2.size(), X1, X2, opts, &centroids[(k1 - 1) * M], &centroids[(k2 - 1) * M], &covmats[(k1 - 1) * M * M], &covmats[(k2 - 1) * M * M]);
                free(X1);
                free(X2);
            }
            if (do_merge) {
                for (bigint i = 0; i < (bigint)inds2.size(); i++) {
                    new_labels[inds2[i]] = k1;
                }
                *total_num_label_changes += inds2.size();
                clusters_changed_vec[k1 - 1] = 1;
                clusters_changed_vec[k2 - 1] = 1;
            }
            else {
                //redistribute
                bool something_was_redistributed = false;
                for (bigint i = 0; i < (bigint)inds1.size(); i++) {
                    if (L12[i] == 2) {
                        new_labels[inds1[i]] = k2;
                        (*total_num_label_changes)++;
                        something_was_redistributed = true;
                    }
                }
                for (bigint i = 0; i < (bigint)inds2.size(); i++) {
                    if (L12[inds1.size() + i] == 1) {
                        new_labels[inds2[i]] = k1;
                        (*total_num_label_changes)++;
                        something_was_redistributed = true;
                    }
                }
                if (something_was_redistributed) {
                    clusters_changed_vec[k1 - 1] = 1;
                    clusters_changed_vec[k2 - 1] = 1;
                }
            }
        }
    }
    clusters_changed->clear();
    for (int k = 0; k < Kmax; k++)
        if (clusters_changed_vec[k])
            clusters_changed->push_back(k + 1);
    for (bigint i = 0; i < N; i++)
        labels[i] = new_labels[i];
    free(new_labels);
}
}
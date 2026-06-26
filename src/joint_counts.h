#ifndef ACOR_JOINT_COUNTS_H
#define ACOR_JOINT_COUNTS_H

#include <vector>

// ============================================================================
// Joint-count statistics over a compressed bivariate sample (Xc, Yc).
// Shared by H_bar_vec_v2_cpp and kernel_agc_v2_cpp.
// ============================================================================

struct JointCounts {
  std::vector<int> count_both_less;    // #{j : x_j <  x_i, y_j <  y_i}
  std::vector<int> count_x_less_y_eq;  // #{j : x_j <  x_i, y_j == y_i}
  std::vector<int> count_x_eq_y_less;  // #{j : x_j == x_i, y_j <  y_i, j!=i}
  std::vector<int> count_both_eq;      // size of i's (X,Y) tie group (incl. i)
};

// Sort-based coordinate compression: maps doubles to consecutive ints 1..M.
std::vector<int> compress(const std::vector<double>& vals, int& M);

JointCounts compute_joint_counts(
    const std::vector<int>& Xc, int Mx,
    const std::vector<int>& Yc, int My, int n);

#endif  // ACOR_JOINT_COUNTS_H

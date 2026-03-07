#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <map>
using namespace Rcpp;

// ============================================================================
// Fenwick Tree (Binary Indexed Tree)
// Supports point updates and prefix sum queries in O(log n)
// ============================================================================

class FenwickTree {
  std::vector<int> tree;
  int n;
public:
  FenwickTree(int n) : tree(n + 1, 0), n(n) {}
  
  void update(int i, int delta = 1) {
    for (; i <= n; i += i & (-i))
      tree[i] += delta;
  }
  
  int query(int i) {
    int s = 0;
    for (; i > 0; i -= i & (-i))
      s += tree[i];
    return s;
  }
};

// ============================================================================
// Coordinate compression
// Maps arbitrary doubles to consecutive integers 1..M
// ============================================================================

std::vector<int> compress(const std::vector<double>& vals, int& M) {
  int n = vals.size();
  std::vector<double> sorted_unique = vals;
  std::sort(sorted_unique.begin(), sorted_unique.end());
  sorted_unique.erase(
    std::unique(sorted_unique.begin(), sorted_unique.end()),
    sorted_unique.end()
  );
  M = sorted_unique.size();
  std::vector<int> compressed(n);
  for (int i = 0; i < n; i++) {
    compressed[i] = std::lower_bound(
      sorted_unique.begin(), sorted_unique.end(), vals[i]
    ) - sorted_unique.begin() + 1;  // 1-indexed
  }
  return compressed;
}

// ============================================================================
// Group-by helper
// Returns groups as vector of vectors of indices, sorted by key (1..M)
// ============================================================================

std::vector<std::vector<int>> group_by_sorted(
    const std::vector<int>& keys, int M) {
  std::vector<std::vector<int>> groups(M);
  for (int i = 0; i < (int)keys.size(); i++) {
    groups[keys[i] - 1].push_back(i);
  }
  return groups;
}

// ============================================================================
// Shared Fenwick core: count_both_less, count_x_less_y_eq,
//                      count_x_eq_y_less, count_both_eq
// Used by both H_bar_vec_v2 and kernel_agc_v2
// ============================================================================

struct JointCounts {
  std::vector<int> count_both_less;
  std::vector<int> count_x_less_y_eq;
  std::vector<int> count_x_eq_y_less;
  std::vector<int> count_both_eq;
};

JointCounts compute_joint_counts(
    const std::vector<int>& Xc, int Mx,
    const std::vector<int>& Yc, int My, int n) {
  
  std::vector<std::vector<int>> X_groups = group_by_sorted(Xc, Mx);
  
  JointCounts out;
  out.count_both_less.assign(n, 0);
  out.count_x_less_y_eq.assign(n, 0);
  out.count_x_eq_y_less.assign(n, 0);
  out.count_both_eq.assign(n, 0);
  
  // Pass 1: Fenwick tree for count_both_less and count_x_less_y_eq
  FenwickTree tree(My);
  for (int g = 0; g < Mx; g++) {
    for (int idx : X_groups[g]) {
      int yc = Yc[idx];
      if (yc > 1)
        out.count_both_less[idx] = tree.query(yc - 1);
      out.count_x_less_y_eq[idx] = tree.query(yc) - tree.query(yc - 1);
    }
    for (int idx : X_groups[g]) {
      tree.update(Yc[idx]);
    }
  }
  
  // Pass 2: within-group counts for count_x_eq_y_less and count_both_eq
  for (int g = 0; g < Mx; g++) {
    std::vector<int>& indices = X_groups[g];
    int gsize = indices.size();
    if (gsize == 0) continue;
    
    // Sort within group by Y compressed rank
    std::vector<int> sorted_idx(indices);
    std::sort(sorted_idx.begin(), sorted_idx.end(),
              [&](int a, int b) { return Yc[a] < Yc[b]; });
    
    // count_x_eq_y_less: how many elements in same X-group have smaller Y
    int cumcount = 0;
    int prev_yr = -1;
    int same_count = 0;
    for (int k = 0; k < gsize; k++) {
      int idx = sorted_idx[k];
      int yr = Yc[idx];
      if (yr == prev_yr) {
        same_count++;
      } else {
        cumcount += same_count;
        same_count = 1;
        prev_yr = yr;
      }
      out.count_x_eq_y_less[idx] = cumcount;
    }
    
    // count_both_eq: how many elements share same (X,Y)
    std::map<int, int> pair_counts;
    for (int idx : indices) {
      pair_counts[Yc[idx]]++;
    }
    for (int idx : indices) {
      out.count_both_eq[idx] = pair_counts[Yc[idx]];
    }
  }
  
  return out;
}


// ============================================================================
// Shared concordant/discordant pair counting via Fenwick tree -- O(n log n)
//
// Counts pairs (i,j) with i<j where X[i]<X[j] & Y[i]<Y[j] (concordant)
// and X[i]<X[j] & Y[i]>Y[j] (discordant). X-tied pairs are excluded.
//
// Returns: {concordant, discordant}
// ============================================================================

struct PairCounts {
  double concordant;
  double discordant;
};

PairCounts count_concordant_discordant(NumericVector X, NumericVector Y) {
  int n = X.size();
  
  // Compress Y to 1..M
  std::map<double, int> y_map;
  for (int i = 0; i < n; i++) y_map[Y[i]];
  int rank = 1;
  for (auto& p : y_map) p.second = rank++;
  int M = y_map.size();
  
  std::vector<int> Yc(n);
  for (int i = 0; i < n; i++) Yc[i] = y_map[Y[i]];
  
  // Sort indices by (X, Y)
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
    if (X[a] != X[b]) return X[a] < X[b];
    return Y[a] < Y[b];
  });
  
  // Count concordant/discordant via Fenwick tree, processing by X-groups
  FenwickTree tree(M);
  double concordant = 0;
  double discordant = 0;
  
  int i = 0;
  while (i < n) {
    int j = i;
    while (j < n && X[ord[j]] == X[ord[i]]) j++;
    
    // Query before inserting (so X-tied pairs are not counted)
    for (int k = i; k < j; k++) {
      int yc = Yc[ord[k]];
      concordant += tree.query(yc - 1);
      discordant += tree.query(M) - tree.query(yc);
    }
    
    // Insert all elements of this X-group
    for (int k = i; k < j; k++)
      tree.update(Yc[ord[k]]);
    
    i = j;
  }
  
  return {concordant, discordant};
}


// ============================================================================
// Compute proportion of tied pairs from frequency counts -- O(n)
//
// For a vector with frequency table {n_1, n_2, ...}, the proportion of
// tied pairs is sum(n_k * (n_k - 1)) / (n * (n - 1)).
// This equals 1 - kendall_tau_a(v, v), but without the O(n log n) sort.
// ============================================================================

double pair_tie_proportion(NumericVector V) {
  int n = V.size();
  std::map<double, int> freq;
  for (int i = 0; i < n; i++) freq[V[i]]++;
  
  double n_tied_pairs = 0;
  for (auto& p : freq) {
    double nk = p.second;
    n_tied_pairs += nk * (nk - 1);
  }
  return n_tied_pairs / ((double)n * (n - 1));
}


// ============================================================================
// kendall_tau_sign_cpp
//
// Computes the asymmetric Kendall tau statistic using a Fenwick tree.
// Returns: list(tau, expectation)
// ============================================================================

// [[Rcpp::export]]
List kendall_tau_sign_cpp(NumericVector X, NumericVector Y) {
  int n = X.size();
  double num_pairs = (double)n * (n - 1) / 2.0;
  
  PairCounts pc = count_concordant_discordant(X, Y);
  
  double sum_sign = pc.concordant - pc.discordant;
  double expectation = sum_sign / num_pairs;
  double p_tie_y = pair_tie_proportion(Y);
  double scale_factor = 1.0 - p_tie_y;
  double tau = (scale_factor > 1e-10) ? expectation / scale_factor : 0.0;
  
  return List::create(
    Named("tau") = tau,
    Named("expectation") = expectation
  );
}


// ============================================================================
// kendall_tau_a
//
// Kendall's tau-a (no tie correction): (C - D) / n(n-1)/2
// ============================================================================

// [[Rcpp::export]]
double kendall_tau_a(NumericVector X, NumericVector Y) {
  int n = X.size();
  double num_pairs = (double)n * (n - 1) / 2.0;
  PairCounts pc = count_concordant_discordant(X, Y);
  return (pc.concordant - pc.discordant) / num_pairs;
}


// ============================================================================
// kendall_tau_b
//
// Kendall's tau-b (pair-based tie correction).
// Uses O(n) tie proportion instead of calling kendall_tau_a(V, V).
// ============================================================================

// [[Rcpp::export]]
double kendall_tau_b(NumericVector X, NumericVector Y) {
  double tau_a_xy = kendall_tau_a(X, Y);
  
  // tau_a(V, V) = 1 - pair_tie_proportion(V), so:
  //   sqrt(tau_a(X,X) * tau_a(Y,Y))
  // = sqrt((1 - p_x) * (1 - p_y))
  double p_x = pair_tie_proportion(X);
  double p_y = pair_tie_proportion(Y);
  double denom = std::sqrt((1.0 - p_x) * (1.0 - p_y));
  
  return (denom > 1e-10) ? tau_a_xy / denom : 0.0;
}


// ============================================================================
// goodman_kruskal_gamma
//
// Gamma = (C - D) / (C + D)
// ============================================================================

// [[Rcpp::export]]
double goodman_kruskal_gamma(NumericVector X, NumericVector Y) {
  PairCounts pc = count_concordant_discordant(X, Y);
  double denom = pc.concordant + pc.discordant;
  return (denom > 1e-10) ? (pc.concordant - pc.discordant) / denom : 0.0;
}


// ============================================================================
// kendall_tau_b_mod
//
// Modified tau-b using triple-based tie correction.
// ============================================================================

// [[Rcpp::export]]
double kendall_tau_b_mod(NumericVector X, NumericVector Y) {
  double tau_a_val = kendall_tau_a(X, Y);
  int n = X.size();
  
  // Count triple ties for X
  std::map<double, int> freq_x;
  for (int i = 0; i < n; i++) freq_x[X[i]]++;
  double triple_ties_x = 0;
  for (auto& p : freq_x) {
    double nk = p.second;
    triple_ties_x += nk * (nk - 1) * (nk - 2);
  }
  triple_ties_x /= (double)n * (n - 1) * (n - 2);
  
  // Count triple ties for Y
  std::map<double, int> freq_y;
  for (int i = 0; i < n; i++) freq_y[Y[i]]++;
  double triple_ties_y = 0;
  for (auto& p : freq_y) {
    double nk = p.second;
    triple_ties_y += nk * (nk - 1) * (nk - 2);
  }
  triple_ties_y /= (double)n * (n - 1) * (n - 2);
  
  double denom = std::sqrt((1.0 - triple_ties_x) * (1.0 - triple_ties_y));
  return (denom > 1e-10) ? tau_a_val / denom : 0.0;
}


// ============================================================================
// H_bar_vec_v2_cpp
//
// Computes H_bar for all observations using Fenwick tree -- O(n log n)
// Returns: numeric vector of length n
// ============================================================================

// [[Rcpp::export]]
NumericVector H_bar_vec_v2_cpp(NumericVector X, NumericVector Y) {
  int n = X.size();
  
  std::vector<double> Xv(X.begin(), X.end());
  std::vector<double> Yv(Y.begin(), Y.end());
  int Mx, My;
  std::vector<int> Xc = compress(Xv, Mx);
  std::vector<int> Yc = compress(Yv, My);
  
  JointCounts jc = compute_joint_counts(Xc, Mx, Yc, My, n);
  
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = (jc.count_both_less[i] +
      0.5 * jc.count_x_eq_y_less[i] +
      0.5 * jc.count_x_less_y_eq[i] +
      0.25 * jc.count_both_eq[i]) / n;
  }
  return result;
}


// ============================================================================
// kernel_agc_v2_cpp
//
// Computes the AGC kernel influence function values using Fenwick tree.
// Inputs: x_rank, y_rank (average ranks), rho (Spearman-style correlation)
// Returns: numeric vector k_p of length N
// ============================================================================

// [[Rcpp::export]]
NumericVector kernel_agc_v2_cpp(NumericVector x_rank, NumericVector y_rank,
                                double rho) {
  int N = x_rank.size();
  
  // Grade functions
  NumericVector G_x(N), G_y(N);
  for (int i = 0; i < N; i++) {
    G_x[i] = (x_rank[i] - 0.5) / N;
    G_y[i] = (y_rank[i] - 0.5) / N;
  }
  
  // Compress ranks
  std::vector<double> xr(x_rank.begin(), x_rank.end());
  std::vector<double> yr(y_rank.begin(), y_rank.end());
  int Mx, My;
  std::vector<int> Xc = compress(xr, Mx);
  std::vector<int> Yc = compress(yr, My);
  
  // Joint counts via shared core
  JointCounts jc = compute_joint_counts(Xc, Mx, Yc, My, N);
  
  // H_bar(x_i, y_i) for all i
  std::vector<double> Hbar(N);
  for (int i = 0; i < N; i++) {
    Hbar[i] = (jc.count_both_less[i] +
      0.5 * jc.count_x_eq_y_less[i] +
      0.5 * jc.count_x_less_y_eq[i] +
      0.25 * jc.count_both_eq[i]) / N;
  }
  
  // Frequency tables for Y
  std::vector<int> n_y(My, 0);
  for (int i = 0; i < N; i++) n_y[Yc[i] - 1]++;
  std::vector<int> cum_n_y(My + 1, 0);
  for (int m = 0; m < My; m++) cum_n_y[m + 1] = cum_n_y[m] + n_y[m];
  
  std::vector<int> n_y_eq(N), n_y_lt(N), n_y_gt(N);
  for (int i = 0; i < N; i++) {
    n_y_eq[i] = n_y[Yc[i] - 1];
    n_y_lt[i] = cum_n_y[Yc[i] - 1];
    n_y_gt[i] = N - n_y_lt[i] - n_y_eq[i];
  }
  
  // Frequency tables for X
  std::vector<int> n_x(Mx, 0);
  for (int i = 0; i < N; i++) n_x[Xc[i] - 1]++;
  std::vector<int> cum_n_x(Mx + 1, 0);
  for (int m = 0; m < Mx; m++) cum_n_x[m + 1] = cum_n_x[m] + n_x[m];
  
  std::vector<int> n_x_eq(N), n_x_lt(N), n_x_gt(N);
  for (int i = 0; i < N; i++) {
    n_x_eq[i] = n_x[Xc[i] - 1];
    n_x_lt[i] = cum_n_x[Xc[i] - 1];
    n_x_gt[i] = N - n_x_lt[i] - n_x_eq[i];
  }
  
  // Weights
  std::vector<double> w_full(N), w_half(N);
  for (int i = 0; i < N; i++) {
    w_full[i] = n_y_gt[i] + 0.5 * n_y_eq[i];
    w_half[i] = 0.5 * n_y_gt[i] + 0.25 * n_y_eq[i];
  }
  
  // Group processing for g_2 (S_i)
  std::vector<std::vector<int>> X_groups = group_by_sorted(Xc, Mx);
  
  // w_half group sums
  std::vector<double> w_half_group_sum(Mx, 0.0);
  for (int g = 0; g < Mx; g++)
    for (int idx : X_groups[g])
      w_half_group_sum[g] += w_half[idx];
  
  std::vector<double> sum_xeq_whalf(N);
  for (int i = 0; i < N; i++)
    sum_xeq_whalf[i] = w_half_group_sum[Xc[i] - 1];
  
  // Prefix sums for w_full (x_j < x_i)
  std::vector<double> sum_xlt_wfull(N, 0.0);
  double running_sum = 0.0;
  for (int g = 0; g < Mx; g++) {
    for (int idx : X_groups[g])
      sum_xlt_wfull[idx] = running_sum;
    for (int idx : X_groups[g])
      running_sum += w_full[idx];
  }
  
  // S_i = (1/N) * (sum_xlt_wfull + sum_xeq_whalf)
  std::vector<double> S_i(N);
  for (int i = 0; i < N; i++)
    S_i[i] = (1.0 / N) * (sum_xlt_wfull[i] + sum_xeq_whalf[i]);
  
  // T_m computation via prefix sums over Y groups
  std::vector<double> sum_nxgt_by_ygroup(My, 0.0);
  std::vector<double> sum_nxeq_by_ygroup(My, 0.0);
  for (int i = 0; i < N; i++) {
    sum_nxgt_by_ygroup[Yc[i] - 1] += n_x_gt[i];
    sum_nxeq_by_ygroup[Yc[i] - 1] += n_x_eq[i];
  }
  
  std::vector<double> cum_nxgt(My + 1, 0.0);
  std::vector<double> cum_nxeq(My + 1, 0.0);
  for (int m = 0; m < My; m++) {
    cum_nxgt[m + 1] = cum_nxgt[m] + sum_nxgt_by_ygroup[m];
    cum_nxeq[m + 1] = cum_nxeq[m] + sum_nxeq_by_ygroup[m];
  }
  
  std::vector<double> T_m(My);
  for (int m = 0; m < My; m++) {
    T_m[m] = (cum_nxgt[m] + 0.5 * cum_nxeq[m] +
      0.5 * sum_nxgt_by_ygroup[m] +
      0.25 * sum_nxeq_by_ygroup[m]) / N;
  }
  
  // Final kernel values
  NumericVector result(N);
  for (int i = 0; i < N; i++) {
    double g_1 = T_m[Yc[i] - 1] / N;
    double g_2 = S_i[i] / N;
    result[i] = 4.0 * (g_1 + g_2 + G_x[i] * G_y[i] - G_y[i] - G_x[i])
      + 1.0 - rho;
  }
  
  return result;
}
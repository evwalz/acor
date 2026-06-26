#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
using namespace Rcpp;

// ============================================================================
// Coordinate compression (doubles -> consecutive integers 1..M)
// ============================================================================

static std::vector<int> compress_d(const std::vector<double>& vals, int& M) {
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
    ) - sorted_unique.begin() + 1;
  }
  return compressed;
}


// ============================================================================
// agc_ij_cpp -- AGC IJ (sample-ratio chain rule)
//
// Used by acor.test(variance = "ij"). Differentiates the sample ratio agc = A / B
// directly w.r.t. observation weights, where
//
//   A = sum_i tQ_i tR_i,   B = sum_i tR_i^2
//
// with mid-ranks R_i, Q_i and tR_i, tQ_i their centred versions. The IF
// is assembled via the quotient rule with rank-shift corrections collected
// as prefix sums (SRx, SQy, SRy):
//
//   dA_i = tQ_i tR_i + SRx_i + SQy_i
//   dB_i = tR_i^2   + 2 SRy_i
//   IC_i = n * ( dA_i / B  -  (A / B^2) dB_i )
//   var_ij = (1/n) sum_i IC_i^2
//
// Complexity: O(n log n) via three prefix-sum sweeps over rank groups.
//
// Returns: list(agc, ic, var_ij)
// ============================================================================

// [[Rcpp::export]]
List agc_ij_cpp(NumericVector X, NumericVector Y) {
  int n = X.size();

  // ------------------------------------------------------------------
  // Step 0: mid-ranks and centred mid-ranks
  // ------------------------------------------------------------------
  std::vector<double> yvals(Y.begin(), Y.end());
  std::vector<double> xvals(X.begin(), X.end());

  auto avg_rank = [&](const std::vector<double>& v) -> std::vector<double> {
    int m = v.size();
    std::vector<int> ord(m);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
      return v[a] < v[b];
    });
    std::vector<double> ranks(m);
    int i = 0;
    while (i < m) {
      int j = i;
      while (j < m && v[ord[j]] == v[ord[i]]) j++;
      double avg = 0.5 * (i + 1 + j);
      for (int k = i; k < j; k++) ranks[ord[k]] = avg;
      i = j;
    }
    return ranks;
  };

  std::vector<double> R = avg_rank(yvals);
  std::vector<double> Q = avg_rank(xvals);

  double mean_rank = 0.5 * (n + 1);

  std::vector<double> tR(n), tQ(n);
  for (int i = 0; i < n; i++) {
    tR[i] = R[i] - mean_rank;
    tQ[i] = Q[i] - mean_rank;
  }

  // ------------------------------------------------------------------
  // Step 1: A (numerator) and B (denominator) of AGC
  // ------------------------------------------------------------------
  double A = 0.0, B = 0.0;
  for (int i = 0; i < n; i++) {
    A += tQ[i] * tR[i];
    B += tR[i] * tR[i];
  }
  double agc = A / B;

  // ------------------------------------------------------------------
  // Step 2: rank-weighted sign sums via prefix sums
  //
  //   SRx_i = (1/2) sum_j tR[j] * s^x_{ji}
  //   SQy_i = (1/2) sum_j tQ[j] * s^y_{ji}
  //   SRy_i = (1/2) sum_j tR[j] * s^y_{ji}
  //
  // with the convention s^z_{ji} = sgn(z_j - z_i) for j != i, and self-term
  // s^z_{ii} absorbed into the "+ tZ[i]" piece below.
  // ------------------------------------------------------------------
  int Mx;
  std::vector<int> Xc = compress_d(xvals, Mx);
  std::vector<std::vector<int>> X_groups(Mx);
  for (int i = 0; i < n; i++) X_groups[Xc[i] - 1].push_back(i);

  std::vector<double> SRx(n, 0.0);
  double cum_tR = 0.0;
  double total_tR = 0.0;
  for (int i = 0; i < n; i++) total_tR += tR[i];

  for (int g = 0; g < Mx; g++) {
    double group_tR = 0.0;
    for (int idx : X_groups[g]) group_tR += tR[idx];

    for (int idx : X_groups[g]) {
      double sum_gt = total_tR - cum_tR - group_tR;
      double sum_lt = cum_tR;
      SRx[idx] = 0.5 * (sum_gt - sum_lt + tR[idx]);
    }
    cum_tR += group_tR;
  }

  int My;
  std::vector<int> Yc = compress_d(yvals, My);
  std::vector<std::vector<int>> Y_groups(My);
  for (int i = 0; i < n; i++) Y_groups[Yc[i] - 1].push_back(i);

  std::vector<double> SQy(n, 0.0);
  std::vector<double> SRy(n, 0.0);
  double cum_tQ = 0.0;
  double cum_tR_y = 0.0;
  double total_tQ = 0.0;
  for (int i = 0; i < n; i++) total_tQ += tQ[i];

  for (int g = 0; g < My; g++) {
    double group_tQ = 0.0;
    double group_tR = 0.0;
    for (int idx : Y_groups[g]) {
      group_tQ += tQ[idx];
      group_tR += tR[idx];
    }

    for (int idx : Y_groups[g]) {
      double sum_gt_Q = total_tQ - cum_tQ - group_tQ;
      double sum_lt_Q = cum_tQ;
      SQy[idx] = 0.5 * (sum_gt_Q - sum_lt_Q + tQ[idx]);

      double sum_gt_R = total_tR - cum_tR_y - group_tR;
      double sum_lt_R = cum_tR_y;
      SRy[idx] = 0.5 * (sum_gt_R - sum_lt_R + tR[idx]);
    }
    cum_tQ += group_tQ;
    cum_tR_y += group_tR;
  }

  // ------------------------------------------------------------------
  // Step 3: assemble dA, dB, IC, var_ij
  //
  //   dA_i = tQ_i tR_i + SRx_i + SQy_i
  //   dB_i = tR_i^2   + 2 * SRy_i
  //   IC_i = n * ( dA_i / B  -  (A / B^2) dB_i )
  // ------------------------------------------------------------------
  NumericVector ic(n);
  double sum_ic2 = 0.0;
  double invB  = 1.0 / B;
  double AoB2  = A / (B * B);

  for (int i = 0; i < n; i++) {
    double dA_i = tQ[i] * tR[i] + SRx[i] + SQy[i];
    double dB_i = tR[i] * tR[i] + 2.0 * SRy[i];
    double ic_i = (double)n * (dA_i * invB - AoB2 * dB_i);
    ic[i] = ic_i;
    sum_ic2 += ic_i * ic_i;
  }
  double var_ij = sum_ic2 / n;

  return List::create(
    Named("agc")    = agc,
    Named("ic")     = ic,
    Named("var_ij") = var_ij
  );
}


// ============================================================================
// agc_ij_variance_cpp
//
// Simplified interface: returns just AGC and its IJ variance.
// ============================================================================

// [[Rcpp::export]]
List agc_ij_variance_cpp(NumericVector X, NumericVector Y) {
  List full = agc_ij_cpp(X, Y);
  return List::create(
    Named("agc") = full["agc"],
    Named("var") = full["var_ij"]
  );
}


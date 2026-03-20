#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
using namespace Rcpp;

// ============================================================================
// Double-valued Fenwick Tree
// Supports point updates and prefix sum queries for double values -- O(log n)
// ============================================================================

class FenwickTreeD {
  std::vector<double> tree;
  int n;
public:
  FenwickTreeD(int n) : tree(n + 1, 0.0), n(n) {}
  
  void update(int i, double delta) {
    for (; i <= n; i += i & (-i))
      tree[i] += delta;
  }
  
  double query(int i) {
    double s = 0.0;
    for (; i > 0; i -= i & (-i))
      s += tree[i];
    return s;
  }
  
  double query_range(int lo, int hi) {
    if (lo > hi) return 0.0;
    return query(hi) - (lo > 1 ? query(lo - 1) : 0.0);
  }
};


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
// agc_ij_cpp
//
// Computes AGC point estimate, IJ influence function values, and IJ variance
// using two Fenwick-tree sweeps.
//
// Inputs:
//   X  -- numeric predictor vector (raw values, not ranks)
//   Y  -- numeric outcome vector (raw values, not ranks)
//
// Returns: list with
//   agc      -- AGC point estimate
//   ic       -- IJ influence function values (length n)
//   var_ij   -- IJ variance = (1/n) * sum(ic^2)
//   dA       -- partial A / partial w_i values (for debugging)
//   dB       -- partial B / partial w_i values (for debugging)
//
// Complexity: O(n log n)
// ============================================================================

// [[Rcpp::export]]
List agc_ij_cpp(NumericVector X, NumericVector Y) {
  int n = X.size();
  
  // ------------------------------------------------------------------
  // Step 0: Compute mid-ranks and centered mid-ranks
  // ------------------------------------------------------------------
  
  // Compute average ranks for Y
  std::vector<double> yvals(Y.begin(), Y.end());
  std::vector<double> xvals(X.begin(), X.end());
  
  // Sort indices by value to compute average ranks
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
      double avg = 0.5 * (i + 1 + j);  // average of ranks i+1 .. j
      for (int k = i; k < j; k++) ranks[ord[k]] = avg;
      i = j;
    }
    return ranks;
  };
  
  std::vector<double> R = avg_rank(yvals);  // mid-ranks of Y
  std::vector<double> Q = avg_rank(xvals);  // mid-ranks of X
  
  double mean_rank = 0.5 * (n + 1);
  
  // Centered mid-ranks
  std::vector<double> tR(n), tQ(n);
  for (int i = 0; i < n; i++) {
    tR[i] = R[i] - mean_rank;
    tQ[i] = Q[i] - mean_rank;
  }
  
  // ------------------------------------------------------------------
  // Step 1: Compute A (numerator) and B (denominator) of AGC
  // ------------------------------------------------------------------
  
  double A = 0.0, B = 0.0;
  for (int i = 0; i < n; i++) {
    A += tQ[i] * tR[i];
    B += tR[i] * tR[i];
  }
  double agc = A / B;
  
  // ------------------------------------------------------------------
  // Step 2: Rank-weighted sign sums via Fenwick trees
  //
  // We need three quantities for each observation i:
  //
  //   SRx_i = (1/2) * sum_j tR[j] * sgn(x_j - x_i)
  //         = sum_{j: x_j > x_i} tR[j] - sum_{j: x_j < x_i} tR[j]
  //           (with half-weight for ties), divided by 2
  //         ... but note: sum_{x_j > x_i} tR[j] + sum_{x_j < x_i} tR[j]
  //           + sum_{x_j == x_i} tR[j] = 0 (since tR sums to 0)
  //
  //   SQy_i = (1/2) * sum_j tQ[j] * sgn(y_j - y_i)
  //
  //   SRy_i = (1/2) * sum_j tR[j] * sgn(y_j - y_i)
  //           (needed for dB)
  //
  // For SRx: sweep sorted by X, accumulate tR values in Fenwick tree
  //          keyed by... actually we just need prefix sums by X-groups.
  //
  // For SQy, SRy: sweep sorted by Y, accumulate tQ and tR values.
  // ------------------------------------------------------------------
  
  // --- Sweep 1: sorted by X, accumulate tR values ---
  // For each i, compute sum of tR[j] over {j : x_j < x_i} and {j : x_j == x_i}
  
  // Compress X values
  int Mx;
  std::vector<int> Xc = compress_d(xvals, Mx);
  std::vector<std::vector<int>> X_groups(Mx);
  for (int i = 0; i < n; i++) X_groups[Xc[i] - 1].push_back(i);
  
  // Process X-groups in sorted order (prefix sum approach, no Fenwick needed)
  std::vector<double> SRx(n, 0.0);
  double cum_tR = 0.0;       // sum of tR[j] for x_j < x_i
  double total_tR = 0.0;     // sum of all tR = 0, but compute for clarity
  for (int i = 0; i < n; i++) total_tR += tR[i];
  
  for (int g = 0; g < Mx; g++) {
    // sum of tR within this X-group
    double group_tR = 0.0;
    for (int idx : X_groups[g]) group_tR += tR[idx];
    
    for (int idx : X_groups[g]) {
      // sum_{x_j < x_i} tR[j] = cum_tR
      // sum_{x_j == x_i} tR[j] = group_tR  (includes self)
      // sum_{x_j > x_i} tR[j] = total_tR - cum_tR - group_tR
      //
      // sgn sum = (sum_{x_j > x_i} - sum_{x_j < x_i}) tR[j]
      //         = total_tR - 2*cum_tR - group_tR
      // For ties: sgn(x_j - x_i) = 0 when x_j == x_i, so tie group
      //   contributes 0 to the sign sum (except self-term handled separately).
      //
      // Actually: sgn(x_j - x_i) for x_j == x_i (and j != i) is 0.
      // Self-term s_{ii} = sgn(0) + 1 = 1 (from our definition).
      //
      // So: (1/2) sum_j tR[j] * s^x_{ji}
      //   = (1/2) [sum_{x_j>x_i} tR[j] - sum_{x_j<x_i} tR[j] + tR[i]]
      
      double sum_gt = total_tR - cum_tR - group_tR;
      double sum_lt = cum_tR;
      SRx[idx] = 0.5 * (sum_gt - sum_lt + tR[idx]);
    }
    cum_tR += group_tR;
  }
  
  // --- Sweep 2: sorted by Y, accumulate tQ and tR values ---
  
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
      // SQy: (1/2) sum_j tQ[j] * s^y_{ji}
      double sum_gt_Q = total_tQ - cum_tQ - group_tQ;
      double sum_lt_Q = cum_tQ;
      SQy[idx] = 0.5 * (sum_gt_Q - sum_lt_Q + tQ[idx]);
      
      // SRy: (1/2) sum_j tR[j] * s^y_{ji} -- note: same as sum_j tR[j] * sgn(y_j - y_i) / 2 + tR[i]/2
      // but we can also express this: since sum tR = 0,
      //   SRy = 0.5 * (sum_{y_j > y_i} tR[j] - sum_{y_j < y_i} tR[j] + tR[i])
      double sum_gt_R = total_tR - cum_tR_y - group_tR;
      double sum_lt_R = cum_tR_y;
      SRy[idx] = 0.5 * (sum_gt_R - sum_lt_R + tR[idx]);
    }
    cum_tQ += group_tQ;
    cum_tR_y += group_tR;
  }
  
  // ------------------------------------------------------------------
  // Step 3: Assemble dA/dw_i and dB/dw_i, then IC_i
  // ------------------------------------------------------------------
  
  NumericVector dA_vec(n), dB_vec(n), ic(n);
  double sum_ic2 = 0.0;
  
  for (int i = 0; i < n; i++) {
    // dA/dw_i = tQ_i * tR_i + SRx_i + SQy_i
    double dA_i = tQ[i] * tR[i] + SRx[i] + SQy[i];
    
    // dB/dw_i = tR_i^2 + 2 * SRy_i
    //   Derivation:
    //   B = sum_j tR_j^2
    //   dB/dw_i = tR_i^2 + sum_j 2 * tR_j * (dtR_j/dw_i)
    //           = tR_i^2 + 2 * sum_j tR_j * (1/2 s^y_{ji} - c_R(i))
    //           = tR_i^2 + sum_j tR_j * s^y_{ji}    (centering drops out)
    //           = tR_i^2 + 2 * SRy_i
    double dB_i = tR[i] * tR[i] + 2.0 * SRy[i];
    
    // Raw derivative: d(AGC)/dw_i = dA_i/B - (A/B^2)*dB_i
    // This is O(1/n) since A,B are O(n).
    //
    // The IJ influence function is IC_i = n * d(AGC)/dw_i, which is O(1).
    // The asymptotic variance is then: sigma^2 = (1/n) * sum IC_i^2
    // This matches the convention of Sigma_agc_v2 where Var(AGC) = sigma^2/n.
    double raw_deriv = dA_i / B - (A / (B * B)) * dB_i;
    double ic_i = (double)n * raw_deriv;
    
    dA_vec[i] = dA_i;
    dB_vec[i] = dB_i;
    ic[i] = ic_i;
    sum_ic2 += ic_i * ic_i;
  }
  
  // Asymptotic variance: sigma^2 = (1/n) sum IC_i^2
  // so that Var(AGC) = sigma^2 / n
  double var_ij = sum_ic2 / n;
  
  return List::create(
    Named("agc") = agc,
    Named("ic") = ic,
    Named("var_ij") = var_ij,
    Named("dA") = dA_vec,
    Named("dB") = dB_vec
  );
}


// ============================================================================
// agc_ij_variance_cpp
//
// Simplified interface: returns just AGC and its IJ variance.
// For use in acor.test() as a drop-in.
// ============================================================================

// [[Rcpp::export]]
List agc_ij_variance_cpp(NumericVector X, NumericVector Y) {
  List full = agc_ij_cpp(X, Y);
  return List::create(
    Named("agc") = full["agc"],
                       Named("var") = full["var_ij"]
  );
}



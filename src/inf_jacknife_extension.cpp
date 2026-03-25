#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <vector>

using namespace Rcpp;

// ============================================================================
// Shared utilities (rho_b, tau_b, gamma IJ)
// ============================================================================

class FenwickTreeIJ {
  std::vector<int> tree;
  int n;

 public:
  explicit FenwickTreeIJ(int n_in) : tree(n_in + 1, 0), n(n_in) {}
  void update(int i, int delta = 1) {
    for (; i <= n; i += i & (-i)) tree[i] += delta;
  }
  int query(int i) {
    int s = 0;
    for (; i > 0; i -= i & (-i)) s += tree[i];
    return s;
  }
};

static std::vector<int> compress_ext(const std::vector<double>& vals, int& M) {
  int nv = vals.size();
  std::vector<double> sorted_unique = vals;
  std::sort(sorted_unique.begin(), sorted_unique.end());
  sorted_unique.erase(std::unique(sorted_unique.begin(), sorted_unique.end()),
                      sorted_unique.end());
  M = sorted_unique.size();
  std::vector<int> compressed(nv);
  for (int i = 0; i < nv; i++) {
    compressed[i] = static_cast<int>(
        std::lower_bound(sorted_unique.begin(), sorted_unique.end(), vals[i]) -
        sorted_unique.begin() + 1);
  }
  return compressed;
}

static void fill_na_ic(NumericVector ic, double* sum_ic2) {
  int n = ic.size();
  for (int i = 0; i < n; i++) ic[i] = NA_REAL;
  *sum_ic2 = NA_REAL;
}

// ============================================================================
// rho_b_ij_cpp
// ============================================================================

// [[Rcpp::export]]
List rho_b_ij_cpp(NumericVector X, NumericVector Y) {
  int n = X.size();
  std::vector<double> yvals(Y.begin(), Y.end());
  std::vector<double> xvals(X.begin(), X.end());

  auto avg_rank = [&](const std::vector<double>& v) -> std::vector<double> {
    int m = v.size();
    std::vector<int> ord(m);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(),
              [&](int a, int b) { return v[a] < v[b]; });
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

  double A = 0.0, BX = 0.0, BY = 0.0;
  for (int i = 0; i < n; i++) {
    A += tQ[i] * tR[i];
    BX += tQ[i] * tQ[i];
    BY += tR[i] * tR[i];
  }
  double sqrt_BXBY0 = std::sqrt(BX * BY);
  const double tol_rb = 1e-10;
  double rho_b = (sqrt_BXBY0 > tol_rb && BX > tol_rb && BY > tol_rb)
                     ? A / sqrt_BXBY0
                     : 0.0;

  int Mx, My;
  std::vector<int> Xc = compress_ext(xvals, Mx);
  std::vector<int> Yc = compress_ext(yvals, My);

  std::vector<std::vector<int>> X_groups(Mx);
  for (int i = 0; i < n; i++) X_groups[Xc[i] - 1].push_back(i);

  std::vector<std::vector<int>> Y_groups(My);
  for (int i = 0; i < n; i++) Y_groups[Yc[i] - 1].push_back(i);

  double total_tR = 0.0, total_tQ = 0.0;
  for (int i = 0; i < n; i++) {
    total_tR += tR[i];
    total_tQ += tQ[i];
  }

  std::vector<double> SRx(n, 0.0);
  std::vector<double> SQx(n, 0.0);
  double cum_tR_x = 0.0;
  double cum_tQ_x = 0.0;

  for (int g = 0; g < Mx; g++) {
    double group_tR = 0.0;
    double group_tQ = 0.0;
    for (int idx : X_groups[g]) {
      group_tR += tR[idx];
      group_tQ += tQ[idx];
    }
    for (int idx : X_groups[g]) {
      double sum_gt_R = total_tR - cum_tR_x - group_tR;
      double sum_lt_R = cum_tR_x;
      SRx[idx] = 0.5 * (sum_gt_R - sum_lt_R + tR[idx]);
      double sum_gt_Q = total_tQ - cum_tQ_x - group_tQ;
      double sum_lt_Q = cum_tQ_x;
      SQx[idx] = 0.5 * (sum_gt_Q - sum_lt_Q + tQ[idx]);
    }
    cum_tR_x += group_tR;
    cum_tQ_x += group_tQ;
  }

  std::vector<double> SQy(n, 0.0);
  std::vector<double> SRy(n, 0.0);
  double cum_tQ_y = 0.0;
  double cum_tR_y = 0.0;

  for (int g = 0; g < My; g++) {
    double group_tQ = 0.0;
    double group_tR = 0.0;
    for (int idx : Y_groups[g]) {
      group_tQ += tQ[idx];
      group_tR += tR[idx];
    }
    for (int idx : Y_groups[g]) {
      double sum_gt_Q = total_tQ - cum_tQ_y - group_tQ;
      double sum_lt_Q = cum_tQ_y;
      SQy[idx] = 0.5 * (sum_gt_Q - sum_lt_Q + tQ[idx]);
      double sum_gt_R = total_tR - cum_tR_y - group_tR;
      double sum_lt_R = cum_tR_y;
      SRy[idx] = 0.5 * (sum_gt_R - sum_lt_R + tR[idx]);
    }
    cum_tQ_y += group_tQ;
    cum_tR_y += group_tR;
  }

  NumericVector dA_vec(n), dBX_vec(n), dBY_vec(n), ic(n);
  double sum_ic2 = 0.0;
  double sqrt_BXBY = sqrt_BXBY0;

  const double tol = tol_rb;
  if (sqrt_BXBY < tol || BX < tol || BY < tol) {
    fill_na_ic(ic, &sum_ic2);
    for (int i = 0; i < n; i++) {
      dA_vec[i] = NA_REAL;
      dBX_vec[i] = NA_REAL;
      dBY_vec[i] = NA_REAL;
    }
  } else {
    for (int i = 0; i < n; i++) {
      double dA_i = tQ[i] * tR[i] + SRx[i] + SQy[i];
      double dBX_i = tQ[i] * tQ[i] + 2.0 * SQx[i];
      double dBY_i = tR[i] * tR[i] + 2.0 * SRy[i];
      double raw_deriv = dA_i / sqrt_BXBY -
                         (A * dBX_i) / (2.0 * BX * sqrt_BXBY) -
                         (A * dBY_i) / (2.0 * BY * sqrt_BXBY);
      double ic_i = static_cast<double>(n) * raw_deriv;
      dA_vec[i] = dA_i;
      dBX_vec[i] = dBX_i;
      dBY_vec[i] = dBY_i;
      ic[i] = ic_i;
      sum_ic2 += ic_i * ic_i;
    }
  }

  double var_ij = (R_IsNA(sum_ic2)) ? NA_REAL : (sum_ic2 / n);

  return List::create(
      Named("rho_b") = rho_b, Named("ic") = ic, Named("var_ij") = var_ij,
      Named("dA") = dA_vec, Named("dBX") = dBX_vec, Named("dBY") = dBY_vec);
}

// ============================================================================
// tau_b_ij_cpp
// ============================================================================

// [[Rcpp::export]]
List tau_b_ij_cpp(NumericVector X, NumericVector Y) {
  int n = X.size();
  std::vector<double> xvals(X.begin(), X.end());
  std::vector<double> yvals(Y.begin(), Y.end());

  int Mx, My;
  std::vector<int> Xc = compress_ext(xvals, Mx);
  std::vector<int> Yc = compress_ext(yvals, My);

  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
    if (xvals[a] != xvals[b]) return xvals[a] < xvals[b];
    return yvals[a] < yvals[b];
  });

  std::vector<double> H_A(n, 0.0);
  std::vector<double> H_BY(n, 0.0);
  std::vector<double> H_BX(n, 0.0);

  {
    FenwickTreeIJ tree(My);
    int i = 0;
    while (i < n) {
      int j = i;
      while (j < n && xvals[ord[j]] == xvals[ord[i]]) j++;
      int total_inserted = (i > 0) ? tree.query(My) : 0;
      for (int k = i; k < j; k++) {
        int idx = ord[k];
        int yc = Yc[idx];
        int below = (yc > 1) ? tree.query(yc - 1) : 0;
        int above = total_inserted - tree.query(yc);
        H_A[idx] += static_cast<double>(below - above);
        H_BY[idx] += static_cast<double>(below + above);
      }
      for (int k = i; k < j; k++) tree.update(Yc[ord[k]]);
      i = j;
    }
  }

  {
    FenwickTreeIJ tree(My);
    int i = n;
    while (i > 0) {
      int j = i;
      while (j > 0 && xvals[ord[j - 1]] == xvals[ord[i - 1]]) j--;
      int total_inserted = tree.query(My);
      for (int k = j; k < i; k++) {
        int idx = ord[k];
        int yc = Yc[idx];
        int below = (yc > 1) ? tree.query(yc - 1) : 0;
        int above = total_inserted - tree.query(yc);
        H_A[idx] += static_cast<double>(above - below);
        H_BY[idx] += static_cast<double>(below + above);
      }
      for (int k = j; k < i; k++) tree.update(Yc[ord[k]]);
      i = j;
    }
  }

  {
    std::vector<std::vector<int>> X_groups(Mx);
    for (int i = 0; i < n; i++) X_groups[Xc[i] - 1].push_back(i);
    for (int g = 0; g < Mx; g++) {
      int gsize = X_groups[g].size();
      if (gsize <= 1) continue;
      std::map<int, int> y_freq;
      for (int idx : X_groups[g]) y_freq[Yc[idx]]++;
      for (int idx : X_groups[g]) {
        int same_y = y_freq[Yc[idx]] - 1;
        int diff_y = (gsize - 1) - same_y;
        H_BY[idx] += static_cast<double>(diff_y);
      }
    }
  }

  {
    std::map<double, int> x_freq;
    for (int i = 0; i < n; i++) x_freq[xvals[i]]++;
    for (int i = 0; i < n; i++) H_BX[i] = static_cast<double>(n - x_freq[xvals[i]]);
  }

  double num_pairs = static_cast<double>(n) * (n - 1) / 2.0;
  double sum_HA = 0.0, sum_HBX = 0.0, sum_HBY = 0.0;
  for (int i = 0; i < n; i++) {
    sum_HA += H_A[i];
    sum_HBX += H_BX[i];
    sum_HBY += H_BY[i];
  }

  double A = sum_HA / (2.0 * num_pairs);
  double BX = sum_HBX / (2.0 * num_pairs);
  double BY = sum_HBY / (2.0 * num_pairs);
  double sqrt_BXBY = std::sqrt(BX * BY);

  const double tol = 1e-10;
  double tau_b = (sqrt_BXBY > tol) ? A / sqrt_BXBY : 0.0;

  NumericVector ic(n);
  NumericVector dA_vec(n), dBX_vec(n), dBY_vec(n);
  double sum_ic2 = 0.0;

  if (sqrt_BXBY < tol || BX < tol || BY < tol) {
    fill_na_ic(ic, &sum_ic2);
    for (int i = 0; i < n; i++) {
      dA_vec[i] = NA_REAL;
      dBX_vec[i] = NA_REAL;
      dBY_vec[i] = NA_REAL;
    }
  } else {
    for (int i = 0; i < n; i++) {
      double n_dA_i = 2.0 * (H_A[i] / (n - 1) - A);
      double n_dBX_i = 2.0 * (H_BX[i] / (n - 1) - BX);
      double n_dBY_i = 2.0 * (H_BY[i] / (n - 1) - BY);
      double ic_i = n_dA_i / sqrt_BXBY -
                     (A / (2.0 * BX * sqrt_BXBY)) * n_dBX_i -
                     (A / (2.0 * BY * sqrt_BXBY)) * n_dBY_i;
      dA_vec[i] = n_dA_i;
      dBX_vec[i] = n_dBX_i;
      dBY_vec[i] = n_dBY_i;
      ic[i] = ic_i;
      sum_ic2 += ic_i * ic_i;
    }
  }

  double var_ij = (R_IsNA(sum_ic2)) ? NA_REAL : (sum_ic2 / n);

  return List::create(Named("tau_b") = tau_b, Named("ic") = ic,
                      Named("var_ij") = var_ij, Named("dA") = dA_vec,
                      Named("dBX") = dBX_vec, Named("dBY") = dBY_vec);
}

// ============================================================================
// gamma_ij_cpp
// ============================================================================

// [[Rcpp::export]]
List gamma_ij_cpp(NumericVector X, NumericVector Y) {
  int n = X.size();
  std::vector<double> xvals(X.begin(), X.end());
  std::vector<double> yvals(Y.begin(), Y.end());

  int My;
  std::vector<int> Yc = compress_ext(yvals, My);

  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
    if (xvals[a] != xvals[b]) return xvals[a] < xvals[b];
    return yvals[a] < yvals[b];
  });

  std::vector<double> H_A(n, 0.0);

  {
    FenwickTreeIJ tree(My);
    int i = 0;
    while (i < n) {
      int j = i;
      while (j < n && xvals[ord[j]] == xvals[ord[i]]) j++;
      int total_inserted = (i > 0) ? tree.query(My) : 0;
      for (int k = i; k < j; k++) {
        int idx = ord[k];
        int yc = Yc[idx];
        int below = (yc > 1) ? tree.query(yc - 1) : 0;
        int above = total_inserted - tree.query(yc);
        H_A[idx] += static_cast<double>(below - above);
      }
      for (int k = i; k < j; k++) tree.update(Yc[ord[k]]);
      i = j;
    }
  }

  {
    FenwickTreeIJ tree(My);
    int i = n;
    while (i > 0) {
      int j = i;
      while (j > 0 && xvals[ord[j - 1]] == xvals[ord[i - 1]]) j--;
      int total_inserted = tree.query(My);
      for (int k = j; k < i; k++) {
        int idx = ord[k];
        int yc = Yc[idx];
        int below = (yc > 1) ? tree.query(yc - 1) : 0;
        int above = total_inserted - tree.query(yc);
        H_A[idx] += static_cast<double>(above - below);
      }
      for (int k = j; k < i; k++) tree.update(Yc[ord[k]]);
      i = j;
    }
  }

  std::vector<double> H_nu(n, 0.0);
  std::map<double, int> freq_x;
  for (int i = 0; i < n; i++) freq_x[xvals[i]]++;
  std::map<double, int> freq_y;
  for (int i = 0; i < n; i++) freq_y[yvals[i]]++;
  std::map<std::pair<double, double>, int> freq_xy;
  for (int i = 0; i < n; i++)
    freq_xy[std::make_pair(xvals[i], yvals[i])]++;

  for (int i = 0; i < n; i++) {
    int ties_x = freq_x[xvals[i]] - 1;
    int ties_y = freq_y[yvals[i]] - 1;
    int ties_xy = freq_xy[std::make_pair(xvals[i], yvals[i])] - 1;
    H_nu[i] = static_cast<double>(ties_x + ties_y - ties_xy);
  }

  double num_pairs = static_cast<double>(n) * (n - 1) / 2.0;
  double sum_HA = 0.0, sum_Hnu = 0.0;
  for (int i = 0; i < n; i++) {
    sum_HA += H_A[i];
    sum_Hnu += H_nu[i];
  }

  double A = sum_HA / (2.0 * num_pairs);
  double nu = sum_Hnu / (2.0 * num_pairs);
  double D = 1.0 - nu;

  const double tol = 1e-10;
  double gamma_val = (D > tol) ? A / D : 0.0;

  NumericVector ic(n);
  NumericVector dA_vec(n), dnu_vec(n);
  double sum_ic2 = 0.0;

  if (D < tol) {
    fill_na_ic(ic, &sum_ic2);
    for (int i = 0; i < n; i++) {
      dA_vec[i] = NA_REAL;
      dnu_vec[i] = NA_REAL;
    }
  } else {
    for (int i = 0; i < n; i++) {
      double n_dA_i = 2.0 * (H_A[i] / (n - 1) - A);
      double n_dnu_i = 2.0 * (H_nu[i] / (n - 1) - nu);
      double ic_i = n_dA_i / D + (A / (D * D)) * n_dnu_i;
      dA_vec[i] = n_dA_i;
      dnu_vec[i] = n_dnu_i;
      ic[i] = ic_i;
      sum_ic2 += ic_i * ic_i;
    }
  }

  double var_ij = (R_IsNA(sum_ic2)) ? NA_REAL : (sum_ic2 / n);

  return List::create(Named("gamma") = gamma_val, Named("ic") = ic,
                      Named("var_ij") = var_ij, Named("dA") = dA_vec,
                      Named("dnu") = dnu_vec);
}

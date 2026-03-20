#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <map>
#include <cmath>
using namespace Rcpp;

// ============================================================================
  // Fenwick Tree (int-valued, for counting)
// ============================================================================
  
  class FenwickTreeI {
    std::vector<int> tree;
    int n;
    public:
      FenwickTreeI(int n) : tree(n + 1, 0), n(n) {}
    void update(int i, int delta = 1) {
      for (; i <= n; i += i & (-i)) tree[i] += delta;
    }
    int query(int i) {
      int s = 0;
      for (; i > 0; i -= i & (-i)) s += tree[i];
      return s;
    }
  };


// ============================================================================
  // Coordinate compression (doubles -> consecutive integers 1..M)
// ============================================================================
  
  static std::vector<int> compress_ij(const std::vector<double>& vals, int& M) {
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
  // akc_ij_cpp
//
  // Computes AKC point estimate, IJ influence function values, and IJ variance.
//
  // AKC = tau_XY / (1 - p_Y) = A / B  where
  //   A = (2 / n(n-1)) * sum_{i<j} sgn(x_i - x_j) * sgn(y_i - y_j)
  //   B = (2 / n(n-1)) * sum_{i<j} (1 - 1(y_i == y_j))
  //     = 1 - p_Y
  //
    // Both A and B are degree-2 U-statistics with normalisation 2/(n(n-1)).
  //
    // The weighted versions under F_w are:
    //   A(w) = sum_{i<j} w_i w_j sgn(x_i-x_j) sgn(y_i-y_j) / sum_{i<j} w_i w_j
    //   B(w) = sum_{i<j} w_i w_j (1 - 1(y_i==y_j))          / sum_{i<j} w_i w_j
    //
      // IJ influence: IC_i = n * d(AKC)/dw_i |_{w=1}
      //
        // For a ratio A/B of two degree-2 U-statistics sharing the same denominator
      // sum_{i<j} w_i w_j, the derivative simplifies to:
        //
        //   d(A/B)/dw_i = (1/B) dA/dw_i - (A/B^2) dB/dw_i
        //
          // where dA/dw_i and dB/dw_i each involve two terms:
          //   (a) the sum over j of the kernel evaluated at (z_i, z_j)
        //   (b) a centering term from differentiating the normalisation
        //
          // Specifically, for U = sum_{i<j} w_i w_j h(z_i,z_j) / sum_{i<j} w_i w_j:
            //   dU/dw_i = [sum_j h(z_i,z_j) * w_j - U * sum_j w_j] / sum_{i<j} w_i w_j
            // At w=1:
              //   dU/dw_i = [H_i - U*(n-1)] / (n(n-1)/2)
              // where H_i = sum_{j != i} h(z_i, z_j).
              //
                // Therefore:
                //   n * dU/dw_i = 2 * [H_i/(n-1) - U] / (1)  ... simplifying
                //   n * dU/dw_i = 2 * (H_i_bar - U)
                // where H_i_bar = H_i / (n-1) is the leave-one-in conditional mean.
                //
                  // This gives the degree-2 U-statistic Hoeffding projection multiplier of 2.
                //
                  // Complexity: O(n log n) via Fenwick tree for concordance sums.
                // ============================================================================
                  
                  // [[Rcpp::export]]
                List akc_ij_cpp(NumericVector X, NumericVector Y) {
                  int n = X.size();
                  
                  std::vector<double> xvals(X.begin(), X.end());
                  std::vector<double> yvals(Y.begin(), Y.end());
                  
                  // ------------------------------------------------------------------
                    // Step 0: Compute AKC point estimate
                  // ------------------------------------------------------------------
                    
                    // Compress Y for Fenwick tree
                  int My;
                  std::vector<int> Yc = compress_ij(yvals, My);
                  
                  // Sort indices by (X, Y) for concordance counting
                  std::vector<int> ord(n);
                  std::iota(ord.begin(), ord.end(), 0);
                  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
                    if (xvals[a] != xvals[b]) return xvals[a] < xvals[b];
                    return yvals[a] < yvals[b];
                  });
                  
                  // Group by X value
                  // For each observation i, compute:
                    //   conc_i = number of j != i with x_j < x_i, y_j < y_i  (concordant)
                  //   disc_i = number of j != i with x_j < x_i, y_j > y_i  (discordant)
                  // Using Fenwick tree: process X-groups in order, query before insert.
                  
                  // We need per-observation concordance and discordance counts.
                  // conc_i counts j with (x_j < x_i and y_j < y_i) 
                  // disc_i counts j with (x_j < x_i and y_j > y_i)
                  // For the sign sum: sum_{j!=i} sgn(x_i-x_j)*sgn(y_i-y_j)
                  //   = (conc_i + conc_i') - (disc_i + disc_i')
                  // where conc_i' counts j with x_j > x_i, y_j > y_i (= conc from j's perspective)
//
  // Actually, it's simpler: for degree-2 U-stat with symmetric-in-pair kernel,
  // H_i = sum_{j != i} h(z_i, z_j) where h(z_i,z_j) = sgn(x_i-x_j)*sgn(y_i-y_j).
  // We can compute this as:
  //   H_i = (conc_below_i - disc_below_i) + (conc_above_i - disc_above_i)
  // where "below" means x_j < x_i and "above" means x_j > x_i.
  //
  // By antisymmetry: conc_above_i = disc_below_from_above = (number of j with
  //   x_j > x_i, y_j > y_i). We can get this by doing two passes or by noting
  //   that the total concordant pairs involving i from both sides equals:
  //   H_i = sum_{j: x_j < x_i} [conc(y_j < y_i) - disc(y_j > y_i)]
  //        + sum_{j: x_j > x_i} [conc(y_j > y_i) - disc(y_j < y_i)]
  //   But by symmetry of the kernel, we can compute:
  //     pass 1 (forward): for each i, query Fenwick for y_j < y_i (conc) and y_j > y_i (disc) among x_j < x_i
  //     pass 2 (backward): same for x_j > x_i
  
  // --- Forward pass: count conc/disc from observations with smaller x ---
  std::vector<double> H_A(n, 0.0);  // sum_{j!=i} sgn(x_i-x_j)*sgn(y_i-y_j)
  std::vector<double> H_B(n, 0.0);  // sum_{j!=i} (1 - 1(y_i==y_j))
  
  {
    FenwickTreeI tree(My);
    int i = 0;
    while (i < n) {
      int j = i;
      // Find the X-group boundary
      while (j < n && xvals[ord[j]] == xvals[ord[i]]) j++;
      
      // Query for each element in this group (x_j < x_i for all previous groups)
      int total_inserted = (i > 0) ? tree.query(My) : 0;
      for (int k = i; k < j; k++) {
        int idx = ord[k];
        int yc = Yc[idx];
        int below = (yc > 1) ? tree.query(yc - 1) : 0;  // y_j < y_i
        int above = total_inserted - tree.query(yc);       // y_j > y_i
        
        // sgn(x_i-x_j)=+1 since x_j < x_i.
        // y_j < y_i: sgn(y_i-y_j)=+1, product=+1 (concordant)
        // y_j > y_i: sgn(y_i-y_j)=-1, product=-1 (discordant)
        H_A[idx] += (double)(below - above);
        
        // 1 - 1(y_i==y_j): count pairs with y_j != y_i
        H_B[idx] += (double)(below + above);
      }
      
      // Insert all elements of this group
      for (int k = i; k < j; k++) {
        tree.update(Yc[ord[k]]);
      }
      i = j;
    }
  }
  
  // --- Backward pass: count conc/disc from observations with larger x ---
  {
    FenwickTreeI tree(My);
    int i = n;
    while (i > 0) {
      int j = i;
      // Find the X-group boundary (going backwards)
      while (j > 0 && xvals[ord[j-1]] == xvals[ord[i-1]]) j--;
      
      // Query for each element in this group (x_j > x_i for all later groups)
      int total_inserted = tree.query(My);
      for (int k = j; k < i; k++) {
        int idx = ord[k];
        int yc = Yc[idx];
        int below = (yc > 1) ? tree.query(yc - 1) : 0;  // y_j < y_i
        int above = total_inserted - tree.query(yc);       // y_j > y_i
        
        // sgn(x_i-x_j)=-1 since x_j > x_i.
        // y_j < y_i: sgn(y_i-y_j)=+1, product=-1 (discordant)
        // y_j > y_i: sgn(y_i-y_j)=-1, product=+1 (concordant)
        H_A[idx] += (double)(above - below);
        
        // 1 - 1(y_i==y_j): count pairs with y_j != y_i
        H_B[idx] += (double)(below + above);
      }
      
      // Insert all elements of this group
      for (int k = j; k < i; k++) {
        tree.update(Yc[ord[k]]);
      }
      i = j;
    }
  }
  
  // Now also add the x_j == x_i (same X-group, j != i) contribution to H_B.
  // For H_A: sgn(x_i - x_j) = 0 when x_j == x_i, so no contribution.
  // For H_B: need to count j in same X-group with y_j != y_i.
  {
    // Compress X
    int Mx;
    std::vector<int> Xc = compress_ij(xvals, Mx);
    
    // Group by X
    std::vector<std::vector<int>> X_groups(Mx);
    for (int i = 0; i < n; i++) X_groups[Xc[i] - 1].push_back(i);
    
    for (int g = 0; g < Mx; g++) {
      int gsize = X_groups[g].size();
      if (gsize <= 1) continue;
      
      // Count Y-value frequencies within this X-group
      std::map<int, int> y_freq;
      for (int idx : X_groups[g]) y_freq[Yc[idx]]++;
      
      for (int idx : X_groups[g]) {
        int same_y = y_freq[Yc[idx]] - 1;  // exclude self
        int diff_y = (gsize - 1) - same_y;
        H_B[idx] += (double)diff_y;
      }
    }
  }
  
  // ------------------------------------------------------------------
  // Step 1: Compute A, B, AKC
  // ------------------------------------------------------------------
  
  double num_pairs = (double)n * (n - 1) / 2.0;
  
  // A = (1/num_pairs) * sum_{i<j} sgn(x_i-x_j)*sgn(y_i-y_j)
  //   = (1/num_pairs) * 0.5 * sum_i H_A[i]  (since each pair counted from both sides)
  double sum_HA = 0.0, sum_HB = 0.0;
  for (int i = 0; i < n; i++) {
    sum_HA += H_A[i];
    sum_HB += H_B[i];
  }
  
  double A = sum_HA / (2.0 * num_pairs);  // tau_XY
  double B = sum_HB / (2.0 * num_pairs);  // 1 - p_Y
  double akc = (B > 1e-10) ? A / B : 0.0;
  
  // ------------------------------------------------------------------
  // Step 2: Compute IJ influence function
  // ------------------------------------------------------------------
  //
  // For U = sum_{i<j} w_i w_j h(z_i,z_j) / sum_{i<j} w_i w_j,
  // at w = 1:
  //   n * dU/dw_i = 2 * (H_i/(n-1) - U)
  //
  // where H_i = sum_{j!=i} h(z_i,z_j) and U is the U-statistic value.
  // The factor 2 is the degree-2 Hoeffding projection multiplier.
  //
  // Then by quotient rule on AKC = A/B:
  //   IC_i = n * d(A/B)/dw_i = (1/B) * n*dA/dw_i - (A/B^2) * n*dB/dw_i
  
  NumericVector ic(n);
  NumericVector dA_vec(n), dB_vec(n);
  double sum_ic2 = 0.0;
  
  for (int i = 0; i < n; i++) {
    // H_A[i] / (n-1) is the conditional mean of h_A given z_i
    double n_dA_i = 2.0 * (H_A[i] / (n - 1) - A);  // n * dA/dw_i
    double n_dB_i = 2.0 * (H_B[i] / (n - 1) - B);  // n * dB/dw_i
    
    double ic_i = n_dA_i / B - (A / (B * B)) * n_dB_i;
    
    dA_vec[i] = n_dA_i;
    dB_vec[i] = n_dB_i;
    ic[i] = ic_i;
    sum_ic2 += ic_i * ic_i;
  }
  
  // Asymptotic variance: sigma^2 = (1/n) sum IC_i^2
  double var_ij = sum_ic2 / n;
  
  return List::create(
    Named("akc") = akc,
    Named("ic") = ic,
    Named("var_ij") = var_ij,
    Named("dA") = dA_vec,
    Named("dB") = dB_vec
  );
}


// ============================================================================
// akc_ij_variance_cpp
//
// Simplified interface: returns just AKC and its IJ variance.
// ============================================================================

// [[Rcpp::export]]
List akc_ij_variance_cpp(NumericVector X, NumericVector Y) {
  List full = akc_ij_cpp(X, Y);
  return List::create(
    Named("akc") = full["akc"],
    Named("var") = full["var_ij"]
  );
}
#include <vector>
#include <iostream>
#include <Eigen/SparseCore>

#ifdef __cplusplus
extern "C" {
#endif
  void compute_Q(int size, double *ret, double *entries_sp_mat,
                int *i_sp_mat, int *j_sp_mat, int n_nonzero){
          typedef Eigen::Triplet<double> Trip;
                        std::vector<Trip> trp_sp_mat;
                        int k;

                        Eigen::SparseMatrix<double> sp_mat(size,size);

                        for(k = 0; k < n_nonzero; k++){
                                trp_sp_mat.push_back(Trip(i_sp_mat[k],j_sp_mat[k],entries_sp_mat[k]));
                        }

                        sp_mat.setFromTriplets(trp_sp_mat.begin(), trp_sp_mat.end());

                        std::cout << sp_mat << std::endl;     

    
  }
#ifdef __cplusplus
}
#endif
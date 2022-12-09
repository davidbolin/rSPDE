
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

extern "C" void compute_Q(int size, double *entries_C, int *i_C, int *j_C, 
                        double *entries_G, int *i_G, int *j_G,
                        double *entries_B_kappa, double *entries_B_tau,
                        int ncol_B, int rspde_order, 
                        double *theta_entries, double *rat_p, double *rat_r, 
                        double rat_k, int m_alpha,
                        double *Q_out, int *i_Q, int *j_Q,
                        int *graph_i, int *graph_j, int M);

void test_function(double *entries, int *i, int *j, int size) {
    typedef Eigen::Triplet<int> Trip;
    std::vector<Trip> trp;
    int k;

    // Eigen::MatrixXd A(size,size);

    Eigen::SparseMatrix<double> A(size, size);

    for(k = 0; k < size; k++){
        trp.push_back(Trip(i[k],j[k],entries[k]));
    }

    for(k = 0; k < size-1; k++){
        trp.push_back(Trip(i[k+1],j[k],entries[k]));
    }

    
    A.setFromTriplets(trp.begin(), trp.end());

    // for(k = 0; k < size; k++){
    //     A(k,k) = entries[k];
    // }

    Eigen::VectorXd v_tmp(size), w_tmp(size);

    // v_tmp = A.diagonal();

    for(k = 0; k < size; k++){
        v_tmp(k) = entries[k];
    }

    w_tmp = v_tmp.array().log();

    std::cout << A << std::endl;
    
    // w_tmp = A * v_tmp;


    // std::cout << w_tmp << std::endl;

    // for (int k=0; k < A.outerSize(); ++k)
    // {
    //     for (Eigen::SparseMatrix<int>::InnerIterator it(A,k); it; ++it)
    //     {
    //         entries[k] = it.value();
    //     }
    // }
}

// int main()
// {


    // double entries[4];
    // int i[4], j[4];

    // int k;

    // i[0] = 0;
    // i[1] = 1;
    // i[2] = 2;
    // i[3] = 3;

    // j[0] = 0;
    // j[1] = 1;
    // j[2] = 2;
    // j[3] = 3;

    // entries[0] = 1;
    // entries[1] = 10;
    // entries[2] = 100;
    // entries[3] = 1000;



    // test_function(entries, i, j, 4);

    // for(k = 0; k < 4; k++){
    //     printf("Entries[%d] = %f\n", k+1, entries[k]);
    // }
    // for(k = 0; k<4; k++){
    //     std::cout << entries[k] << std::endl;
    // }



    // typedef Eigen::Triplet<int> Trip;
    // std::vector<Trip> trp;
    
    // trp.push_back(Trip(1,1,10)); // (index, index, value)
    // trp.push_back(Trip(2,2,100));
    // trp.push_back(Trip(0,0,1));
    // trp.push_back(Trip(3,3,1000));

    // int rows, cols;
    // rows = cols = 4;
    // Eigen::SparseMatrix<int> A(rows,cols);

    // A.setFromTriplets(trp.begin(), trp.end());
    // std::cout << A << std::endl;



    // std::cout << "Row\tCol\tVal" <<std::endl;
    // for (int k=0; k < A.outerSize(); ++k)
    // {
    //     for (Eigen::SparseMatrix<int>::InnerIterator it(A,k); it; ++it)
    //     {
    //         std::cout << it.row() << "\t"; // row index
    //         std::cout << it.col() << "\t"; // col index (here it is equal to k)
    //         std::cout << it.value() << std::endl;
    //     }
    // }
//     return 0;
// }

// clang++ -I/opt/homebrew/opt/eigen/include/eigen3 -stdlib=libc++ -O test.cpp -o test
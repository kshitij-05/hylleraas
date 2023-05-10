//
// Created by Kshitij Surjuse on 5/9/23.
//

#ifndef HYLLERAAS_EIGEN_SOVE_H
#define HYLLERAAS_EIGEN_SOVE_H

#include <Eigen/Dense>

namespace hylleraas {

    typedef long double doubleT;
    typedef Eigen::Matrix<doubleT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

    doubleT solve_HC_SCE(const Matrix& H,const Matrix& S, bool verbose = false){
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver(H,S);
        auto gs_energy = solver.eigenvalues()[0];
        if(verbose){
            std::cout << "ground state energy = " << gs_energy << std::endl;
        }
        return gs_energy;
    }
}


#endif //HYLLERAAS_EIGEN_SOVE_H

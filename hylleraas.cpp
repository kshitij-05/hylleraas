//
// Created by Kshitij Surjuse on 5/9/23.
//

#include <iostream>
#include "basis.h"
#include "integrals.h"
#include "eigen_solve.h"
#include <Eigen/Dense>

using namespace hylleraas;

int main(int argc, char* argv[]){

    std::cout << std::setprecision(15);

    numeric_type alpha = 3.6;  // actually 2 * 1.8
    numeric_type gamma = 0.0;

    std::cout << "N\tbasis size\tenergy (Eh)"<< std::endl;
    for(auto N=0;N<14;++N){
        auto b_set = lambda_N(N,alpha,gamma);
        auto S = compute_integral(b_set, &S_ij);
        auto Vne = compute_integral(b_set, &Vne_ij);
        auto Vee = compute_integral(b_set, &Vee_ij);
        auto Te = compute_integral(b_set, &T_ij);
        auto H = Vne + Vee + Te;
        auto e = solve_HC_SCE(H, S, false);
        std::cout << N << "\t" << b_set.size() <<"\t" << e << std::endl;
    }
    return 0;
}



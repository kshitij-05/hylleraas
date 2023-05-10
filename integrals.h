//
// Created by Kshitij Surjuse on 5/9/23.
//

#ifndef HYLLERAAS_INTEGRALS_H
#define HYLLERAAS_INTEGRALS_H

#include "basis.h"
#include "k.h"
#include <Eigen/Dense>

namespace hylleraas {

    typedef long double numeric_type;
    typedef Eigen::Matrix<numeric_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

    /// computes overlap integral between basis function bf1 anf bf2
    numeric_type S_ij(const basis &bf1, const basis &bf2) {
        return K(bf1.n_ + bf2.n_, bf1.l_ + bf2.l_, bf1.m_ + bf2.m_,
                 bf1.alpha_, bf1.beta_, bf1.gamma_);
    }

    /// computes nu-ele attraction integral between basis function bf1 anf bf2
    numeric_type Vne_ij(const basis &bf1, const basis &bf2) {
        numeric_type one = 1.0;
        auto k1 = K(bf1.n_ + bf2.n_ - one, bf1.l_ + bf2.l_, bf1.m_ + bf2.m_,
                    bf1.alpha_, bf1.beta_, bf1.gamma_);
        auto k2 = K(bf1.n_ + bf2.n_, bf1.l_ + bf2.l_ - one, bf1.m_ + bf2.m_,
                    bf1.alpha_, bf1.beta_, bf1.gamma_);
        return -one * bf1.Z_ * (k1 + k2);
    }

    /// computes ele-ele repulsion integral between basis function bf1 anf bf2
    numeric_type Vee_ij(const basis &bf1, const basis &bf2) {
        numeric_type one = 1.0;
        return K(bf1.n_ + bf2.n_, bf1.l_ + bf2.l_, bf1.m_ + bf2.m_ - one,
                 bf1.alpha_, bf1.beta_, bf1.gamma_);
    }

    /// computes kinetic energy repulsion integral between basis function bf1 anf bf2
    numeric_type T_ij(const basis &bf1, const basis &bf2) {

        auto alpha = bf1.alpha_;
        auto beta = bf1.beta_;
        auto gamma = bf1.gamma_;
        auto nj = bf2.n_;
        auto lj = bf2.l_;
        auto mj = bf2.m_;

        auto screen_prefactor = [](numeric_type prefactor) {
            return prefactor !=0.0;
        };

        auto Knlm = [&](numeric_type prefactor, numeric_type n, numeric_type l, numeric_type m) {
            if(screen_prefactor(prefactor)){
              return prefactor * K(bf1.n_ + bf2.n_ + n, bf1.l_ + bf2.l_ + l, bf1.m_ + bf2.m_ + m, alpha,
                              beta, gamma);
            }
            else{
                numeric_type zero = 0.0;
                return zero;
            }
        };

        numeric_type Tij = 0.;

        numeric_type zero, one, two, half,one_fourth,one_eighth ;
        zero = 0.;
        one = 1.;
        two = 2.;
        half = 0.5;
        one_fourth = 0.25;
        one_eighth = 0.125;

        Tij = -one_eighth* (pow(alpha, two) + pow(beta, two) + pow(gamma, two)) * S_ij(bf1, bf2);

        Tij += Knlm((half * nj * alpha), -one, zero, zero);
        Tij += Knlm((half * lj * beta) ,zero, -one, zero);
        Tij += Knlm((mj * gamma),zero, zero, -one);

        Tij -= Knlm((half * nj * (nj - one)),-two, zero, zero);
        Tij -= Knlm((half * lj * (lj - one)),zero, -two, zero);
        Tij -= Knlm((mj * (mj - one)),zero, zero, -two);

        Tij += Knlm((half * alpha),-one, zero, zero);
        Tij += Knlm((half * alpha),zero, -one, zero);
        Tij += Knlm(gamma, zero, zero, -one);


        Tij -= Knlm(nj,-two, zero,zero);
        Tij -= Knlm(lj, zero, -two, zero);
        Tij -= Knlm((two * mj),zero, zero, -two);

        Tij -= Knlm((one_eighth * alpha * gamma),-one, zero, one)
                + Knlm((one_eighth * alpha * gamma),one, zero, -one)
                - Knlm((one_eighth * alpha * gamma),-one, two, -one);

        Tij -= Knlm((one_eighth * beta * gamma),zero, -one, one)
                + Knlm((one_eighth * beta * gamma),zero, one, -one)
                - Knlm((one_eighth * beta * gamma),two, -one, -one);

        Tij += Knlm((one_fourth * nj * gamma),-two, zero, one)
                + Knlm((one_fourth * nj * gamma),zero, zero, -one)
                - Knlm((one_fourth * nj * gamma),-two, two, -one);

        Tij += Knlm((one_fourth * mj * alpha),-one, zero, zero)
                + Knlm((one_fourth * mj * alpha),one, zero, -two)
                - Knlm((one_fourth * mj * alpha),-one, two, -two);

        Tij -= Knlm((half * nj * mj),-two, zero, zero)
                + Knlm((half * nj * mj),zero, zero, -two)
                - Knlm((half * nj * mj),-two, two, -two);

        Tij += Knlm((one_fourth * lj * gamma),zero, -two, one)
                + Knlm((one_fourth * lj * gamma),zero, zero, -one)
                - Knlm((one_fourth * lj * gamma),two, -two, -one);

        Tij += Knlm((one_fourth * mj * beta),zero, -one , zero)
                + Knlm((one_fourth * mj * beta),zero, one, -two)
                - Knlm((one_fourth * mj * beta),two,-one, -two);

        Tij -= Knlm((half * lj * mj),zero, -two, zero)
                + Knlm((half * lj * mj),zero, zero, -two)
                - Knlm((half * lj * mj),two, -two, -two);

        return Tij;
    }


    /// computes integral matrix given a vector of basis functions
    /// @param2 a reference to the function that computes the elements of the matrix
    template<typename integral_type>
    Matrix compute_integral(const std::vector<basis> &bfn, const integral_type &integral_ij) {
        auto n = bfn.size();
        Matrix matrix(n, n);
        numeric_type zero = 0.0;
        matrix.fill(zero);
        for (auto i = 0; i < n; ++i) {
            for (auto j = 0; j < n; ++j) {
                matrix(i, j) = integral_ij(bfn[i], bfn[j]);
            }
        }
        return matrix;
    }

}

#endif //HYLLERAAS_INTEGRALS_H

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
            return (abs(prefactor) !=0.0);
        };

        auto Knlm = [&](numeric_type n, numeric_type l, numeric_type m) {
            return K(bf1.n_ + bf2.n_ + n, bf1.l_ + bf2.l_ + l, bf1.m_ + bf2.m_ + m, alpha,
                     beta, gamma);
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

        Tij += screen_prefactor(half * nj * alpha) ? half * nj * alpha * Knlm(-one, zero, zero) : zero;
        Tij += screen_prefactor(half * lj * beta) ? half * lj * beta * Knlm(zero, -one, zero) : zero;
        Tij += screen_prefactor(mj * gamma) ? mj * gamma * Knlm(zero, zero, -one) : zero;

        Tij -= screen_prefactor(half * nj * (nj - one)) ? half * nj * (nj - one) * Knlm(-two, zero, zero) : zero;
        Tij -= screen_prefactor(half * lj * (lj - one)) ? half * lj * (lj - one) * Knlm(zero, -two, zero) : zero;
        Tij -= screen_prefactor(mj * (mj - one)) ? mj * (mj - one) * Knlm(zero, zero, -two) : zero;

        Tij += screen_prefactor(half * alpha) ? half * alpha * Knlm(-one, zero, zero) : zero;
        Tij += screen_prefactor(half * alpha) ? half * beta * Knlm(zero, -one, zero) : zero;
        Tij += screen_prefactor(gamma) ? gamma * Knlm(zero, zero, -one) : zero;


        Tij -= screen_prefactor(nj) ? nj * Knlm(-two, zero,zero) : zero;
        Tij -= screen_prefactor(lj) ? lj * Knlm(zero, -two, zero) : zero;
        Tij -= screen_prefactor(two * mj) ? two * mj * Knlm(zero, zero, -two) : zero;

        Tij -= screen_prefactor(one_eighth * alpha * gamma) ? one_eighth * alpha * gamma *
                (Knlm(-one, zero, one) + Knlm(one, zero, -one) - Knlm(-one, two, -one)) : zero;

        Tij -= screen_prefactor(one_eighth * beta * gamma) ? one_eighth * beta * gamma *
                (Knlm(zero, -one, one) + Knlm(zero, one, -one) - Knlm(two, -one, -one)) : zero;

        Tij += screen_prefactor(one_fourth * nj * gamma) ? one_fourth * nj * gamma *
                (Knlm(-two, zero, one) + Knlm(zero, zero, -one) - Knlm(-two, two, -one)) : zero;

        Tij += screen_prefactor(one_fourth * mj * alpha) ? one_fourth * mj * alpha *
                (Knlm(-one, zero, zero) + Knlm(one, zero, -two) - Knlm(-one, two, -two)): zero;

        Tij -= screen_prefactor(half * nj * mj) ? half * nj * mj *
                (Knlm(-two, zero, zero) + Knlm(zero, zero, -two) - Knlm(-two, two, -two)): zero;

        Tij += screen_prefactor(one_fourth * lj * gamma) ? one_fourth * lj * gamma *
                (Knlm(zero, -two, one) + Knlm(zero, zero, -one) - Knlm(two, -two, -one)): zero;

        Tij += screen_prefactor(one_fourth * mj * beta) ? one_fourth * mj * beta *
                (Knlm(zero, -one , zero) + Knlm(zero, one, -two) - Knlm(two,-one, -two)): zero;

        Tij -= screen_prefactor(half * lj * mj) ? half * lj * mj *
                (Knlm(zero, two, zero) + Knlm(zero, zero, -two) - Knlm(two, -two, -two)): zero;

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

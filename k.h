//
// Created by Kshitij Surjuse on 5/9/23.
//

#ifndef HYLLERAAS_K_H
#define HYLLERAAS_K_H

#include "basis.h"
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace boost::math;

namespace hylleraas {

    template<typename intT, typename doubleT>
    doubleT K(const intT &n, const intT &l, const intT &m,
              const doubleT &alpha, const doubleT &beta, const doubleT &gamma) {

        intT one = 1.;
        intT two = 2.;
        doubleT sixteen = 16.;
        const double pi = boost::math::constants::pi<doubleT>();
        auto prefactor = sixteen * pow(pi, two) * factorial<doubleT>(n + one)
                         * factorial<doubleT>(l + one) * factorial<doubleT>(m + one);
        doubleT k = 0.0;
        for (auto a = 0; a <= n + 1; ++a) {
            for (auto b = 0; b <= l + 1; ++b) {
                for (auto c = 0; c <= m + 1; ++c) {
                    auto numerator = binomial_coefficient<doubleT>(l + one - b + a, a)
                                     * binomial_coefficient<doubleT>(m + one - c + b, b)
                                     * binomial_coefficient<doubleT>(n + one - a + c, c);
                    auto denominator = pow(alpha + beta, l - b + a + two)
                                       * pow(alpha + gamma, n - a + c + two)
                                       * pow(beta + gamma, m - c + b + two);
                    k += numerator / denominator;
                }
            }
        }
        return prefactor * k;
    }

}

#endif //HYLLERAAS_K_H

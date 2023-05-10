//
// Created by Kshitij Surjuse on 5/9/23.
//

#ifndef HYLLERAAS_BASIS_H
#define HYLLERAAS_BASIS_H

namespace hylleraas {

    typedef long double numeric_type;

    struct basis {
        numeric_type Z_ = 2.;  // we will be only studying He atom
        numeric_type n_, l_, m_;
        numeric_type alpha_, beta_, gamma_;

        basis(numeric_type n, numeric_type l, numeric_type m, numeric_type alpha,
              numeric_type beta, numeric_type gamma) {
            n_ = n;
            l_ = l;
            m_ = m;
            alpha_ = alpha;
            beta_ = beta;
            gamma_ = gamma;
        }
    };

    inline std::ostream &operator<<(std::ostream &os, const basis &bf) {
        os << "{ n=" << bf.n_ << ", l= " << bf.l_ << ", m=" << bf.m_ << ", alpha = "
           << bf.alpha_ << ", beta= " << bf.beta_ << ", gamma=" << bf.gamma_ << "}" << std::endl;
        return os;
    }

    std::vector<basis> test_basis() {
        numeric_type alpha = 2. * 1.8;
        numeric_type beta = 2. * 1.8;
        numeric_type gamma = 0.0;
        std::vector<basis> test_basis;
        test_basis.emplace_back(0, 0, 0, alpha, beta, gamma);
        test_basis.emplace_back(1, 1, 0, alpha, beta, gamma);
        test_basis.emplace_back(0, 0, 1, alpha, beta, gamma);

        return test_basis;
    }

    std::vector<basis> lambda_N(numeric_type N, numeric_type alpha, numeric_type gamma) {
        std::vector<basis> lambda_n;
        for (auto n = 0; n <= N; ++n) {
            for (auto l = 0; l <= N - n; ++l) {
                for (auto m = 0; m <= N - l - n; ++m) {
                    basis bs = basis(n, l, m, alpha, alpha, gamma);
                    lambda_n.emplace_back(bs);
                }
            }
        }
        return lambda_n;
    }
}

#endif //HYLLERAAS_BASIS_H

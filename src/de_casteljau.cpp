// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "eft.hpp"
#include <tuple>
#include <vector>

namespace de_casteljau {

double basic(double s, const std::vector<double>& coeffs)
{
    std::vector<double> pk = coeffs;

    double r = 1.0 - s;
    size_t degree = coeffs.size() - 1;
    for (size_t k = 0; k < degree; ++k) {
        for (size_t j = 0; j < degree - k; ++j) {
            pk[j] = r * pk[j] + s * pk[j + 1];
        }
    }

    // NOTE: This **assumes** ``degree >= 0``.
    return pk[0];
}

double local_error(
    const std::vector<double>& errors, double rho, double delta_b)
{
    double l_hat = 0.0;
    for (double const& error : errors) {
        l_hat += error;
    }
    l_hat += rho * delta_b;

    return l_hat;
}

double local_error_eft(std::vector<double>& errors, double rho, double delta_b)
{
    size_t num_errs = errors.size();
    errors.resize(num_errs + 1);

    // NOTE: This **assumes** ``num_errs >= 2``.
    double l_hat;
    std::tie(l_hat, errors[0]) = eft::two_sum(errors[0], errors[1]);
    for (size_t j = 2; j < num_errs; ++j) {
        std::tie(l_hat, errors[j - 1]) = eft::two_sum(l_hat, errors[j]);
    }

    double prod;
    std::tie(prod, errors[num_errs - 1]) = eft::two_prod(rho, delta_b);
    std::tie(l_hat, errors[num_errs]) = eft::two_sum(l_hat, prod);

    return l_hat;
}

std::vector<double> compensated_k(
    double s, const std::vector<double>& coeffs, size_t K)
{
    // NOTE: This function **assumes** ``K >= 2``.
    double r, rho;
    std::tie(r, rho) = eft::two_sum(1.0, -s);

    std::vector<double> errors(5 * K - 7);
    std::vector<std::vector<double>> bk(K);
    // NOTE: This **assumes** ``degree >= 0``.
    size_t degree = coeffs.size() - 1;

    // Initialize ``bk``.
    bk[0] = coeffs;
    for (size_t F = 1; F < K; ++F) {
        bk[F].resize(degree + 1, 0.0);
    }

    for (size_t k = 0; k < degree; ++k) {
        for (size_t j = 0; j < degree - k; ++j) {
            // NOTE: The ``size()`` of ``errors`` is important since it is
            //       used by ``local_error_eft()``.
            errors.resize(3);
            double delta_b = bk[0][j];

            // Update the "level 0" stuff.
            double val1, val2;
            std::tie(val1, errors[0]) = eft::two_prod(r, bk[0][j]);
            std::tie(val2, errors[1]) = eft::two_prod(s, bk[0][j + 1]);
            std::tie(bk[0][j], errors[2]) = eft::two_sum(val1, val2);

            for (size_t F = 1; F < K - 1; ++F) {
                double error;
                val1 = local_error_eft(errors, rho, delta_b);
                delta_b = bk[F][j];

                std::tie(val2, error) = eft::two_prod(s, bk[F][j + 1]);
                errors.push_back(error);

                std::tie(val1, error) = eft::two_sum(val1, val2);
                errors.push_back(error);

                std::tie(val2, error) = eft::two_prod(r, bk[F][j]);
                errors.push_back(error);

                std::tie(bk[F][j], error) = eft::two_sum(val1, val2);
                errors.push_back(error);
            }

            // Update the "level 2" stuff.
            val1 = local_error(errors, rho, delta_b);
            bk[K - 1][j] = val1 + s * bk[K - 1][j + 1] + r * bk[K - 1][j];
        }
    }

    std::vector<double> b_hat(K);
    for (size_t j = 0; j < K; ++j) {
        b_hat[j] = bk[j][0];
    }
    return b_hat;
}

double compensated2(double s, const std::vector<double>& coeffs)
{
    std::vector<double> terms = compensated_k(s, coeffs, 2);
    return eft::sum_k(terms, 2);
}

double compensated3(double s, const std::vector<double>& coeffs)
{
    std::vector<double> terms = compensated_k(s, coeffs, 3);
    return eft::sum_k(terms, 3);
}

double compensated4(double s, const std::vector<double>& coeffs)
{
    std::vector<double> terms = compensated_k(s, coeffs, 4);
    return eft::sum_k(terms, 4);
}
}

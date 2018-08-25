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

double local_error(double* errors, size_t num_errs, double rho, double delta_b)
{
    // NOTE: This **assumes** ``errors`` has (at least) ``num_errs`` entries
    //       and that ``num_errs >= 1``.
    double l_hat = errors[0];
    for (size_t j = 1; j < num_errs; ++j) {
        l_hat += errors[j];
    }
    l_hat += rho * delta_b;

    return l_hat;
}

double local_error_eft(
    double* errors, size_t num_errs, double rho, double delta_b)
{
    // NOTE: This **assumes** ``errors`` has (at least) ``num_errs + 1``
    //       entries and that ``num_errs >= 2``.
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
    // NOTE: This **assumes** ``degree >= 0``.
    size_t degree = coeffs.size() - 1;
    std::vector<double> bk((degree + 1) * K, 0.0);

    // Initialize ``bk``.
    for (size_t j = 0; j <= degree; ++j) {
        bk[j] = coeffs[j];
    }

    for (size_t k = 0; k < degree; ++k) {
        for (size_t j = 0; j < degree - k; ++j) {
            double delta_b = bk[j];

            // Update the "level 0" stuff.
            double val1, val2;
            std::tie(val1, errors[0]) = eft::two_prod(r, bk[j]);
            std::tie(val2, errors[1]) = eft::two_prod(s, bk[j + 1]);
            std::tie(bk[j], errors[2]) = eft::two_sum(val1, val2);

            size_t num_errs = 3;
            size_t index_shift = degree + 1;
            for (size_t F = 1; F < K - 1; ++F) {
                val1 = local_error_eft(&errors[0], num_errs, rho, delta_b);
                delta_b = bk[index_shift + j];

                std::tie(val2, errors[num_errs + 1])
                    = eft::two_prod(s, bk[index_shift + j + 1]);

                std::tie(val1, errors[num_errs + 2])
                    = eft::two_sum(val1, val2);

                std::tie(val2, errors[num_errs + 3])
                    = eft::two_prod(r, bk[index_shift + j]);

                std::tie(bk[index_shift + j], errors[num_errs + 4])
                    = eft::two_sum(val1, val2);

                num_errs += 5;
                // Update the index shift for the next iteration.
                index_shift += degree + 1;
            }

            // Update the "level 2" stuff.
            val1 = local_error(&errors[0], num_errs, rho, delta_b);
            bk[index_shift + j]
                = val1 + s * bk[index_shift + j + 1] + r * bk[index_shift + j];
        }
    }

    std::vector<double> b_hat(K);
    for (size_t F = 0; F < K; ++F) {
        b_hat[F] = bk[(degree + 1) * F];
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

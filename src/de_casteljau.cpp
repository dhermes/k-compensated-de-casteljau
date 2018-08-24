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

#include <cmath>
#include <tuple>
#include <utility>
#include <vector>

namespace eft {

std::pair<double, double> two_prod(double val1, double val2)
{
    double product, error;
    product = val1 * val2;
    error = std::fma(val1, val2, -product);
    return std::make_pair(product, error);
}

std::pair<double, double> two_sum(double val1, double val2)
{
    double sum, almost_val2, error;
    sum = val1 + val2;
    almost_val2 = sum - val1;
    error = (val1 - (sum - almost_val2)) + (val2 - almost_val2);
    return std::make_pair(sum, error);
}
}

namespace de_casteljau {

double basic(double s, const std::vector<double>& coeffs)
{
    size_t degree, j, k;
    double r;
    std::vector<double> pk;

    pk = coeffs;
    r = 1.0 - s;
    degree = coeffs.size() - 1;
    for (k = 0; k < degree; ++k) {
        for (j = 0; j < degree - k; ++j) {
            pk[j] = r * pk[j] + s * pk[j + 1];
        }
    }

    // NOTE: This **assumes** ``degree >= 0``.
    return pk[0];
}

double local_error(
    const std::vector<double>& errors, double rho, double delta_b)
{
    double l_hat;

    l_hat = 0.0;
    for (double const& error : errors) {
        l_hat += error;
    }
    l_hat += rho * delta_b;

    return l_hat;
}

double local_error_eft(std::vector<double>& errors, double rho, double delta_b)
{
    size_t num_errs, j;
    double l_hat, prod;

    num_errs = errors.size();
    errors.resize(num_errs + 1);

    // NOTE: This **assumes** ``num_errs >= 2``.
    std::tie(l_hat, errors[0]) = eft::two_sum(errors[0], errors[1]);
    for (j = 2; j < num_errs; ++j) {
        std::tie(l_hat, errors[j - 1]) = eft::two_sum(l_hat, errors[j]);
    }

    std::tie(prod, errors[num_errs - 1]) = eft::two_prod(rho, delta_b);
    std::tie(l_hat, errors[num_errs]) = eft::two_sum(l_hat, prod);

    return l_hat;
}

std::vector<double> compensated_k(
    double s, const std::vector<double>& coeffs, size_t K)
{
    // NOTE: This function **assumes** ``K >= 2``.
    size_t degree, F, j, k;
    std::vector<double> b_hat, errors;
    std::vector<std::vector<double>> bk;
    double r, rho, val1, val2, error, delta_b;

    std::tie(r, rho) = eft::two_sum(1.0, -s);

    errors.reserve(5 * K - 7);
    bk.resize(K);
    // NOTE: This **assumes** ``degree >= 0``.
    degree = coeffs.size() - 1;

    // Initialize ``bk``.
    bk[0] = coeffs;
    for (F = 1; F < K; ++F) {
        bk[F].resize(degree + 1, 0.0);
    }

    for (k = 0; k < degree; ++k) {
        for (j = 0; j < degree - k; ++j) {
            // NOTE: The ``size()`` of ``errors`` is important since it is
            //       used by ``local_error_eft()``.
            errors.resize(3);
            delta_b = bk[0][j];

            // Update the "level 0" stuff.
            std::tie(val1, errors[0]) = eft::two_prod(r, bk[0][j]);
            std::tie(val2, errors[1]) = eft::two_prod(s, bk[0][j + 1]);
            std::tie(bk[0][j], errors[2]) = eft::two_sum(val1, val2);

            for (F = 1; F < K - 1; ++F) {
                val1 = de_casteljau::local_error_eft(errors, rho, delta_b);
                delta_b = bk[F][j];

                std::tie(val2, error) = eft::two_prod(s, bk[F][j + 1]);
                errors.push_back(error);

                std::tie(val1, error) = eft::two_sum(val1, val2);
                errors.push_back(error);

                std::tie(val2, error) = eft::two_prod(r, bk[F][j]);
                errors.push_back(error);

                std::tie(bk[F][j], error) = eft::two_prod(val1, val2);
                errors.push_back(error);
            }

            // Update the "level 2" stuff.
            val1 = de_casteljau::local_error(errors, rho, delta_b);
            bk[K - 1][j] = val1 + s * bk[K - 1][j + 1] + r * bk[K - 1][j];
        }
    }

    b_hat.resize(K);
    for (j = 0; j < K; ++j) {
        b_hat[j] = bk[j][0];
    }
    return b_hat;
}
}

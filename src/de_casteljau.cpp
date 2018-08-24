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

double basic(double s, std::vector<double>& coeffs)
{
    size_t degree, j, k;
    double r;
    std::vector<double> pk(coeffs);

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

double local_error(std::vector<double>& errors, double rho, double delta_b)
{
    double l_hat;

    l_hat = 0.0;
    for (double const& error : errors) {
        l_hat += error;
    }
    l_hat += rho * delta_b;

    return l_hat;
}

std::pair<std::vector<double>, double> local_error_eft(
    std::vector<double>& errors, double rho, double delta_b)
{
    size_t num_errs, j;
    std::vector<double> new_errors;
    double l_hat, prod;

    num_errs = errors.size();
    new_errors.reserve(num_errs + 1);

    // NOTE: This **assumes** ``num_errs >= 2``.
    std::tie(l_hat, new_errors[0]) = eft::two_sum(errors[0], errors[1]);
    for (j = 2; j < num_errs; ++j) {
        std::tie(l_hat, new_errors[j - 1]) = eft::two_sum(l_hat, errors[j]);
    }

    std::tie(prod, new_errors[num_errs - 1]) = eft::two_prod(rho, delta_b);
    std::tie(l_hat, new_errors[num_errs]) = eft::two_sum(l_hat, prod);

    return std::make_pair(new_errors, l_hat);
}
}

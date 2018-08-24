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

    return pk[0];
}
}

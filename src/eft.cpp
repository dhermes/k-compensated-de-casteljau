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
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

namespace eft {

std::pair<double, double> two_prod(double val1, double val2)
{
    double product = val1 * val2;
    double error = std::fma(val1, val2, -product);
    return std::make_pair(product, error);
}

std::pair<double, double> two_sum(double val1, double val2)
{
    double sum = val1 + val2;
    double almost_val2 = sum - val1;
    double error = (val1 - (sum - almost_val2)) + (val2 - almost_val2);
    return std::make_pair(sum, error);
}

void vec_sum(std::vector<double>& vec)
{
    for (size_t i = 1; i < vec.size(); ++i) {
        std::tie(vec[i], vec[i - 1]) = two_sum(vec[i], vec[i - 1]);
    }
}

double sum_k(const std::vector<double>& vec, size_t K)
{
    std::vector<double> workspace = vec;

    for (size_t i = 0; i < K - 1; ++i) {
        vec_sum(workspace);
    }

    double result = std::accumulate(workspace.begin(), workspace.end(), 0.0);
    return result;
}
}

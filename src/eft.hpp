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

#ifndef EFT_HPP
#define EFT_HPP

#include <utility>
#include <vector>

namespace eft {
std::pair<double, double> two_prod(double val1, double val2);
std::pair<double, double> two_sum(double val1, double val2);
double sum_k(const std::vector<double>& vec, size_t K);
}

#endif /* EFT_HPP */

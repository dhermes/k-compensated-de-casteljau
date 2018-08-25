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

#include <vector>

namespace de_casteljau {
double basic(double s, const std::vector<double>& coeffs);
std::vector<double> compensated_k(
    double s, const std::vector<double>& coeffs, size_t K);
double compensated2(double s, const std::vector<double>& coeffs);
double compensated3(double s, const std::vector<double>& coeffs);
double compensated4(double s, const std::vector<double>& coeffs);
}

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

#include "de_casteljau.hpp"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

int main()
{
    int i;
    double s, evaluated;
    std::vector<double> coeffs1, coeffs2, coeffs3;
    std::pair<double, double> eft_result;

    coeffs1 = { 3.0, 2.0 };
    std::cout << "p1(s) = 3(1 - s) + 2s" << std::endl;
    coeffs2 = { 3.0, 2.5, 0.0 };
    std::cout << "p2(s) = 3(1 - s)^2 + 2.5[2(1 - s)s]" << std::endl;
    coeffs3 = { 4.0, 0.0, 0.0, -1.0 };
    std::cout << "p3(s) = 4(1 - s)^3 - s^3" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    for (i = -2; i < 3; ++i) {
        s = 1.5 * i;
        evaluated = de_casteljau::basic(s, coeffs1);
        std::cout << "p1(" << s << ") = " << evaluated << std::endl;
        evaluated = de_casteljau::basic(s, coeffs2);
        std::cout << "p2(" << s << ") = " << evaluated << std::endl;
        evaluated = de_casteljau::basic(s, coeffs3);
        std::cout << "p3(" << s << ") = " << evaluated << std::endl;
        if (i < 2) {
            std::cout << std::string(30, '*') << std::endl;
        }
    }

    std::cout << std::string(60, '=') << std::endl;
    eft_result = eft::two_prod(1.0 - pow(0.5, 27), 1.0 + pow(0.5, 27));
    std::cout << "two_prod(1 - 2^{-27}, 1 + 2^{-27}) = P + pi, where"
              << std::endl;
    std::cout << "P  = " << eft_result.first << std::endl;
    std::cout << "pi = " << eft_result.second * pow(2.0, 54) << " / 2^{54}"
              << std::endl;

    eft_result = eft::two_sum(1.0 + pow(0.5, 52), pow(0.5, 53));
    std::cout << "two_sum(1 + 2^{-52}, 2^{-53}) = S + sigma, where"
              << std::endl;
    std::cout << "S     = 1 + " << (eft_result.first - 1.0) * pow(2.0, 51)
              << " / 2^{51}" << std::endl;
    std::cout << "sigma = " << eft_result.second * pow(2.0, 53) << " / 2^{53}"
              << std::endl;
    return 0;
}
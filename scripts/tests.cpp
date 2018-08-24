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
    size_t i, j, F, K;
    double s, evaluated;
    std::vector<double> compensated;
    std::vector<std::vector<double>> coeffs;
    std::pair<double, double> eft_result;

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

    std::cout << std::string(60, '=') << std::endl;
    coeffs.resize(3);
    coeffs[0] = { 3.0, 2.0 };
    std::cout << "p1(s) = 3(1 - s) + 2s" << std::endl;
    coeffs[1] = { 3.0, 2.5, 0.0 };
    std::cout << "p2(s) = 3(1 - s)^2 + 2.5[2(1 - s)s]" << std::endl;
    coeffs[2] = { 4.0, 0.0, 0.0, -1.0 };
    std::cout << "p3(s) = 4(1 - s)^3 - s^3" << std::endl;

    std::cout << std::string(60, '=') << std::endl;
    std::cout << std::scientific;
    std::cout << "DeCasteljau:" << std::endl;
    for (i = 0; i < 5; ++i) {
        s = 1.5 * i - 3.0;
        for (j = 0; j < 3; ++j) {
            evaluated = de_casteljau::basic(s, coeffs[j]);
            std::cout << "p" << j + 1 << "(" << s << ") = " << evaluated
                      << std::endl;
        }
        std::cout << std::string(30, '*') << std::endl;
    }

    for (K = 2; K <= 4; ++K) {
        std::cout << std::string(60, '=') << std::endl;
        std::cout << "CompDeCasteljau" << K << ":" << std::endl;
        for (i = 0; i < 5; ++i) {
            s = 1.5 * i - 3.0 + pow(0.5, 50);
            for (j = 0; j < 3; ++j) {
                compensated = de_casteljau::compensated_k(s, coeffs[j], K);
                std::cout << "p" << j + 1 << "(" << s
                          << ") = " << compensated[0];
                for (F = 1; F < K; ++F) {
                    std::cout << " + (" << compensated[F] << ")";
                }
                std::cout << std::endl;
            }
            std::cout << std::string(30, '*') << std::endl;
        }
    }

    return 0;
}

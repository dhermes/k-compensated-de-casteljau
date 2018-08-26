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
#include "eft.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

int main()
{
    std::pair<double, double> eft_result
        = eft::two_prod(1.0 - pow(0.5, 27), 1.0 + pow(0.5, 27));
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
    std::vector<std::vector<double>> coeffs(3);
    coeffs[0] = { 3.0, 2.0 };
    std::array<double, 2> coeffs1 = { 3.0, 2.0 };
    std::cout << "p1(s) = 3(1 - s) + 2s" << std::endl;
    coeffs[1] = { 3.0, 2.5, 0.0 };
    std::array<double, 3> coeffs2 = { 3.0, 2.5, 0.0 };
    std::cout << "p2(s) = 3(1 - s)^2 + 2.5[2(1 - s)s]" << std::endl;
    coeffs[2] = { 4.0, 0.0, 0.0, -1.0 };
    std::array<double, 4> coeffs3 = { 4.0, 0.0, 0.0, -1.0 };
    std::cout << "p3(s) = 4(1 - s)^3 - s^3" << std::endl;

    std::cout << std::string(60, '=') << std::endl;
    std::cout << std::scientific;
    std::cout << "DeCasteljau:" << std::endl;
    for (size_t i = 0; i < 5; ++i) {
        double s = 1.5 * i - 3.0;
        for (size_t j = 0; j < 3; ++j) {
            double evaluated = de_casteljau::basic(s, coeffs[j]);
            std::cout << "p" << j + 1 << "(" << s << ") = " << evaluated
                      << std::endl;
        }
        std::cout << std::string(30, '*') << std::endl;
    }

    std::cout << std::string(60, '=') << std::endl;
    std::cout << "CompDeCasteljau:" << std::endl;
    for (size_t i = 0; i < 5; ++i) {
        double s = 1.5 * i - 3.0 + pow(0.5, 50);

        std::array<double, 2> compensated2
            = de_casteljau::compensated<2, 1>(s, coeffs1);
        std::cout << "p1(" << s << ") = " << compensated2[0] << " + ("
                  << compensated2[1] << ")" << std::endl;

        compensated2 = de_casteljau::compensated<2, 2>(s, coeffs2);
        std::cout << "p2(" << s << ") = " << compensated2[0] << " + ("
                  << compensated2[1] << ")" << std::endl;

        compensated2 = de_casteljau::compensated<2, 3>(s, coeffs3);
        std::cout << "p3(" << s << ") = " << compensated2[0] << " + ("
                  << compensated2[1] << ")" << std::endl;

        std::cout << std::string(30, '*') << std::endl;
    }

    std::cout << std::string(60, '=') << std::endl;
    std::cout << "CompDeCasteljau3:" << std::endl;
    for (size_t i = 0; i < 5; ++i) {
        double s = 1.5 * i - 3.0 + pow(0.5, 50);

        std::array<double, 3> compensated3
            = de_casteljau::compensated<3, 1>(s, coeffs1);
        std::cout << "p1(" << s << ") = " << compensated3[0] << " + ("
                  << compensated3[1] << ") + (" << compensated3[2] << ")"
                  << std::endl;

        compensated3 = de_casteljau::compensated<3, 2>(s, coeffs2);
        std::cout << "p2(" << s << ") = " << compensated3[0] << " + ("
                  << compensated3[1] << ") + (" << compensated3[2] << ")"
                  << std::endl;

        compensated3 = de_casteljau::compensated<3, 3>(s, coeffs3);
        std::cout << "p3(" << s << ") = " << compensated3[0] << " + ("
                  << compensated3[1] << ") + (" << compensated3[2] << ")"
                  << std::endl;

        std::cout << std::string(30, '*') << std::endl;
    }

    std::cout << std::string(60, '=') << std::endl;
    std::cout << "CompDeCasteljau4:" << std::endl;
    for (size_t i = 0; i < 5; ++i) {
        double s = 1.5 * i - 3.0 + pow(0.5, 50);

        std::array<double, 4> compensated4
            = de_casteljau::compensated<4, 1>(s, coeffs1);
        std::cout << "p1(" << s << ") = " << compensated4[0] << " + ("
                  << compensated4[1] << ") + (" << compensated4[2] << ") + ("
                  << compensated4[3] << ")" << std::endl;

        compensated4 = de_casteljau::compensated<4, 2>(s, coeffs2);
        std::cout << "p2(" << s << ") = " << compensated4[0] << " + ("
                  << compensated4[1] << ") + (" << compensated4[2] << ") + ("
                  << compensated4[3] << ")" << std::endl;

        compensated4 = de_casteljau::compensated<4, 3>(s, coeffs3);
        std::cout << "p3(" << s << ") = " << compensated4[0] << " + ("
                  << compensated4[1] << ") + (" << compensated4[2] << ") + ("
                  << compensated4[3] << ")" << std::endl;

        std::cout << std::string(30, '*') << std::endl;
    }

    return 0;
}

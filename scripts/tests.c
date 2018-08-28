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

#include "de_casteljau.h"
#include "eft.h"
#include <math.h>
#include <stdio.h>

void print_sep(void)
{
    printf("============================================================\n");
}

void print_small_sep(void) { printf("******************************\n"); }

int main(void)
{
    double product, error;
    two_prod(1.0 - pow(0.5, 27), 1.0 + pow(0.5, 27), &product, &error);
    printf("two_prod(1 - 2^{-27}, 1 + 2^{-27}) = P + pi, where\n");
    printf("P  = %f\n", product);
    printf("pi = %f / 2^{54}\n", error * pow(2.0, 54));

    double sum;
    two_sum(1.0 + pow(0.5, 52), pow(0.5, 53), &sum, &error);
    printf("two_sum(1 + 2^{-52}, 2^{-53}) = S + sigma, where\n");
    printf("S     = 1 + %f / 2^{51}\n", (sum - 1.0) * pow(2.0, 51));
    printf("sigma = %f / 2^{53}\n", error * pow(2.0, 53));

    print_sep();
    double coeffs1[2] = { 3.0, 2.0 };
    printf("p1(s) = 3(1 - s) + 2s\n");
    double coeffs2[3] = { 3.0, 2.5, 0.0 };
    printf("p2(s) = 3(1 - s)^2 + 2.5[2(1 - s)s]\n");
    double coeffs3[4] = { 4.0, 0.0, 0.0, -1.0 };
    printf("p3(s) = 4(1 - s)^3 - s^3\n");

    print_sep();
    printf("DeCasteljau:\n");
    double pk[4];
    for (size_t i = 0; i < 5; ++i) {
        double s = 1.5 * i - 3.0;
        double evaluated = basic(s, coeffs1, pk, 1);
        printf("p1(%e) = %e\n", s, evaluated);
        evaluated = basic(s, coeffs2, pk, 2);
        printf("p2(%e) = %e\n", s, evaluated);
        evaluated = basic(s, coeffs3, pk, 3);
        printf("p3(%e) = %e\n", s, evaluated);
        print_small_sep();
    }

    double errors[13];  // length(errors) >= 5K - 7; K <= 4;
    double bk[16];  // length(bk) >= K(degree + 1); K <= 4; degree <= 3;
    double result[4];  // length(result) >= K; K <= 4

    print_sep();
    printf("CompDeCasteljau:\n");
    for (size_t i = 0; i < 5; ++i) {
        double s = 1.5 * i - 3.0 + pow(0.5, 50);
        compensated(s, coeffs1, 1, 2, errors, bk, result);
        printf("p1(%e) = %e + (%e)\n", s, result[0], result[1]);
        compensated(s, coeffs2, 2, 2, errors, bk, result);
        printf("p2(%e) = %e + (%e)\n", s, result[0], result[1]);
        compensated(s, coeffs3, 3, 2, errors, bk, result);
        printf("p3(%e) = %e + (%e)\n", s, result[0], result[1]);
        print_small_sep();
    }

    print_sep();
    printf("CompDeCasteljau3:\n");
    for (size_t i = 0; i < 5; ++i) {
        double s = 1.5 * i - 3.0 + pow(0.5, 50);
        compensated(s, coeffs1, 1, 3, errors, bk, result);
        printf(
            "p1(%e) = %e + (%e) + (%e)\n", s, result[0], result[1], result[2]);
        compensated(s, coeffs2, 2, 3, errors, bk, result);
        printf(
            "p2(%e) = %e + (%e) + (%e)\n", s, result[0], result[1], result[2]);
        compensated(s, coeffs3, 3, 3, errors, bk, result);
        printf(
            "p3(%e) = %e + (%e) + (%e)\n", s, result[0], result[1], result[2]);
        print_small_sep();
    }

    print_sep();
    printf("CompDeCasteljau4:\n");
    for (size_t i = 0; i < 5; ++i) {
        double s = 1.5 * i - 3.0 + pow(0.5, 50);
        compensated(s, coeffs1, 1, 4, errors, bk, result);
        printf("p1(%e) = %e + (%e) + (%e) + (%e)\n", s, result[0], result[1],
            result[2], result[3]);
        compensated(s, coeffs2, 2, 4, errors, bk, result);
        printf("p2(%e) = %e + (%e) + (%e) + (%e)\n", s, result[0], result[1],
            result[2], result[3]);
        compensated(s, coeffs3, 3, 4, errors, bk, result);
        printf("p3(%e) = %e + (%e) + (%e) + (%e)\n", s, result[0], result[1],
            result[2], result[3]);
        print_small_sep();
    }

    return 0;
}

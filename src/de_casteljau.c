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

#include "eft.h"
#include <stdio.h>

double basic(double s, const double* coeffs, double* pk, size_t degree)
{
    // NOTE: This assumes that ``length(coeffs) == degree + 1``.
    // NOTE: This requires ``length(pk) >= length(coeffs)``.
    for (size_t i = 0; i <= degree; ++i) {
        pk[i] = coeffs[i];
    }

    double r = 1.0 - s;
    for (size_t k = 0; k < degree; ++k) {
        for (size_t j = 0; j < degree - k; ++j) {
            pk[j] = r * pk[j] + s * pk[j + 1];
        }
    }

    // NOTE: This **assumes** ``degree >= 0``.
    return pk[0];
}

double local_error(double* errors, size_t num_errs, double rho, double delta_b)
{
    // NOTE: This **assumes** ``length(errors) >= num_errs`` and that
    //       ``num_errs >= 1``.
    double l_hat = errors[0];
    for (size_t j = 1; j < num_errs; ++j) {
        l_hat += errors[j];
    }
    l_hat += rho * delta_b;

    return l_hat;
}

double local_error_eft(
    double* errors, size_t num_errs, double rho, double delta_b)
{
    // NOTE: This **assumes** ``length(errors) >= num_errs + 1`` and that
    //       ``num_errs >= 2``.
    double l_hat, tmp;
    two_sum(errors[0], errors[1], &l_hat, &tmp);
    errors[0] = tmp;
    for (size_t j = 2; j < num_errs; ++j) {
        two_sum(l_hat, errors[j], &tmp, &errors[j - 1]);
        l_hat = tmp;
    }

    double prod;
    two_prod(rho, delta_b, &prod, &errors[num_errs - 1]);
    two_sum(l_hat, prod, &tmp, &errors[num_errs]);

    return tmp;
}

void compensated(double s, const double* coeffs, size_t degree, size_t K,
    double* errors, double* bk, double* result)
{
    // NOTE: This function **assumes** ``K >= 2`` and ``degree >= 0``.
    // NOTE: This requires ``length(errors) >= 5K - 7``.
    // NOTE: This requires ``length(bk) >= K(degree + 1)``.
    // NOTE: This requires ``length(result) >= K``.
    double r, rho;
    two_sum(1.0, -s, &r, &rho);

    // Initialize ``bk`` (everything after ``F = 0`` is zero).
    for (size_t j = 0; j <= degree; ++j) {
        bk[j] = coeffs[j];
    }
    for (size_t j = degree + 1; j < K * (degree + 1); ++j) {
        bk[j] = 0.0;
    }

    for (size_t k = 0; k < degree; ++k) {
        for (size_t j = 0; j < degree - k; ++j) {
            double delta_b = bk[j];

            // Update the "level 0" stuff.
            double val1, val2, val3;
            two_prod(r, bk[j], &val1, &errors[0]);
            two_prod(s, bk[j + 1], &val2, &errors[1]);
            two_sum(val1, val2, &bk[j], &errors[2]);

            size_t num_errs = 3;
            size_t index_shift = degree + 1;
            for (size_t F = 1; F < K - 1; ++F) {
                val1 = local_error_eft(&errors[0], num_errs, rho, delta_b);
                delta_b = bk[index_shift + j];

                two_prod(
                    s, bk[index_shift + j + 1], &val2, &errors[num_errs + 1]);
                two_sum(val1, val2, &val3, &errors[num_errs + 2]);
                val1 = val3;
                two_prod(r, bk[index_shift + j], &val2, &errors[num_errs + 3]);
                two_sum(
                    val1, val2, &bk[index_shift + j], &errors[num_errs + 4]);

                num_errs += 5;
                // Update the index shift for the next iteration.
                index_shift += degree + 1;
            }

            // Update the "level 2" stuff.
            val1 = local_error(&errors[0], num_errs, rho, delta_b);
            bk[index_shift + j]
                = val1 + s * bk[index_shift + j + 1] + r * bk[index_shift + j];
        }
    }

    for (size_t F = 0; F < K; ++F) {
        result[F] = bk[(degree + 1) * F];
    }
}

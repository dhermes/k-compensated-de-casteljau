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

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

#ifndef DE_CASTELJAU_H
#define DE_CASTELJAU_H

#include <stdio.h>

#if defined(__cplusplus)
extern "C" {
#endif

double basic(double s, const double* coeffs, double* pk, size_t degree);
void compensated(double s, const double* coeffs, size_t degree, size_t K,
    double* errors, double* bk, double* result);

#if defined(__cplusplus)
}
#endif

#endif /* DE_CASTELJAU_H */

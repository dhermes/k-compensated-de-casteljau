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

#include <math.h>
#include <stdio.h>

void two_prod(double val1, double val2, double* product, double* error)
{
    *product = val1 * val2;
    *error = fma(val1, val2, -*product);
}

void two_sum(double val1, double val2, double* sum, double* error)
{
    *sum = val1 + val2;
    double almost_val2 = *sum - val1;
    *error = (val1 - (*sum - almost_val2)) + (val2 - almost_val2);
}

void vec_sum(double* vec, size_t n)
{
    // NOTE: This assumes that ``length(vec) == n``.
    for (size_t i = 1; i < n; ++i) {
        double sum, error;
        two_sum(vec[i], vec[i - 1], &sum, &error);
        vec[i] = sum;
        vec[i - 1] = error;
    }
}

double sum_k(double* vec, double* workspace, size_t n, size_t K)
{
    // NOTE: This assumes that ``length(vec) == n``.
    // NOTE: This requires ``length(workspace) >= length(vec)``.
    for (size_t i = 0; i < n; ++i) {
        workspace[i] = vec[i];
    }

    for (size_t i = 0; i < K - 1; ++i) {
        vec_sum(workspace, n);
    }

    double result = workspace[0];
    for (size_t i = 1; i < n; ++i) {
        result += workspace[i];
    }

    return result;
}

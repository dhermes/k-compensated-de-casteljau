# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

r"""Evaluates a polynomial in Bernstein form via the `VS method`_.

.. _VS method: https://doi.org/10.1016/0167-8396(86)90018-X
.. _survey: http://dx.doi.org/10.1016/j.amc.2015.08.086

This method has been described in the 1D case in a `survey`_
paper that followed the original formulation.
"""

import fractions
import math

import eft


def binomial(n, k):
    numerator = math.factorial(n)
    denominator = math.factorial(k) * math.factorial(n - k)
    result = fractions.Fraction(numerator, denominator)
    if float(result) != result:
        raise ValueError(n, k)
    return float(result)


def basic(s, coeffs):
    n = len(coeffs) - 1
    r = 1.0 - s

    result = coeffs[0]

    s_pow = 1.0
    for j in range(1, n + 1):
        s_pow = s * s_pow
        binom_val = binomial(n, j)
        result = r * result + binom_val * s_pow * coeffs[j]

    return result


def compensated(s, coeffs):
    n = len(coeffs) - 1
    r, rho = eft.add_eft(1.0, -s)

    pk = coeffs[0]
    dpk = 0.0

    s_pow = 1.0
    ds = 0.0
    for j in range(1, n + 1):
        # Update ``s^k``, using
        #   s^k = s_pow + ds ==> s^{k + 1} = s(s_pow) + s(ds)
        s_pow, pi_s = eft.multiply_eft(s, s_pow)
        ds = s * ds + pi_s
        # Now, update ``pk`` and ``dpk``.
        P1, pi1 = eft.multiply_eft(r, pk)
        local_err = pi1 + rho * pk
        P2, pi2 = eft.multiply_eft(coeffs[j], binomial(n, j))
        local_err += P2 * ds + pi2 * s_pow
        P3, pi3 = eft.multiply_eft(P2, s_pow)
        local_err += pi3
        pk, sigma4 = eft.add_eft(P1, P3)
        local_err += sigma4
        dpk = r * dpk + local_err

    return pk + dpk

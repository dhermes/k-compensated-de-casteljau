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

"""Performs de Casteljau's method.

de Casteljau's method evaluates a function in Bernstein-Bezier form

.. math::

    p(s) = p_0 (1 - s)^n + \cdots + p_j \binom{n}{j} (1 - s)^{n - j} s^j +
        \cdots + p_n s^n

by progressively reducing the control points

.. math::

    \begin{align*}
    p_j^{(0)} &= p_j \\
    p_j^{(k + 1)} &= (1 - s) p_j^{(k)} + p_{j + 1}^{(k)} \\
    p(s) &= p_0^{(n)}.
    \end{align*}

This module provides both the standard version and a compensated version.

.. note::

   This assumes throughout that ``coeffs`` is ordered from
   :math:`p_n` to :math:`p_0`.
"""


import eft


def basic(s, coeffs):
    r = 1.0 - s

    degree = len(coeffs) - 1
    pk = list(coeffs)
    for k in range(degree):
        new_pk = []
        for j in range(degree - k):
            new_pk.append(r * pk[j] + s * pk[j + 1])
        # Update the "current" values.
        pk = new_pk

    return pk[0]


def _compensated(s, coeffs):
    r, rho = eft.add_eft(1.0, -s)

    degree = len(coeffs) - 1
    pk = list(coeffs)
    ek = [0.0] * (degree + 1)
    for k in range(degree):
        new_pk = []
        new_ek = []

        for j in range(degree - k):
            prod1, pi = eft.multiply_eft(r, pk[j])
            prod2, sigma = eft.multiply_eft(s, pk[j + 1])
            sum_, beta = eft.add_eft(prod1, prod2)
            new_pk.append(sum_)
            # Now update the error.
            w = pi + sigma + beta + pk[j] * rho
            new_ek.append(w + s * ek[j + 1] + r * ek[j])

        # Update the "current" values.
        pk = new_pk
        ek = new_ek

    return pk[0], ek[0]


def compensated(x, coeffs):
    p, e = _compensated(x, coeffs)
    return p + e

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

"""Performs Horner's method.

Horner's method computes

.. math::

    p(x) = a_n x^n + \cdots a_1 x + a_0

via

.. math::

    \begin{align*}
    p_n &= a_n \\
    p_k &= p_{k + 1} x + a_k \\
    p(x) &= p_0
    \end{align*}

This module provides both the standard version and a compensated version.

.. note::

   This assumes throughout that ``coeffs`` is ordered from
   :math:`a_n` to :math:`a_0`.
"""

import eft


def basic(x, coeffs):
    p = coeffs[0]
    for coeff in coeffs[1:]:
        p = p * x + coeff

    return p


def _compensated(x, coeffs):
    p = coeffs[0]
    e = 0.0
    for coeff in coeffs[1:]:
        prod, e1 = eft.multiply_eft(p, x)
        p, e2 = eft.add_eft(prod, coeff)
        e = x * e + (e1 + e2)

    return p, e


def compensated(x, coeffs):
    p, e = _compensated(x, coeffs)
    return p + e

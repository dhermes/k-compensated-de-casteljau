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


def compensated(s, coeffs):
    p, e = _compensated(s, coeffs)
    return p + e


def _compensated3(s, coeffs):
    r, rho = eft.add_eft(1.0, -s)

    degree = len(coeffs) - 1
    pk = list(coeffs)
    eAk = [0.0] * (degree + 1)
    eBk = [0.0] * (degree + 1)
    for k in range(degree):
        new_pk = []
        new_eAk = []
        new_eBk = []

        for j in range(degree - k):
            # Update the "level 0" stuff.
            P1, pi1 = eft.multiply_eft(r, pk[j])
            P2, pi2 = eft.multiply_eft(s, pk[j + 1])
            S3, sigma3 = eft.add_eft(P1, P2)
            new_pk.append(S3)
            # Update the "level 1" stuff.
            S4, sigma4 = eft.add_eft(pi1, pi2)
            S5, sigma5 = eft.add_eft(S4, sigma3)
            P6, pi6 = eft.multiply_eft(rho, pk[j])
            wA, sigma7 = eft.add_eft(S5, P6)
            P8, pi8 = eft.multiply_eft(s, eAk[j + 1])
            S9, sigma9 = eft.add_eft(P8, wA)
            P10, pi10 = eft.multiply_eft(r, eAk[j])
            S11, sigma11 = eft.add_eft(S9, P10)
            new_eAk.append(S11)
            # Update the "level 2" stuff.
            wB = (
                sigma4
                + sigma5
                + pi6
                + sigma7
                + pi8
                + sigma9
                + pi10
                + sigma11
                + rho
                * eAk[j]
            )
            new_eBk.append(wB + s * eBk[j + 1] + r * eBk[j])

        # Update the "current" values.
        pk = new_pk
        eAk = new_eAk
        eBk = new_eBk

    return pk[0], eAk[0], eBk[0]


def compensated3(s, coeffs):
    p, eA, eB = _compensated3(s, coeffs)
    return (p + eA) + eB


def _compensated4(s, coeffs):
    r, rho = eft.add_eft(1.0, -s)

    degree = len(coeffs) - 1
    pk = list(coeffs)
    eAk = [0.0] * (degree + 1)
    eBk = [0.0] * (degree + 1)
    eCk = [0.0] * (degree + 1)
    for k in range(degree):
        new_pk = []
        new_eAk = []
        new_eBk = []
        new_eCk = []

        for j in range(degree - k):
            # Update the "level 0" stuff.
            P1, pi1 = eft.multiply_eft(r, pk[j])
            P2, pi2 = eft.multiply_eft(s, pk[j + 1])
            S3, sigma3 = eft.add_eft(P1, P2)
            new_pk.append(S3)
            # Update the "level 1" stuff.
            S4, sigma4 = eft.add_eft(pi1, pi2)
            S5, sigma5 = eft.add_eft(S4, sigma3)
            P6, pi6 = eft.multiply_eft(rho, pk[j])
            wA, sigma7 = eft.add_eft(S5, P6)
            P8, pi8 = eft.multiply_eft(s, eAk[j + 1])
            S9, sigma9 = eft.add_eft(P8, wA)
            P10, pi10 = eft.multiply_eft(r, eAk[j])
            S11, sigma11 = eft.add_eft(S9, P10)
            new_eAk.append(S11)
            # Update the "level 2" stuff.
            S12, sigma12 = eft.add_eft(sigma4, sigma5)
            S14, sigma14 = eft.add_eft(S12, pi6)
            S15, sigma15 = eft.add_eft(S14, sigma7)
            S16, sigma16 = eft.add_eft(S15, pi8)
            S17, sigma17 = eft.add_eft(S16, sigma9)
            S18, sigma18 = eft.add_eft(S17, pi10)
            S19, sigma19 = eft.add_eft(S18, sigma11)
            P13, pi13 = eft.multiply_eft(rho, eAk[j])
            wB, sigma20 = eft.add_eft(S19, P13)
            P21, pi21 = eft.multiply_eft(s, eBk[j + 1])
            S22, sigma22 = eft.add_eft(P21, wB)
            P23, pi23 = eft.multiply_eft(r, eBk[j])
            S24, sigma24 = eft.add_eft(P23, S22)
            new_eBk.append(S24)
            # Update the "level 3" stuff.
            wC = (
                sigma12
                + pi13
                + sigma14
                + sigma15
                + sigma16
                + sigma17
                + sigma18
                + sigma19
                + sigma20
                + pi21
                + sigma22
                + pi23
                + sigma24
                + rho
                * eBk[j]
            )
            new_eCk.append(wC + s * eCk[j + 1] + r * eCk[j])

        # Update the "current" values.
        pk = new_pk
        eAk = new_eAk
        eBk = new_eBk
        eCk = new_eCk

    return pk[0], eAk[0], eBk[0], eCk[0]


def compensated4(s, coeffs):
    p, eA, eB, eC = _compensated4(s, coeffs)
    return ((p + eA) + eB) + eC

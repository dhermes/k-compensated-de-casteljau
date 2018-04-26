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
    eAk = [0.0] * (degree + 1)
    for k in range(degree):
        new_pk = []
        new_eAk = []

        for j in range(degree - k):
            P1, pi1 = eft.multiply_eft(r, pk[j])
            P2, pi2 = eft.multiply_eft(s, pk[j + 1])
            S3, sigma3 = eft.add_eft(P1, P2)
            new_pk.append(S3)
            # Now update the error.
            wA = pi1 + pi2 + sigma3 + rho * pk[j]
            new_eAk.append(wA + s * eAk[j + 1] + r * eAk[j])

        # Update the "current" values.
        pk = new_pk
        eAk = new_eAk

    return pk[0], eAk[0]


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
            S4, sigma4 = eft.add_eft(pi1, pi2)  # wA = pi1 + pi2 + ...
            S5, sigma5 = eft.add_eft(S4, sigma3)  # ... + sigma3 + ...
            P6, pi6 = eft.multiply_eft(rho, pk[j])
            wA, sigma7 = eft.add_eft(S5, P6)  # ... + rho * pk[j]
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
            # wB = sigma4 + sigma5 + ...
            S12, sigma12 = eft.add_eft(sigma4, sigma5)
            S14, sigma14 = eft.add_eft(S12, pi6)  # ... + pi6 + ...
            S15, sigma15 = eft.add_eft(S14, sigma7)  # ... + sigma7 + ...
            S16, sigma16 = eft.add_eft(S15, pi8)  # ... + pi8 + ...
            S17, sigma17 = eft.add_eft(S16, sigma9)  # ... + sigma9 + ...
            S18, sigma18 = eft.add_eft(S17, pi10)  # ... + pi10 + ...
            S19, sigma19 = eft.add_eft(S18, sigma11)  # ... + sigma11 + ...
            P13, pi13 = eft.multiply_eft(rho, eAk[j])
            wB, sigma20 = eft.add_eft(S19, P13)  # ... + rho eAk[j] + ...
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


def _compensated5(s, coeffs):
    r, rho = eft.add_eft(1.0, -s)

    degree = len(coeffs) - 1
    pk = list(coeffs)
    eAk = [0.0] * (degree + 1)
    eBk = [0.0] * (degree + 1)
    eCk = [0.0] * (degree + 1)
    eDk = [0.0] * (degree + 1)
    for k in range(degree):
        new_pk = []
        new_eAk = []
        new_eBk = []
        new_eCk = []
        new_eDk = []

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
            # wC = sigma12 + pi13 + ...
            S25, sigma25 = eft.add_eft(sigma12, pi13)
            S26, sigma26 = eft.add_eft(S25, sigma14)  # ... + sigma14 + ...
            S27, sigma27 = eft.add_eft(S26, sigma15)  # ... + sigma15 + ...
            S28, sigma28 = eft.add_eft(S27, sigma16)  # ... + sigma16 + ...
            S29, sigma29 = eft.add_eft(S28, sigma17)  # ... + sigma17 + ...
            S30, sigma30 = eft.add_eft(S29, sigma18)  # ... + sigma18 + ...
            S31, sigma31 = eft.add_eft(S30, sigma19)  # ... + sigma19 + ...
            S32, sigma32 = eft.add_eft(S31, sigma20)  # ... + sigma20 + ...
            S33, sigma33 = eft.add_eft(S32, pi21)  # ... + pi21 + ...
            S34, sigma34 = eft.add_eft(S33, sigma22)  # ... + sigma22 + ...
            S35, sigma35 = eft.add_eft(S34, pi23)  # ... + pi23 + ...
            S36, sigma36 = eft.add_eft(S35, sigma24)  # ... + sigma24 + ...
            P37, pi37 = eft.multiply_eft(rho, eBk[j])
            wC, sigma38 = eft.add_eft(S36, P37)  # ... + rho * eBk[j]
            P39, pi39 = eft.multiply_eft(s, eCk[j + 1])
            S40, sigma40 = eft.add_eft(P39, wC)
            P41, pi41 = eft.multiply_eft(r, eCk[j])
            S42, sigma42 = eft.add_eft(P41, S40)
            new_eCk.append(S42)
            # Update the "level 4" stuff.
            wD = (
                sigma25
                + sigma26
                + sigma27
                + sigma28
                + sigma29
                + sigma30
                + sigma31
                + sigma32
                + sigma33
                + sigma34
                + sigma35
                + sigma36
                + pi37
                + sigma38
                + pi39
                + sigma40
                + pi41
                + sigma42
                + rho
                * eCk[j]
            )
            new_eDk.append(wD + s * eDk[j + 1] + r * eDk[j])

        # Update the "current" values.
        pk = new_pk
        eAk = new_eAk
        eBk = new_eBk
        eCk = new_eCk
        eDk = new_eDk

    return pk[0], eAk[0], eBk[0], eCk[0], eDk[0]


def compensated5(s, coeffs):
    p, eA, eB, eC, eD = _compensated5(s, coeffs)
    return (((p + eA) + eB) + eC) + eD

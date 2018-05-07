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

"""Verify table in manuscript.

This performs compensated de Casteljau (``K = 2``) in both IEEE-754
double precision and with exact arithmetic and verifies that the terms
in the table are correct.
"""

import fractions

import eft


F = fractions.Fraction
K = 1001
ROOT = 0.5
U = F(1, 2 ** 53)
MACH = 0.5 ** 53
S = ROOT + K * MACH
R_HAT, RHO = eft.add_eft(1.0, -S)


def exact_de_castlejau(b4_vals):
    b = {}
    b[4] = tuple(map(F, b4_vals))
    s = F(S)

    for k in (3, 2, 1, 0):
        prev_b = b[k + 1]
        b[k] = tuple(
            (1 - s) * b1 + s * b2 for b1, b2 in zip(prev_b, prev_b[1:])
        )

    return b


def stage1(b4_vals, exact_b):
    b04, b14, b24, b34, b44 = b4_vals

    P1, pi1 = eft.multiply_eft(R_HAT, b04)
    P2, pi2 = eft.multiply_eft(S, b14)
    b03, sigma3 = eft.add_eft(P1, P2)
    e03 = pi1 + pi2 + sigma3 + RHO * b04

    P1, pi1 = eft.multiply_eft(R_HAT, b14)
    P2, pi2 = eft.multiply_eft(S, b24)
    b13, sigma3 = eft.add_eft(P1, P2)
    e13 = pi1 + pi2 + sigma3 + RHO * b14

    P1, pi1 = eft.multiply_eft(R_HAT, b24)
    P2, pi2 = eft.multiply_eft(S, b34)
    b23, sigma3 = eft.add_eft(P1, P2)
    e23 = pi1 + pi2 + sigma3 + RHO * b24

    P1, pi1 = eft.multiply_eft(R_HAT, b34)
    P2, pi2 = eft.multiply_eft(S, b44)
    b33, sigma3 = eft.add_eft(P1, P2)
    e33 = pi1 + pi2 + sigma3 + RHO * b34

    # Verify columns 1 (bhat) and 2 (dbhat).
    assert F(b03) == F(1, 8) - 7 * (K * U) / 4 - U / 4
    assert F(b13) == -F(1, 8) + 5 * (K * U) / 4 + U / 4
    assert F(b23) == F(1, 8) - 3 * (K * U) / 4
    assert F(b33) == -F(1, 8) + (K * U) / 4
    assert F(e03) == U / 4
    assert F(e13) == -U / 4
    assert F(e23) == 0
    assert F(e33) == 0
    # Verify column 2 (d2b).
    exact_b03, exact_b13, exact_b23, exact_b33 = exact_b[3]
    assert exact_b03 - F(b03) - F(e03) == 0
    assert exact_b13 - F(b13) - F(e13) == 0
    assert exact_b23 - F(b23) - F(e23) == 0
    assert exact_b33 - F(b33) - F(e33) == 0

    b_vals = (b03, b13, b23, b33)
    e_vals = (e03, e13, e23, e33)
    return b_vals, e_vals


def stage2(b3_vals, e3_vals, exact_b):
    b03, b13, b23, b33 = b3_vals
    e03, e13, e23, e33 = e3_vals

    P1, pi1 = eft.multiply_eft(R_HAT, b03)
    P2, pi2 = eft.multiply_eft(S, b13)
    b02, sigma3 = eft.add_eft(P1, P2)
    ell02 = pi1 + pi2 + sigma3 + RHO * b03
    e02 = ell02 + S * e13 + R_HAT * e03

    P1, pi1 = eft.multiply_eft(R_HAT, b13)
    P2, pi2 = eft.multiply_eft(S, b23)
    b12, sigma3 = eft.add_eft(P1, P2)
    ell12 = pi1 + pi2 + sigma3 + RHO * b13
    e12 = ell12 + S * e23 + R_HAT * e13

    P1, pi1 = eft.multiply_eft(R_HAT, b23)
    P2, pi2 = eft.multiply_eft(S, b33)
    b22, sigma3 = eft.add_eft(P1, P2)
    ell22 = pi1 + pi2 + sigma3 + RHO * b23
    e22 = ell22 + S * e33 + R_HAT * e23

    # Verify columns 1 (bhat) and 2 (dbhat).
    assert F(b02) == -(K * U) / 2
    assert F(b12) == (K * U) / 2 + U / 8
    assert F(b22) == -(K * U) / 2
    assert F(e02) == 3 * (K * U) ** 2
    assert F(e12) == -U / 8 - 2 * (K * U) ** 2
    assert F(e22) == (K * U) ** 2
    # Verify column 2 (d2b).
    exact_b02, exact_b12, exact_b22 = exact_b[2]
    assert exact_b02 - F(b02) - F(e02) == 0
    assert exact_b12 - F(b12) - F(e12) == 0
    assert exact_b22 - F(b22) - F(e22) == 0

    b_vals = (b02, b12, b22)
    e_vals = (e02, e12, e22)
    return b_vals, e_vals


def stage3(b2_vals, e2_vals, exact_b):
    b02, b12, b22 = b2_vals
    e02, e12, e22 = e2_vals

    P1, pi1 = eft.multiply_eft(R_HAT, b02)
    P2, pi2 = eft.multiply_eft(S, b12)
    b01, sigma3 = eft.add_eft(P1, P2)
    ell01 = pi1 + pi2 + sigma3 + RHO * b02
    e01 = ell01 + S * e12 + R_HAT * e02

    P1, pi1 = eft.multiply_eft(R_HAT, b12)
    P2, pi2 = eft.multiply_eft(S, b22)
    b11, sigma3 = eft.add_eft(P1, P2)
    ell11 = pi1 + pi2 + sigma3 + RHO * b12
    e11 = ell11 + S * e22 + R_HAT * e12

    # Verify columns 1 (bhat) and 2 (dbhat).
    assert F(b01) == U / 16 + (K * U) ** 2 + 239 * U ** 2
    assert F(b11) == U / 16 - (K * U) ** 2 - 239 * U ** 2
    assert F(e01) == -U / 16 + (K * U) ** 2 / 2 - 239 * U ** 2
    assert F(e11) == -U / 16 - (K * U) ** 2 / 2 + 239 * U ** 2
    # Verify column 2 (d2b).
    exact_b01, exact_b11 = exact_b[1]
    assert exact_b01 - F(b01) - F(e01) == -5 * (K * U) ** 3
    assert exact_b11 - F(b11) - F(e11) == 3 * (K * U) ** 3

    b_vals = (b01, b11)
    e_vals = (e01, e11)
    return b_vals, e_vals


def stage4(b1_vals, e1_vals, exact_b):
    b01, b11 = b1_vals
    e01, e11 = e1_vals

    P1, pi1 = eft.multiply_eft(R_HAT, b01)
    P2, pi2 = eft.multiply_eft(S, b11)
    b00, sigma3 = eft.add_eft(P1, P2)
    ell00 = pi1 + pi2 + sigma3 + RHO * b01
    e00 = ell00 + S * e11 + R_HAT * e01

    # Verify columns 1 (bhat) and 2 (dbhat).
    assert F(b00) == U / 16
    assert F(e00) == -U / 16
    # Verify column 2 (d2b).
    exact_b00, = exact_b[0]
    assert exact_b00 - F(b00) - F(e00) == -4 * (K * U) ** 3 + 8 * (K * U) ** 4


def main():
    # First, verify S, R_HAT and RHO.
    assert F(S) == F(1, 2) + K * U
    assert F(RHO) == 0
    assert F(R_HAT) == F(1, 2) - K * U

    # p(s) = (2s - 1)^3 (s - 1)
    b4_vals = (1.0, -0.75, 0.5, -0.25, 0.0)
    exact_b = exact_de_castlejau(b4_vals)

    b3_vals, e3_vals = stage1(b4_vals, exact_b)
    b2_vals, e2_vals = stage2(b3_vals, e3_vals, exact_b)
    b1_vals, e1_vals = stage3(b2_vals, e2_vals, exact_b)
    stage4(b1_vals, e1_vals, exact_b)


if __name__ == "__main__":
    main()

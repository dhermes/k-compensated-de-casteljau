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

"""Run many experiments to confirm the flop counts for algorithms."""

from __future__ import print_function

import de_casteljau
import eft
import horner
import operation_count


SEPARATOR = "-" * 80


def count_add_eft():
    parent = operation_count.Computation()
    val1 = operation_count.Float(1.5, parent)
    val2 = operation_count.Float(0.5 + 0.5 ** 52, parent)
    sum_, error = eft.add_eft(val1, val2)
    assert sum_.value == 2.0
    assert error.value == 0.5 ** 52
    print("     add_eft(): {}".format(parent.display))


def count__split():
    parent = operation_count.Computation()
    val = operation_count.Float(1.0 + 0.5 ** 27, parent)
    high, low = eft._split(val)
    assert high.value == 1.0
    assert low.value == 0.5 ** 27
    print("      _split(): {}".format(parent.display))


def count_multiply_eft():
    print("multiply_eft():")
    for use_fma in (True, False):
        parent = operation_count.Computation()
        val1 = operation_count.Float(1.0 + 0.5 ** 40, parent)
        val2 = operation_count.Float(1.0 - 0.5 ** 40, parent)
        product, error = eft.multiply_eft(val1, val2, use_fma=use_fma)
        assert product.value == 1.0
        assert error.value == -0.5 ** 80
        if use_fma:
            description = "with FMA:    "
        else:
            description = "w / out FMA: "
        print("  {} {}".format(description, parent.display))


def count__vec_sum():
    print("_vec_sum() (6(|p| - 1)):")
    for size_p in range(1, 9 + 1):
        parent = operation_count.Computation()
        p = [operation_count.Float(1.0, parent)] * size_p
        eft._vec_sum(p)

        assert p[size_p - 1].value == float(size_p)
        print("  |p| = {}:     {}".format(size_p, parent.display))


def count_sum_k():
    print("sum_k() ((6K - 5)(|p| - 1)):")
    for k in (2, 3, 4, 5):
        print("  K = {}".format(k))
        for size_p in range(1, 9 + 1):
            parent = operation_count.Computation()
            p = [operation_count.Float(1.0, parent)] * size_p
            total = eft.sum_k(p, k)

            assert total.value == float(size_p)
            print("    |p| = {}:    {}".format(size_p, parent.display))


def count_horner_basic():
    print("horner.basic() (2n):")
    for degree in range(1, 9 + 1):
        parent = operation_count.Computation()
        x = operation_count.Float(2.0, parent)
        coeffs = (operation_count.Float(1.0, parent),) * (degree + 1)
        p = horner.basic(x, coeffs)
        assert p.value == 2.0 ** (degree + 1) - 1
        print("  degree {}:     {}".format(degree, parent.display))


def count_horner_compensated():
    print("horner.compensated() (11n + 1):")
    for degree in range(1, 9 + 1):
        parent = operation_count.Computation()
        x = operation_count.Float(2.0, parent)
        coeffs = (operation_count.Float(1.0, parent),) * (degree + 1)
        p = horner.compensated(x, coeffs)
        assert p.value == 2.0 ** (degree + 1) - 1
        print("  degree {}:     {}".format(degree, parent.display))


def count_de_casteljau_basic():
    print("de_casteljau.basic() ((3n^2 + 3n + 2) / 2 = 3 T_n + 1):")
    for degree in range(1, 9 + 1):
        parent = operation_count.Computation()
        x = operation_count.Float(0.25, parent)
        coeffs = tuple(
            operation_count.Float((-1.0) ** k, parent)
            for k in range(degree + 1)
        )
        p = de_casteljau.basic(x, coeffs)
        assert p.value == 0.5 ** degree
        print("  degree {}:     {}".format(degree, parent.display))


def count_de_casteljau_compensated():
    print("de_casteljau.compensated() (9n^2 + 9n + 7 = 18 T_n + 7):")
    for degree in range(1, 9 + 1):
        parent = operation_count.Computation()
        x = operation_count.Float(0.25, parent)
        coeffs = tuple(
            operation_count.Float((-1.0) ** k, parent)
            for k in range(degree + 1)
        )
        p = de_casteljau.compensated(x, coeffs)
        assert p.value == 0.5 ** degree
        print("  degree {}:     {}".format(degree, parent.display))


def count_de_casteljau_compensated3():
    print("de_casteljau.compensated3() ((59n^2 + 59n + 16) / 2 = 59 T_n + 8):")
    for degree in range(1, 9 + 1):
        parent = operation_count.Computation()
        x = operation_count.Float(0.25, parent)
        coeffs = tuple(
            operation_count.Float((-1.0) ** k, parent)
            for k in range(degree + 1)
        )
        p = de_casteljau.compensated3(x, coeffs)
        assert p.value == 0.5 ** degree
        print("  degree {}:     {}".format(degree, parent.display))


def main():
    count_add_eft()
    print(SEPARATOR)
    count__split()
    print(SEPARATOR)
    count_multiply_eft()
    print(SEPARATOR)
    count__vec_sum()
    print(SEPARATOR)
    count_sum_k()
    print(SEPARATOR)
    count_horner_basic()
    print(SEPARATOR)
    count_horner_compensated()
    print(SEPARATOR)
    count_de_casteljau_basic()
    print(SEPARATOR)
    count_de_casteljau_compensated()
    print(SEPARATOR)
    count_de_casteljau_compensated3()


if __name__ == "__main__":
    main()

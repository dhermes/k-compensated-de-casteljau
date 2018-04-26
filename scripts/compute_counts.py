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

import eft
import horner
import operation_count


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


def count_horner_basic():
    print("horner.basic() (2n):")
    for degree in range(1, 9 + 1):
        parent = operation_count.Computation()
        x = operation_count.Float(0.0, parent)
        coeffs = (operation_count.Float(0.0, parent),) * (degree + 1)
        horner.basic(x, coeffs)
        print("  degree {}:     {}".format(degree, parent.display))


def count_horner_compensated():
    print("horner.compensated() (11n + 1):")
    for degree in range(1, 9 + 1):
        parent = operation_count.Computation()
        x = operation_count.Float(0.0, parent)
        coeffs = (operation_count.Float(0.0, parent),) * (degree + 1)
        horner.compensated(x, coeffs)
        print("  degree {}:     {}".format(degree, parent.display))


def main():
    count_add_eft()
    count__split()
    count_multiply_eft()
    count_horner_basic()
    count_horner_compensated()


if __name__ == "__main__":
    main()

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

"""Helpers for counting flops."""


import fractions


_DISPLAY_TEMPLATE = (
    "{:4d} flops ({:4d} add, {:4d} sub, {:3d} multiply, {:3d} FMA)"
)


class Computation(object):
    """Stateful manager of the number of flops."""

    def __init__(self):
        self.add_count = 0
        self.sub_count = 0
        self.mul_count = 0
        self.fma_count = 0

    @property
    def count(self):
        return (
            self.add_count + self.sub_count + self.mul_count + self.fma_count
        )

    @property
    def display(self):
        return _DISPLAY_TEMPLATE.format(
            self.count,
            self.add_count,
            self.sub_count,
            self.mul_count,
            self.fma_count,
        )


class Float(object):
    """A ``float``-like type that will increment a flop count.

    Args:
        value (float): The current value.
        computation (.Computation): The current computation
            in progress.
    """

    def __init__(self, value, computation):
        self.value = value
        self.computation = computation

    def _get_value(self, other):
        if isinstance(other, Float):
            if other.computation is not self.computation:
                raise ValueError(
                    "Two `Float`s being combined should have the "
                    "same parent computation."
                )
            return other.value
        elif isinstance(other, float):
            return other
        elif isinstance(other, int):
            other_float = float(other)
            if other == other_float:
                return other_float
            else:
                return None
        else:
            return None

    def __add__(self, other):
        value = self._get_value(other)
        if value is None:
            return NotImplemented

        self.computation.add_count += 1
        return Float(self.value + value, self.computation)

    def __radd__(self, other):
        value = self._get_value(other)
        if value is None:
            return NotImplemented

        self.computation.add_count += 1
        return Float(value + self.value, self.computation)

    def __sub__(self, other):
        value = self._get_value(other)
        if value is None:
            return NotImplemented

        self.computation.sub_count += 1
        return Float(self.value - value, self.computation)

    def __rsub__(self, other):
        value = self._get_value(other)
        if value is None:
            return NotImplemented

        self.computation.sub_count += 1
        return Float(value - self.value, self.computation)

    def __neg__(self):
        return Float(-self.value, self.computation)

    def __mul__(self, other):
        value = self._get_value(other)
        if value is None:
            return NotImplemented

        self.computation.mul_count += 1
        return Float(self.value * value, self.computation)

    def __rmul__(self, other):
        value = self._get_value(other)
        if value is None:
            return NotImplemented

        self.computation.mul_count += 1
        return Float(value * self.value, self.computation)

    def __truediv__(self, other):
        value = self._get_value(other)
        if value is None:
            return NotImplemented

        self.computation.mul_count += 1
        return Float(self.value / value, self.computation)

    def fma(self, val1, val2, val3):
        float1 = self._get_value(val1)
        float2 = self._get_value(val2)
        float3 = self._get_value(val3)
        if None in (float1, float2, float3):
            raise TypeError("Only `Float` or `float` allowed in fma")

        self.computation.fma_count += 1
        frac1 = fractions.Fraction(float1)
        frac2 = fractions.Fraction(float2)
        frac3 = fractions.Fraction(float3)
        result = float(frac1 * frac2 + frac3)
        return Float(result, self.computation)

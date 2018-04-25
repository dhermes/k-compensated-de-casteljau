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

"""Collection of error-free transforms."""


def add_eft(val1, val2):
    # See: https://doi.org/10.1137/030601818
    sum_ = val1 + val2
    delta1 = sum_ - val1
    error = (val1 - (sum_ - delta1)) + (val2 - delta1)
    return sum_, error


def _split(val):
    # Helper for ``multiply_eft``.
    scaled = val * 67108865.0  # 67108865 == 2^{26} + 1.
    high_bits = scaled - (scaled - val)
    low_bits = val - high_bits
    return high_bits, low_bits


def multiply_eft(val1, val2):
    # See: https://doi.org/10.1109/TC.2008.215
    product = val1 * val2
    high1, low1 = _split(val1)
    high2, low2 = _split(val2)
    error = (
        low1
        * low2
        - (((product - high1 * high2) - low1 * high2) - high1 * low2)
    )
    return product, error

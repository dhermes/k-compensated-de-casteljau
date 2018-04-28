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

r"""Perform plots of relative error against the condition number.

This uses :math:`p(s) = (s - 1) \left(s - \frac{3}{4}\right)^7` which
has :math:`\widetilde{p}(s) = (s - 1) \left(\frac{s}{2} -
\frac{3}{4}\right)^7`.
"""

import fractions
import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn

import de_casteljau


F = fractions.Fraction
U = F(1, 2 ** 53)
# p(s) = (s - 1) (s - 3/4)^7 = SUM_j b_j B_{j, 8}(s)
BEZIER_COEFFS = (
    2187.0 / 16384.0,
    -5103.0 / 131072.0,
    729.0 / 65536.0,
    -405.0 / 131072.0,
    27.0 / 32768.0,
    -27.0 / 131072.0,
    3.0 / 65536.0,
    -1.0 / 131072.0,
    0.0,
)


def get_path(filename):
    """Get a file path in the ``images/`` directory.

    This assumes the script is currently in the ``scripts/``
    directory.
    """
    curr_dir = os.path.abspath(os.path.dirname(__file__))
    root_dir = os.path.dirname(curr_dir)
    images_dir = os.path.join(root_dir, "images")
    return os.path.join(images_dir, filename)


def main(filename=None):
    # n = 8
    gamma2n = 16 * U / (1 - 16 * U)
    gamma3n = 24 * U / (1 - 24 * U)

    cond_nums = []
    forward_errs1 = []
    forward_errs2 = []
    a_priori_bound2 = []
    for j in range(-5, -44 - 1, -1):
        s = 0.75 - 1.3 ** j
        exact_s = F(s)

        # Compute the condition number.
        exact_p = (exact_s - 1) * (4 * exact_s - 3) ** 7 / 16384
        # p_tilde(s) = SUM_j |b_j| B_{j, 8}(s) = (s - 1) (s/2 - 3/4)^7
        exact_p_tilde = (exact_s - 1) * (2 * exact_s - 3) ** 7 / 16384
        exact_cond = abs(exact_p_tilde / exact_p)
        cond_nums.append(float(exact_cond))

        # Update the lines of the a priori bounds.
        a_priori_bound2.append(float(U + 2 * gamma3n ** 2 * exact_cond))

        # Compute the forward error for uncompensated de Casteljau.
        b1, db, _ = de_casteljau._compensated3(s, BEZIER_COEFFS)
        exact_b1 = F(b1)
        exact_forward_err1 = abs((exact_b1 - exact_p) / exact_p)
        forward_errs1.append(float(exact_forward_err1))

        # Compute the forward error for compensated de Casteljau.
        b2 = b1 + db
        exact_b2 = F(b2)
        exact_forward_err2 = abs((exact_b2 - exact_p) / exact_p)
        forward_errs2.append(float(exact_forward_err2))

    # Set a tight ``x``-limit.
    min_exp = np.log(min(cond_nums))
    max_exp = np.log(max(cond_nums))
    delta_exp = max_exp - min_exp
    min_x = np.exp(min_exp - 0.01 * delta_exp)
    max_x = np.exp(max_exp + 0.01 * delta_exp)

    figure = plt.figure()
    ax = figure.gca()
    ax.loglog(
        cond_nums,
        forward_errs1,
        marker="v",
        linestyle="none",
        zorder=2,
        label=r"$\mathtt{DeCasteljau}$",
    )
    ax.loglog(
        cond_nums,
        forward_errs2,
        marker="d",
        linestyle="none",
        zorder=2,
        label=r"$\mathtt{CompDeCasteljau}$",
    )
    # Plot the lines of the a priori error bounds.
    ax.loglog(
        [min_x, max_x],
        [float(gamma2n * min_x), float(gamma2n * max_x)],
        color="black",
        zorder=1,
    )
    ax.loglog(
        [min_x] + cond_nums,
        [float(U + 2 * gamma3n ** 2 * min_x)] + a_priori_bound2,
        color="black",
        zorder=1,
    )
    # Figure out the bounds before adding the ``1/u`` and ``1/u^2`` lines.
    min_y, max_y = ax.get_ylim()
    delta_y = max_y - min_y
    ax.loglog(
        [1.0 / float(U), 1.0 / float(U)],
        [min_y - 0.05 * delta_y, max_y + 0.05 * delta_y],
        color="black",
        linestyle="dashed",
        zorder=1,
    )
    ax.loglog(
        [1.0 / float(U) ** 2, 1.0 / float(U) ** 2],
        [min_y - 0.05 * delta_y, max_y + 0.05 * delta_y],
        color="black",
        linestyle="dashed",
        zorder=1,
    )

    # Make sure the ``y``-limit stays set (the bounds lines exceed).
    ax.set_ylim(min_y, 1.0)
    ax.set_xlim(min_x, max_x)
    # Add the legend.
    ax.legend(loc="lower right", framealpha=1.0, frameon=True)
    # Set "nice" ticks.
    ax.set_xticks([10.0 ** n for n in range(5, 35 + 5, 5)])
    ax.set_yticks([10.0 ** n for n in range(-18, 0 + 2, 2)])
    # Set special ``xticks`` for ``1/u`` and ``1/u^2``.
    ax.set_xticks([1.0 / float(U), 1.0 / float(U) ** 2], minor=True)
    ax.set_xticklabels([r"$1/\mathbf{u}$", r"$1/\mathbf{u}^2$"], minor=True)
    ax.tick_params(
        axis="x",
        which="minor",
        direction="out",
        top=1,
        bottom=0,
        labelbottom=0,
        labeltop=1,
    )
    # Label the axes.
    ax.set_xlabel("Condition Number")
    ax.set_ylabel("Relative Forward Error")

    if filename is None:
        plt.show()
    else:
        path = get_path(filename)
        figure.savefig(path, bbox_inches="tight")
        print("Saved {}".format(filename))
        plt.close(figure)


if __name__ == "__main__":
    seaborn.set()
    main(filename="de_casteljau_rel_error.pdf")

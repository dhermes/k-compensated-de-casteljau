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

"""This is a configuration file for running ``nox`` on this project.

To determine the supported actions run ``nox --list-sessions`` from the
project root.
"""

import os

import nox
import py.path


NOX_DIR = os.path.abspath(os.path.dirname(__file__))
SINGLE_INTERP = "python3.6"
MATPLOTLIB_DEPS = (
    "cycler == 0.10.0",
    "kiwisolver == 1.0.1",
    "matplotlib == 2.2.2",
    "numpy == 1.14.3",
    "pyparsing == 2.2.0",
    "python-dateutil == 2.7.3",
    "pytz == 2018.4",
    "six == 1.11.0",
)
DEPS = {
    "matplotlib": MATPLOTLIB_DEPS,
    "seaborn": MATPLOTLIB_DEPS
    + ("pandas == 0.23.0", "scipy == 1.1.0", "seaborn == 0.8.1"),
}


def get_path(*names):
    return os.path.join(NOX_DIR, *names)


class Remove(object):

    def __init__(self, prefix, extensions):
        self.prefix = prefix
        self.extensions = extensions

    def __call__(self):
        for extension in self.extensions:
            path = "{}.{}".format(self.prefix, extension)
            os.remove(path)


def build_tex_file(session, base, new_id, extensions=(), with_bibtex=False):
    # NOTE: This assumes that ``session.chdir(get_path('doc'))``
    #       has been called.
    modify_id = get_path("scripts", "modify_pdf_id.py")

    if with_bibtex:
        session.run("pdflatex", base)
        session.run("bibtex", base)
        session.run("pdflatex", base)
        session.run("bibtex", base)
        session.run("pdflatex", base)
        session.run("pdflatex", base)
    else:
        session.run("pdflatex", base)
        session.run("pdflatex", base)
        session.run("pdflatex", base)

    path = get_path("doc", base)
    remove = Remove(path, extensions)
    session.run(remove)
    session.run("python", modify_id, "--base", path, "--id", new_id)


@nox.session
def build_tex(session):
    session.interpreter = SINGLE_INTERP

    if py.path.local.sysfind("pdflatex") is None:
        session.skip("`pdflatex` must be installed")

    if py.path.local.sysfind("bibtex") is None:
        session.skip("`bibtex` must be installed")

    # No need to create a virtualenv.
    session.virtualenv = False

    session.chdir(get_path("doc"))

    build_tex_file(
        session,
        "paper",
        "F092359D979FDC08931DA1922F3E123E",
        extensions=("aux", "bbl", "blg", "log", "out", "spl"),
        with_bibtex=True,
    )

    build_tex_file(
        session,
        "cover_letter",
        "BACA8D659970198BDF7D11B67FEA6299",
        extensions=("aux", "log", "out"),
    )

    build_tex_file(
        session,
        "tikz_local_err",
        "EF7ADBEFE6118EFEE506836A7AFF7C9E",
        extensions=("aux", "log", "out"),
    )

    build_tex_file(
        session,
        "tikz_filtration",
        "52F169AADDA4C4C85C6D5038361816C9",
        extensions=("aux", "log", "out"),
    )


@nox.session
def flop_counts(session):
    # No need to create a virtualenv.
    session.virtualenv = False

    env = {"PYTHONPATH": get_path("src")}
    compute_counts = get_path("scripts", "compute_counts.py")
    session.run("python", compute_counts, env=env)


@nox.session
def verify_table(session):
    # No need to create a virtualenv.
    session.virtualenv = False

    env = {"PYTHONPATH": get_path("src")}
    script = get_path("scripts", "verify_table.py")
    session.run("python", script, env=env)


@nox.session
def make_images(session):
    session.interpreter = SINGLE_INTERP
    # Install all dependencies.
    session.install(*DEPS["seaborn"])
    # Run the script(s).
    # Make sure
    # - Custom ``matplotlibrc`` is used
    # - Code in ``src/`` is importable
    # - PDFs have deterministic ``CreationDate``
    env = {
        "MATPLOTLIBRC": get_path("images"),
        "PYTHONPATH": get_path("src"),
        "SOURCE_DATE_EPOCH": "0",
    }
    names = (
        "error_against_cond.py",
        "smooth_drawing.py",
        "horner_inferior.py",
        "compensated_insufficient.py",
    )
    for name in names:
        script = get_path("scripts", name)
        session.run("python", script, env=env)

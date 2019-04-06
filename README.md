# ``K``-compensated de Casteljau

| Cite paper           | Cite code      |
| -------------------- | -------------- |
| [![Paper DOI][7]][5] | [![DOI][8]][9] |

This is the repository for my ``K``-compensated de Casteljau
paper. If you'd just like to read the [paper][1], feel
free.

This repository is laid out in a manner described in
[Good Enough Practices in Scientific Computing][2].

The content itself has been uploaded to the [arXiv][3] and was submitted to
the journal [AMC][4] in May 2018. The paper has been accepted and was
[published][5] [on][6] April 5, 2019.

## Abstract

In computer aided geometric design a polynomial is usually represented in
Bernstein form. This paper presents a family of compensated algorithms to
accurately evaluate a polynomial in Bernstein form with floating point
coefficients. The principle is to apply error-free transformations to
improve the traditional de Casteljau algorithm. At each stage of computation,
round-off error is passed on to first order errors, then to second order
errors, and so on. After the computation has been "filtered" `(K - 1)`
times via this process, the resulting output is as accurate as the de Casteljau
algorithm performed in `K` times the working precision. Forward error
analysis and numerical experiments illustrate the accuracy of this family
of algorithms.

## Citation

To cite this paper:

```
@article{Hermes2019,
  doi = {10.1016/j.amc.2019.03.047},
  url = {https://doi.org/10.1016/j.amc.2019.03.047},
  year = {2019},
  month = {Sep},
  publisher = {Elsevier {BV}},
  volume = {357},
  pages = {57--74},
  author = {Danny Hermes},
  title = {Compensated de {C}asteljau algorithm in {$K$} times the working precision},
  journal = {Applied Mathematics and Computation}
}
```

To cite the code in this repository:

```
@misc{KCompensatedGitHub,
  doi = {10.5281/zenodo.1405259},
  url = {https://zenodo.org/record/1405259},
  author = {Danny Hermes},
  title = {dhermes/k-compensated-de-casteljau: 2018.08.28},
  publisher = {Zenodo},
  year = {2018}
}
```

## Implementation

The ``K``-compensated de Casteljau algorithm is implemented in C, C++ and
Python in this repository. The implementations are contained in the
following source files:

### C

```
src/
├── de_casteljau.c
├── de_casteljau.h
├── eft.c
└── eft.h
```

### C++

```
src/
├── de_casteljau.cpp
├── de_casteljau.hpp
├── eft.cpp
└── eft.hpp
```

### Python

```
src/
├── de_casteljau.py
└── eft.py
```

## Installation

The code used to build the manuscript, generate images and verify
computations is written in Python. To run the code, Python 3.6
should be installed, along with ``nox``:

```
python -m pip install --upgrade 'nox >= 2018.10.17'
```

Once installed, the various build jobs can be listed. For example:

```
$ nox --list-sessions
Available sessions:
* build_tex
* flop_counts
* verify_table
* make_images
* update_requirements
* verify_cpp
* verify_c
```

To run ``nox -s build_tex`` (i.e. to build the PDF), ``pdflatex`` and
``bibtex`` are required.

## Plots

The plots can be (re)generated via ``nox -s make_images``.

## Operation Counts

A "special" numeric type is used to track flops and the actual operation
count for each algorithm is computed and verified via ``nox -s flop_counts``.

## Table of Computation

There is a table in the manuscript that details the **exact** floating point
values during a sample evaluation (at the point ``s = 1/2 + 1001u``).
These table values are verified via ``nox -s verify_table``.

[1]: doc/paper.pdf
[2]: https://arxiv.org/pdf/1609.00037.pdf
[3]: https://arxiv.org/abs/1808.10387
[4]: https://www.journals.elsevier.com/applied-mathematics-and-computation
[5]: https://doi.org/10.1016/j.amc.2019.03.047
[6]: doc/1-s2.0-S0096300319302541-main.pdf
[7]: https://img.shields.io/badge/DOI-10.1016%2Fj.amc.2019.03.047-blue.svg
[8]: https://zenodo.org/badge/131072021.svg
[9]: https://zenodo.org/badge/latestdoi/131072021

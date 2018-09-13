# ``K``-compensated de Casteljau

[![DOI](https://zenodo.org/badge/131072021.svg)](https://zenodo.org/badge/latestdoi/131072021)

This is the repository for my ``K``-compensated de Casteljau
paper. If you'd just like to read the [paper][1], feel
free.

This repository is laid out in a manner described in
[Good Enough Practices in Scientific Computing][2].

The content itself has been uploaded to the [arXiv][3] and submitted to
the journal [AMC][4] in May 2018.

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
python -m pip install --upgrade nox
```

Once installed, the various build jobs can be listed. For example:

```
$ nox --list-sessions
Available sessions:
* build_tex
* flop_counts
* verify_table
* make_images-3.6
* update_requirements-3.6
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

# ``K``-compensated de Casteljau

This is the repository for my ``K``-compensated de Casteljau
paper. If you'd just like to read the [paper][1], feel
free.

The paper demonstrates how to evaluate de Casteljau's Algorithm in
``K``-times the working precision.

This repository is laid out in a manner described in
[Good Enough Practices in Scientific Computing][2].

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

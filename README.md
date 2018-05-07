# ``K``-compensated de Casteljau

This is the repository for my ``K``-compensated de Casteljau
paper. If you'd just like to read the [paper][1], feel
free.

The paper demonstrates how to evaluate de Casteljau's Algorithm in
``K``-times the working precision.

## Installation

The code used to build the manuscript, generate images and verify
computations is written in Python. To run the code, Python 3.6
should be installed, along with ``nox-automation``:

```
python -m pip install --upgrade nox-automation
```

Once installed, the various build jobs can be listed. For example:

```
$ nox --list-sessions
Available sessions:
* build_tex
* flop_counts
* verify_table
* make_images
```

To run ``nox -s build_tex`` (i.e. to build the PDF), ``pdflatex`` is required.

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

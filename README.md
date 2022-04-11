# kimmy

Galactic chemical evolution in python

[![image](http://img.shields.io/pypi/v/kimmy.svg)](https://pypi.python.org/pypi/kimmy/)
[![Binder](http://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/jobovy/kimmy/main)

## Overview

``kimmy`` contains simple tools to study chemical evolution in galaxies.

## Author

Jo Bovy (University of Toronto): bovy - at - astro - dot - utoronto - dot - ca

## Installation

Install the latest release with ``pip``
```
pip install kimmy
```
or install the latest version by cloning/forking/downloading the repository and installing using
```
sudo python setup.py install
```
or locally using
```
python setup.py install --user
```

## Usage

For an example of usage, see the [example notebook](kimmy-example.ipynb). You can also launch a [Binder](https://mybinder.org/v2/gh/jobovy/kimmy/main) instance and directly play around with this notebook. _Or_ you can use [this](kimmy-example-pyodide.ipynb) [pyodide](https://github.com/pyodide/pyodide)-compatible version of the same notebook in [JupyterLite](https://github.com/jupyterlite/jupyterlite).

Currently, the only implemented feature is a simple one-zone chemical model with two elements ``O`` (for oxygen) and ``Fe`` (for iron). Initialize this model as
```
import kimmy
oz= kimmy.OneZone()
```
then for example compute the evolution of the default model and plot the [O/Fe] vs. [Fe/H] sequence
```
ts= numpy.linspace(0.001,10.,1001)*u.Gyr
plot(oz.Fe_H(ts),oz.O_Fe(ts))
```
To compute the distribution of [Fe/H], do for example,
```
FeHs= numpy.linspace(-1.525,1.225,56)
FeH_dist= [oz.Fe_H_DF(f) for f in FeHs]
```
and similar for the distribution of [O/H] and [O/Fe]. You can directly update the main parameters of the model and the model will be re-computed. For example, to set the outflow mass-loading parameter to one and plot the [O/Fe] vs. [Fe/H] sequence, do
```
ts= numpy.linspace(0.001,10.,1001)*u.Gyr
oz.eta= 1.
plot(oz.Fe_H(ts),oz.O_Fe(ts))
```
Keep in mind that once you change a parameter, it remains changed in the model. If you want to go back to the initial set of parameters that you used to initialize the instance, use ``oz.initial()``; if you want to go back to the default set of parameters, use ``oz.default()``. If you want to print the model you are using at any time, do
```
print(oz)
```
which prints a nicely formatted list of all of the model parameters.

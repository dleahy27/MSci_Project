#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from cubicspline import CubicSpline

# latex font configs
plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['text.usetex'] = True


# test function
def cubic(x):
    return x**3

def wavycubic(x):
    return np.piecewise(x, [x < 6, x > 5], [lambda y : y**3, lambda y : 250 - (10-y)**3])

# test knot coordinates and function values
ts_cubic = np.array([1, 2, 3, 4, 5])
ys_cubic = cubic(ts_cubic)

ts_wavycubic = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
ys_wavycubic = wavycubic(ts_wavycubic)

ts_irregular = np.array([1, 3, 4, 7, 10])
ys_irregular = cubic(ts_irregular)

# test splines
spline_cubic = CubicSpline(ts_cubic, ys_cubic)
spline_cubic.plotSpline(scatter=True, func=cubic, fname='PPEplots/cubictestplots/cubictest_cubic')

spline_wavycubic = CubicSpline(ts_wavycubic, ys_wavycubic)
spline_wavycubic.plotSpline(scatter=True, func=wavycubic, fname='PPEplots/cubictestplots/cubictest_wavycubic')

spline_irregular = CubicSpline(ts_irregular, ys_irregular)
spline_irregular.plotSpline(scatter=True, func=cubic, fname='PPEplots/cubictestplots/cubictest_irregular')
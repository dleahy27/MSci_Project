#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from bicubicspline import BicubicSpline

# latex font configs
plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['text.usetex'] = True


# test functions

## simple trig function
def trig(x, y):
    return np.cos(x) * np.sin(y) + 1.5

## edge and corner second derivatives of simple trig function
def deriv_trig(x, y, param='xx'):

    if param == 'xx':
        return - np.cos(x) * np.sin(y)
    
    elif param == 'yy':
        return - np.cos(x) * np.sin(y)
    
    elif param == 'xxyy':
        return np.cos(x) * np.sin(y)

## modulated trig function
## cos and sin arguments such that there is a decent amount
## of oscillation; offset so that the relative differences don't
## cause any zero division errors; added Gaussian to suppress the
## function at the edges so that it goes to zero and the imposed
## zero curvature is natural and causes minimal errors at the
## boundaries
def gauss(x, y):
    return (np.cos(3*y) * np.sin(3*x) + 1.5) * np.exp(-((x-3)**2 + (y-3)**2) / 6)


# constructors

## fill 2D array with function values at each input knot
def constructor(xs, ys, func):

    zs = np.empty([ys.size, xs.size])

    for iy, y in enumerate(ys):
        for ix, x in enumerate(xs):
            zs[iy, ix] = func(x, y)

    return zs

## fill array with edge and corner derivatives
def constructor_deriv(xs, ys, dfunc, edge='left', param='x'):

    if edge == 'bottom':
        
        dzs = np.empty([xs.size])
        for ix, x in enumerate(xs):
            dzs[ix] = dfunc(x, ys[0], param=param)
        return dzs
    
    elif edge == 'top':
        
        dzs = np.empty([xs.size])
        for ix, x in enumerate(xs):
            dzs[ix] = dfunc(x, ys[-1], param=param)
        return dzs

    elif edge == 'left':

        dzs = np.empty([ys.size])
        for iy, y in enumerate(ys):
            dzs[iy] = dfunc(xs[0], y, param=param)
        return dzs
    
    elif edge == 'right':

        dzs = np.empty([ys.size])
        for iy, y in enumerate(ys):
            dzs[iy] = dfunc(xs[-1], y, param=param)
        return dzs
    
    elif edge == 'corners':
        dzs = np.empty([4])
        dzs[0] = dfunc(xs[0], ys[0], param=param)
        dzs[1] = dfunc(xs[-1], ys[0], param=param)
        dzs[2] = dfunc(xs[-1], ys[-1], param=param)
        dzs[3] = dfunc(xs[0], ys[-1], param=param)
        return dzs


# create a control plot using true function values
def controlPlot(xs, ys, func, num=50, scatter=False, fname='PPEplots/bicubictestplots/control'):

    xs_plot = np.linspace(xs[0], xs[-1], num=num)
    ys_plot = np.linspace(ys[0], ys[-1], num=num)

    zs_plot = constructor(xs_plot, ys_plot, func)

    Xs_plot, Ys_plot = np.meshgrid(xs_plot, ys_plot)

    if scatter:
        Xs, Ys = np.meshgrid(xs, ys)


    plt.figure(figsize=(10, 8))

    plt.contourf(Xs_plot, Ys_plot, zs_plot, levels=12, cmap='plasma')

    cbar = plt.colorbar()
    cbar.ax.set_ylabel('$z$', fontsize=20, labelpad=15)

    if scatter:
        plt.scatter(Xs, Ys, color='white')

    plt.xlabel('$x$', fontsize=20, labelpad=15)
    plt.ylabel('$y$', fontsize=20, labelpad=15)
    plt.savefig(fname=fname+'.pdf')


# test coordinates

## test knot coordinates
coords_trig = np.array(np.linspace(-1, 7, 20))
coords_gauss = np.hstack(([-1.01,], np.linspace(-1, 7, 20), [7.01,]))
coords_irregular = np.array([-1, -0.5, 0, 1, 2, 2.1, 2.7, 3.1, 3.9, 4, 4.2, 5.4, 5.45, 6.1, 7])

xs_trig = coords_trig
ys_trig = coords_trig

xs_gauss = coords_gauss
ys_gauss = coords_gauss

xs_irregular = coords_irregular
ys_irregular = coords_irregular

## function values at test knots
zs_trig = constructor(xs_trig, ys_trig, trig)
zs_gauss = constructor(xs_gauss, ys_gauss, gauss)
zs_irregular = constructor(xs_irregular, ys_irregular, gauss)

## edge derivatives for simple trig spline
zs_xx_left_trig = constructor_deriv(xs_trig, ys_trig, deriv_trig, edge='left', param='xx')
zs_xx_right_trig = constructor_deriv(xs_trig, ys_trig, deriv_trig, edge='right', param='xx')
zs_yy_bottom_trig = constructor_deriv(xs_trig, ys_trig, deriv_trig, edge='bottom', param='yy')
zs_yy_top_trig = constructor_deriv(xs_trig, ys_trig, deriv_trig, edge='top', param='yy')
zs_xxyy_corners_trig = constructor_deriv(xs_trig, ys_trig, deriv_trig, edge='corners', param='xxyy')


# test splines

bispl_gauss = BicubicSpline(xs_gauss, ys_gauss, zs_gauss)
bispl_gauss.plotSpline(scatter=True, func=gauss, fname='pytest')
controlPlot(xs_gauss, ys_gauss, gauss, fname='pytest_cntrl')

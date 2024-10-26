#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from particle import Particle

## for LaTeX font on plots
## first need to install latex packages
#[in terminal] sudo apt install -q cm-super dvipng texlive-latex-extra texlive-latex-recommended
plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['text.usetex'] = True



def evaluatePatchPDF(Xs, Qs, fls, func):
    '''
    Evaluates pdf at given x, Q/Q^2 values, for given flavours fl
    and returns a filled array of function values.

    params: Xs:   meshgrid of x values
            Qs:   meshgrid of Q/Q^2 values
            fls:  pdf flavour(s) (pid)
            func: lhapdf method (pdf.xfxQ or pdf.xfxQ2)
    '''

    # evaluate pdf over meshgrid
    if type(fls) == np.ndarray:

        xfs = np.empty([fls.size, Xs[0].size, Qs[0].size])
        for ifl, fl in enumerate(fls):
            for ix in range(Xs[0].size):
                for iq in range(Qs[0].size):
                    xfs[ifl, ix, iq] = func(fl, Xs[ix, iq], Qs[ix, iq])
    
    else:

        xfs = np.empty([Xs[0].size, Qs[0].size])
        for ix in range(Xs[0].size):
            for iq in range(Qs[0].size):
                xfs[ix, iq] = func(fls, Xs[ix, iq], Qs[ix, iq])

    # return pdf values
    return xfs



def plotPatchPDF1DAllFlavours(pdf, func, axis=0, index_x=0, index_q=0, xstart=-3, xend=0, qstart=0, qend=8, num=50, fname='PPEplots/patchpdfplots/patchpdf_1D', axes=['$x$', '$Q^2$', '$xf(x, Q^2)$']):
    '''
    Creates 1D plot of pdf for all flavours based on patch algorithm.

    params: pdf:       lhapdf object
            func:      lhapdf method (pdf.xfxQ or pdf.xfxQ2)

    kwargs: axis:      axis to plot along (0:x; 1:Q/Q^2)
            index_x:   x index to plot for if axis=1
            index_q:   Q/Q^2 index to plot for if axis=0
            xstart:    start of x range (for axis=0)
            xend:      end of x range (for axis=0)
            qstart:    start of Q/Q^2 range (for axis=1)
            qend:      end of Q/Q^2 range (for axis=1)
            num:       number of points to sample within given range
            fname:     file path for plots, without the file type
            axes:      name of the plot axes; default 'x', 'Q', 'xf(x, Q)'
    '''

    # store all flavours in array
    fls = np.array([pid for pid in pdf.flavors()])

    # coordinates
    mesh_xs = np.array(np.logspace(xstart, xend, num=num))
    mesh_qs = np.array(np.logspace(qstart, qend, num=num))

    # create a meshgrid
    Xs, Qs = np.meshgrid(mesh_xs, mesh_qs)
    
    # evaluate pdf of each flavour at each input coordinate
    xfs = evaluatePatchPDF(Xs, Qs, fls, func)

    # initialise figure
    plt.figure(figsize=(10, 8))

    # set x, Q/Q^2 to plot at
    ix = index_x
    iq = index_q

    # list of nice colours for plotting
    colors = ['r', 'orange', 'yellow', 'green', 'cyan', 'b', 'violet', 'pink', 'brown', 'gray', 'k']

    # loop over all flavours
    for ifl, fl in enumerate(fls):
        
        # get name of particle for legend
        p = Particle.from_pdgid(fl)
        p_name = getattr(p, 'latex_name')

        # plot for given flavour along specified axis, distinguishing gluons
        if axis == 0:

            if fl == 21:
                plt.plot(mesh_xs, xfs[ifl, iq], color=colors[ifl], linestyle='--', label='${0}$'.format(p_name))
            else:
                plt.plot(mesh_xs, xfs[ifl, iq], color=colors[ifl], label='${0}$'.format(p_name))
        
        elif axis == 1:

            if fl == 21:
                plt.plot(mesh_qs, xfs[ifl, :, ix], color=colors[ifl], linestyle='--', label='${0}$'.format(p_name))
            else:
                plt.plot(mesh_qs, xfs[ifl, :, ix], color=colors[ifl], label='${0}$'.format(p_name))
        
        plt.xscale('log')
        plt.ylim(0, 1)

    ## plot configs
    plt.xlabel(axes[axis], fontsize=20, labelpad=15)
    plt.ylabel(axes[2], fontsize=20, labelpad=15)
    plt.legend(fontsize=15)
    plt.savefig(fname=fname+'.pdf')



def plotPatchPDF(fl, func, xstart=-3, xend=0, qstart=0, qend=8, num=50, fname='PPEplots/patchpdfplots/patchpdf', axes=['$x$', '$Q^2$', '$xf(x, Q^2)$']):
    '''
    Creates plot of pdf based on patch algorithm.

    params: fl:        pdf flavour (pid)
            func:      lhapdf method (pdf.xfxQ or pdf.xfxQ2)
    
    kwargs: xstart:    start of x range
            xend:      end of x range
            qstart:    start of Q/Q^2 range
            qend:      end of Q/Q^2 range
            num:       number of points to sample within given range
            fname:     file path for plots, without the file type
            axes:      name of the plot axes; default 'x', 'Q', 'xf(x, Q)'
    '''

    # cooridnates
    mesh_xs = np.array(np.logspace(xstart, xend, num=num))
    mesh_qs = np.array(np.logspace(qstart, qend, num=num))

    # create a meshgrid
    Xs, Qs = np.meshgrid(mesh_xs, mesh_qs)
    
    # evaluate pdf of specified flavour at each input coordinate
    xfs = evaluatePatchPDF(Xs, Qs, fl, func)

    ## initialise figure
    plt.figure(figsize=(10, 8))

    ## plot for given flavour
    plt.contourf(Xs, Qs, xfs, levels=12, cmap='plasma')
    plt.xscale('log')
    plt.yscale('log')

    ## set colorbar
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(axes[2], fontsize=20, labelpad=15)

    ## plot configs
    plt.xlabel(axes[0], fontsize=20, labelpad=15)
    plt.ylabel(axes[1], fontsize=20, labelpad=15)
    plt.savefig(fname=fname+'.pdf')



def plotPatchPDFDerivative(fl, func, order=1, xstart=-3, xend=0, qstart=0, qend=8, num=50, fname='PPEplots/patchpdfplots/patchpdf', axes=['$x$', '$Q^2$', '$xf(x, Q^2)$'], plot1d=False, index_x=0, index_y=0):
    '''
    Creates plot of pdf derivative based on patch algorithm.

    params: fl:        pdf flavour (pid)
            func:      lhapdf method (pdf.xfxQ or pdf.xfxQ2)

    kwargs: order:     order of derivative (1 or 2)
            xstart:    start of x range
            xend:      end of x range
            qstart:    start of Q/Q^2 range
            qend:      end of Q/Q^2 range
            num:       number of points to sample within given range
            fname:     file path for plots, without the file type
            axes:      name of the plot axes; default 'x', 'Q', 'xf(x, Q)'
            plot1d:    create additional 1D plots along given indices 
                       in x and y
            index_x:   plot spline y derivative along x=index_x
            index_y:   plot spline x derivative along y=index_y
    '''

    # coordinates
    mesh_xs = np.array(np.logspace(xstart, xend, num=num))
    mesh_qs = np.array(np.logspace(qstart, qend, num=num))

    # create a meshgrid
    Xs, Qs = np.meshgrid(mesh_xs, mesh_qs)

    # evaluate pdf of specified flavour at each input coordinate
    xfs = evaluatePatchPDF(Xs, Qs, fl, func)

    # calculate second derivatives of spline at each xy coordinate
    if order == 1:
        xfs_xx = np.gradient(xfs, mesh_xs, axis=1)
        xfs_yy = np.gradient(xfs, mesh_qs, axis=0)
    
    elif order == 2:
        xfs_xx = np.gradient(np.gradient(xfs, mesh_xs, axis=1), mesh_xs, axis=1)
        xfs_yy = np.gradient(np.gradient(xfs, mesh_qs, axis=0), mesh_qs, axis=0)

    else: raise ValueError('Derivatives available only up to 2nd order.')


    # heatmap of second derivatives along x

    ## initialise figure
    plt.figure(figsize=(10, 8))

    ## plot for given flavour
    plt.contourf(Xs, Qs, xfs_xx, levels=12, cmap='plasma')
    plt.xscale('log')
    plt.yscale('log')

    ## set colorbar
    cbar = plt.colorbar()
    if order == 1:
        cbar.ax.set_ylabel('$\mathrm{d} \;$'+axes[2], fontsize=20, labelpad=15)
    elif order == 2:
        cbar.ax.set_ylabel('$\mathrm{d}^2 \;$'+axes[2], fontsize=20, labelpad=15)

    ## plot configs
    plt.xlabel(axes[0], fontsize=20, labelpad=15)
    plt.ylabel(axes[1], fontsize=20, labelpad=15)

    if order == 1: 
        plt.savefig(fname=fname+'_dx.pdf')
    elif order == 2: 
        plt.savefig(fname=fname+'_d2x.pdf')


    # heatmap of second derivatives along y

    ## initialise figure
    plt.figure(figsize=(10, 8))

    ## plot for given flavour
    plt.contourf(Xs, Qs, xfs_yy, levels=12, cmap='plasma')
    plt.xscale('log')
    plt.yscale('log')

    ## set colorbar
    cbar = plt.colorbar()
    if order == 1:
        cbar.ax.set_ylabel('$\mathrm{d} \;$'+axes[2], fontsize=20, labelpad=15)
    elif order == 2:
        cbar.ax.set_ylabel('$\mathrm{d}^2 \;$'+axes[2], fontsize=20, labelpad=15)

    ## plot configs
    plt.xlabel(axes[0], fontsize=20, labelpad=15)
    plt.ylabel(axes[1], fontsize=20, labelpad=15)

    if order == 1:
        plt.savefig(fname=fname+'_dy.pdf')
    elif order == 2:
        plt.savefig(fname=fname+'_d2y.pdf')


    if plot1d:

        # 1D plot of the second x derivative along an arbitrary grid line

        # relevant slice of derivative values
        xfs_xx_1D = xfs_xx[index_y]
        
        # initialise figure
        plt.figure(figsize=(10, 8))

        # plot spline derivative in x
        plt.plot(mesh_xs, xfs_xx_1D, color='r')

        # plot configs
        plt.xlabel(axes[0], fontsize=20)
        plt.ylabel('$\mathrm{d}^2 \;$'+axes[2], fontsize=20)
        plt.xscale('log')
        
        if order == 1: 
            plt.savefig(fname=fname+'_dx_1D.pdf')
        elif order == 2: 
            plt.savefig(fname=fname+'_d2x_1D.pdf')


        # 1D plot of the second y derivative along an arbitrary grid line

        # relevant slice of derivative values
        xfs_yy_1D = xfs_yy[:, index_x]
        
        # initialise figure
        plt.figure(figsize=(10, 8))

        # plot spline derivative in x
        plt.plot(mesh_qs, xfs_yy_1D, color='r')

        # plot configs
        plt.xlabel(axes[1], fontsize=20)
        plt.ylabel('$\mathrm{d}^2 \;$'+axes[2], fontsize=20)
        plt.xscale('log')

        if order == 1: 
            plt.savefig(fname=fname+'_dy_1D.pdf')
        elif order == 2: 
            plt.savefig(fname=fname+'_d2y_1D.pdf')
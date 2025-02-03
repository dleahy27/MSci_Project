import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats.mstats import mquantiles
from scipy.stats import multivariate_normal

import glob

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'

# Plot data
title_size = 14
cbar_size = 14
axes_size = 14

def dynamicColorbar(z):
    # 5% and 95% quantiles
    minmax = mquantiles(z, prob=[0.1,0.9])
    # colorbar norm, levels and ticks set to each case
    # Bigger negative value
    if ( np.abs(minmax[0]) > np.abs(minmax[1]) ):
        order = int(np.floor(np.log10(np.abs(minmax[0]))))
        
        norm = TwoSlopeNorm(vmin=minmax[0], vcenter=0, vmax=-minmax[0])
        levels = np.linspace(minmax[0], -minmax[0], 100, endpoint=True)
        ticks = np.linspace(minmax[0], -minmax[0], 9, endpoint=True)
        ticks = [np.round(tick, decimals=-order+2) for tick in ticks]
        
    # Bigger positive value    
    elif ( np.abs(minmax[1]) > np.abs(minmax[0]) ):
        # colorbar norm, levels and ticks set to these
        order = int(np.floor(np.log(np.abs(minmax[1]))))
        norm = TwoSlopeNorm(vmin=-minmax[1], vcenter=0, vmax=minmax[1])
        levels = np.linspace(-minmax[1], minmax[1], 100, endpoint=True)
        ticks = np.linspace(-minmax[1], minmax[1], 9, endpoint=True)
        ticks = [np.round(tick, decimals=-order+2) for tick in ticks]      
    else:
        order = 1
        minmax[0] = -0.1
        minmax[1] = 0.1
        norm = TwoSlopeNorm(vmin=minmax[0], vcenter=0, vmax=minmax[1])
        levels = np.linspace(minmax[0], minmax[1], 100, endpoint=True)
        ticks = np.linspace(minmax[0], minmax[1], 9, endpoint=True)
        ticks = [np.round(tick, decimals=-order+2) for tick in ticks]
        
    return norm, levels, ticks

def dataInit(name, plot_nx, plot_ny, x_axes, y_axes, i, j, k, l):
    nx,ny,x,y,z,z_num = np.loadtxt(name, dtype='float64', delimiter=',', unpack=True, skiprows=1)
    nx = nx[0]
    ny = ny[0]
    
    ratio = (z_num/z - 1)*100
    
    median_error = np.median(ratio.flatten())
    
    scaling = False
    if scaling:
        # X Scaling Plots
        if ((plot_nx == i).any() and j == int(2*np.max(plot_ny)/3) ):
            plotGridx(x_axes, k, x, y, z_num, nx, median_error)
            k += 1
            return nx,ny,median_error, k, l
        
        # Y Scaling Plots
        if ((plot_ny == j).any() and i == int(2*np.max(plot_nx)/3)):
            plotGridx(y_axes, l, x, y, z_num, ny, median_error)
            l += 1
            return nx,ny,median_error, k, l
    
    return nx,ny,median_error, k,l

# DO THE SAME THING BUT PLOT Y GRID. X DOESNT SEEM TO CHANGE MUCH
def plotGridx(axes,k,x,y,z, nx, median):
    # imported data is 1D columns
    sizex = 51
    sizey = 51
    X = x.reshape((sizey,sizex))
    Y = y.reshape((sizey,sizex))
    Z = z.reshape((sizey,sizex))
    
    X = X[:,]
    Y = Y[:,]
    Z = Z[:,]

    x = X.flatten()
    y = Y.flatten()
    z = Z.flatten()
    
    data = np.vstack((x, y, z))
    
    cov = np.cov(data)
    cov_inv = np.linalg.inv(cov)
    
    mean = np.mean(data, axis=1)
    
    grid_points = np.vstack([x, y, z]).T
    
    # Compute Mahalanobis distance for each point
    mahalanobis_distances = np.array([
        np.dot(np.dot((point - mean), cov_inv), (point - mean).T) for point in grid_points
    ])

    # Reshape the Mahalanobis distances back to the shape of the grid
    mahalanobis_distances = mahalanobis_distances.reshape(X.shape)
    
    plot1 = axes[k].pcolormesh(X,Y,Z)
    cb1 = plt.colorbar(plot1, ax=axes[k])
    if (k == 3 or k == 7):
        cb1.set_label(r"Relative Error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    
    cont = axes[k].contour(X, Y, mahalanobis_distances, levels=[1], colors='m', linewidths=2)
    CL = plt.clabel(cont, colors = 'm', fontsize=10, fmt=r"1$\sigma$")
    
    # Convert the covariance matrix into LaTeX format
    cov_latex = r"$\Sigma = \begin{bmatrix} " + \
                    r" %.2g & %.2g & %.2g \\" % (cov[0, 0], cov[0, 1], cov[0, 2], ) + \
                    r" %.2g & %.2g & %.2g \\" % (cov[1, 0], cov[1, 1], cov[1, 2], ) + \
                    r" %.2g & %.2g & %.2g " % (cov[2, 0], cov[2, 1], cov[2, 2], ) +r"\end{bmatrix}$"
                   
    textstr = r'$\mu_{1/2}=%.4g\newline$' % (median, ) + cov_latex
    
    axes[k].text(0, 1.01, textstr, transform=axes[k].transAxes,
            ha='left',
            va='bottom',
            fontsize=8, color='black', 
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    
    axes[k].text(0.85, 0.95, r'$n_x = %.2g$' % (nx, ), transform=axes[k].transAxes,
                ha='center',
                va='center',
                fontsize = title_size-2, color='black', 
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    
def plotGridy(axes,l,x,y,z, ny, median):
    # imported data is 1D columns
    sizex = 51
    sizey = 51
    X = x.reshape((sizey,sizex))
    Y = y.reshape((sizey,sizex))
    Z = z.reshape((sizey,sizex))

    x = X.flatten()
    y = Y.flatten()
    z = Z.flatten()
    
    data = np.vstack((x, y, z))
    
    cov = np.cov(data)
    cov_inv = np.linalg.inv(cov)
    
    mean = np.mean(data, axis=1)
    
    grid_points = np.vstack([x, y, z]).T
    
    # Compute Mahalanobis distance for each point
    mahalanobis_distances = np.array([
        np.dot(np.dot((point - mean), cov_inv), (point - mean).T) for point in grid_points
    ])

    # Reshape the Mahalanobis distances back to the shape of the grid
    mahalanobis_distances = mahalanobis_distances.reshape(X.shape)
    
    plot1 = axes[l].pcolormesh(X,Y,Z)
    cb1 = plt.colorbar(plot1, ax=axes[l])
    if (l == 3 or l == 7):
        cb1.set_label(r"Relative Error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    
    cont = axes[l].contour(X, Y, mahalanobis_distances, levels=[1], colors='m', linewidths=2)
    CL = plt.clabel(cont, colors = 'm', fontsize=10, fmt=r"1$\sigma$")
    
    # Convert the covariance matrix into LaTeX format
    cov_latex = r"$\Sigma = \begin{bmatrix} " + \
                    r" %.2g & %.2g & %.2g \\" % (cov[0, 0], cov[0, 1], cov[0, 2], ) + \
                    r" %.2g & %.2g & %.2g \\" % (cov[1, 0], cov[1, 1], cov[1, 2], ) + \
                    r" %.2g & %.2g & %.2g " % (cov[2, 0], cov[2, 1], cov[2, 2], ) +r"\end{bmatrix}$"
                   
    textstr = r'$\mu_{1/2}=%.4g\newline$' % (median, ) + cov_latex
    
    axes[l].text(0, 1.01, textstr, transform=axes[l].transAxes,
            ha='left',
            va='bottom',
            fontsize=8, color='black', 
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    
    axes[l].text(0.85, 0.95, r'$n_y = %.2g$' % (ny, ), transform=axes[l].transAxes,
                ha='center',
                va='center',
                fontsize = title_size-2, color='black', 
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    
if __name__ == '__main__':
    names = sorted(glob.glob("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/grid_n_test/pdf_grid_*.csv"), key=len)
    xs = 20
    ys = 20
    
    plot_xs = np.linspace(0, xs-1, 8, dtype = int)
    plot_ys = np.linspace(0, xs-1, 8, dtype = int)
    
    medians = np.empty((xs,ys))
    nxs = np.empty(xs)
    nys = np.empty(ys)
    
    figx = plt.figure(figsize=(16,12))
    axesx = figx.subplots(2,4).flatten()
    
    figy = plt.figure(figsize=(16,12))
    axesy = figy.subplots(2,4).flatten()
    
    k = 0
    l = 0
    for i in range(xs):
        for j in range(ys):
            nx, ny, median_error, k, l = dataInit(names[j + 20*i], plot_xs, plot_ys, axesx, axesy, i, j, k, l)
            medians[i][j] = median_error
            nxs[i] = nx
            nys[j] = ny

    filenamex = "../outputs/grid_n_test/testx.pdf"
    filenamey = "../outputs/grid_n_test/testy.pdf"
    
    figx.supxlabel(r'$\log(x)$', fontsize=axes_size+4)
    figx.supylabel(r'$\log(Q)$', fontsize=axes_size+4)
    figx.suptitle(r'Relative Error Scaling vs $n_x$ at $n_y = %.2g$ ' % (nys[int(2*np.max(plot_xs)/3)], ), fontsize=title_size+4)
    
    figy.supxlabel(r'$\log(x)$', fontsize=axes_size+4)
    figy.supylabel(r'$\log(Q)$', fontsize=axes_size+4)
    figy.suptitle(r'Relative Error Scaling vs $n_y$ at $n_x = %.2g$ ' % (nxs[int(2*np.max(plot_ys)/3)], ), fontsize=title_size+4)
    
    figx.savefig(filenamex,dpi=200, format="pdf")
    figy.savefig(filenamey,dpi=200, format="pdf")
    
    fig  = plt.figure(figsize=(8,8))
    axes = fig.subplots(1,1)
    
    plot1 = axes.pcolormesh(nxs,nys,medians)
    # axes[0].set_title(r"Median Relative Error", fontsize=title_size)
    axes.set_xlabel(r"X Grid Points")
    axes.set_ylabel(r"Y Grid Points")
    
    cb1 = plt.colorbar(plot1, ax=axes)
    
    cb1.set_label(r"Median Relative Error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    
    fig.savefig("../outputs/grid_n_test/Error_Scale.pdf", dpi=200, format="pdf")
    
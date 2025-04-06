import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats.mstats import mquantiles
from scipy.stats import multivariate_normal
from scipy.optimize import curve_fit

import glob

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'

# Plot data
title_size = 10
cbar_size = 10
axes_size = 10

def error_func_log(m,c,x):
    # straight line as it does not handle powers well
    return m*x + c

def error_func(a,b,x):
    # straight line as it does not handle powers well
    return a*(x**b)

def dataInitMedian(name):
    nx,ny,x,y,z,z_num = np.loadtxt(name, dtype='float64', delimiter=',', unpack=True, skiprows=1)
    nx = nx[0]
    ny = ny[0]
    
    ratio = (z_num/z - 1)*100
    ratio = ratio[ratio != 0]
    # q10, q20, q30, q40 = np.percentile(ratio, [10 ,20, 30, 40])
    
    # return nx, ny, median_error, q10, q20, q30, q40
    q16, q25, median, q75, q84 = np.percentile(ratio, [16, 25, 50, 75, 84], method="averaged_inverted_cdf")
    
    return nx, ny, median, q16, q25, q75, q84

def dataInitMean(name):
    nx,ny,x,y,z,z_num = np.loadtxt(name, dtype='float64', delimiter=',', unpack=True, skiprows=1)
    nx = nx[0]
    ny = ny[0]
    
    ratio = (z_num/z - 1)*100
    
    mean = np.mean(ratio.flatten())

    std_dev = np.std(ratio.flatten())
    
    return nx, ny, mean, std_dev

if __name__ == '__main__':
    names = sorted(glob.glob("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/grid_x_test/pdf_xgrid_*.csv"), key=len)

    xs = len(names)
    
    means = np.empty(xs)
    std_devs = np.empty(xs)
    
    medians = np.empty(xs)
    q25s = np.empty(xs)
    q75s = np.empty(xs)
    q16s = np.empty(xs)
    q84s = np.empty(xs)
    nxs = np.empty(xs)

    for i in range(xs):
            nx, ny, median_error, q16, q25, q75, q84 = dataInitMedian(names[i])
            nx, ny, mean, std_dev = dataInitMean(names[i])
            
            medians[i] = median_error
            q16s[i] = q16
            q25s[i] = q25
            q75s[i] = q75
            q84s[i] = q84
            
            means[i] = mean
            std_devs[i] = std_dev
            
            nxs[i] = nx


    abs_means = np.abs(means)
    test_xs = np.linspace(nxs[0],nxs[-1], 1000)
    
    ########################################################################################################
    filename = "../outputs/grid_x_test/medians.pdf"
    
    fig = plt.figure(figsize=(8,4))
    ax = fig.subplots(1, 2).flatten()
    
    plot = ax[0].plot(nxs, medians, marker = ".", linestyle="", label = r"Median")
    # ax[0].set_xlabel(r"(a)")
    # ax[0].plot(test_xs, y1, color=plot[0].get_color())
    
    plot = ax[0].plot(nxs, q16s, marker = ".", linestyle="", label = r"Median - $\sigma$")
    # ax[0].plot(test_xs, y2, color=plot[0].get_color())
    
    plot = ax[0].plot(nxs, q84s, marker = ".", linestyle="", label = r"Median + $\sigma$")
    # ax[0].plot(test_xs, y3, color=plot[0].get_color())
    
    ax[0].set_xlabel(r"(a)", fontsize=axes_size+6)
    # ax[0].set_title(r"Median")
    ax[0].legend()
    
    plot = ax[1].plot(nxs[nxs>40], medians[nxs>40], marker = ".", linestyle="", label = r"Median")
    # ax[0].plot(test_xs, y1, color=plot[0].get_color())
    
    plot = ax[1].plot(nxs[nxs>40], q16s[nxs>40], marker = ".", linestyle="", label = r"Median - $\sigma$")
    # ax[0].plot(test_xs, y2, color=plot[0].get_color())
    ax[1].set_xlabel(r"(b) $n_u \geq 40$", fontsize = axes_size+6)
    plot = ax[1].plot(nxs[nxs>40], q84s[nxs>40], marker = ".", linestyle="", label = r"Median + $\sigma$")
    # ax[0].plot(test_xs, y3, color=plot[0].get_color())
    #ax[1].set_xlabel(r"(b)")
    
    fig.supxlabel(r'$n_u$', fontsize=axes_size+6)
    fig.supylabel(r'$\Delta xf(u,v)(\%)$', fontsize=axes_size+6)
    fig.suptitle(r'$n_{v} = %.2g$ ' % (ny), fontsize=title_size+6)
    
    fig.tight_layout()
    fig.savefig(filename,dpi=200, format="pdf")
    ######################################################################################################
    filename = "../outputs/grid_x_test/mean.pdf"
    
    fig = plt.figure(figsize=(12,8))
    ax = fig.subplots(1, 2).flatten()
    
    ax[0].plot(nxs, np.abs(means), marker = ".")
    ax[0].set_title(r"Mean")
    
    ax[1].plot(nxs, np.abs(std_devs))
    ax[1].set_title(r"Standard Deviation")
    
    fig.supxlabel(r'$n_x$', fontsize=axes_size+4)
    fig.supylabel(r'$|\Delta_{rel}|(\%)$', fontsize=axes_size+4)
    fig.suptitle(r'$\Delta_{rel}$ Scaling at $n_y = %.2g$ ' % (ny), fontsize=title_size+4)
    
    fig.savefig(filename,dpi=200, format="pdf")
   ###################################################################################################### 
    filename = "../outputs/grid_x_test/mean_cuts.pdf"
    
    fig = plt.figure(figsize=(12,8))
    ax = fig.subplots(1, 2).flatten()
    
    ax[0].plot(nxs[abs_means<10E8], means[abs_means<10E8], marker = ".")
    ax[0].set_title(r"Mean")
    
    ax[1].plot(nxs[abs_means<10E8], std_devs[abs_means<10E8])
    ax[1].set_title(r"Standard Deviation")
    
    fig.supxlabel(r'$n_x$', fontsize=axes_size+4)
    fig.supylabel(r'$\Delta_{rel}(\%)$', fontsize=axes_size+4)
    fig.suptitle(r'$\Delta_{rel}$ Scaling at $n_y = %.2g$ ' % (ny), fontsize=title_size+4)
    
    fig.savefig(filename,dpi=200, format="pdf")
    #######################################################################################################
    filename = "../outputs/grid_x_test/test_2.pdf"
    
    q16s = q16s[nxs>40]
    q25s = q25s[nxs>40]
    q75s = q75s[nxs>40]
    q84s = q84s[nxs>40]
    
    medians = medians[nxs>40]
    means = means[nxs>40]
    abs_means = abs_means[nxs>40]
    std_devs = std_devs[nxs>40]
    
    nxs = nxs[nxs>40]
    test_xs = test_xs[test_xs>40]
    
    fig = plt.figure(figsize=(12,8))
    ax = fig.subplots(1, 2).flatten()
    
    plot = ax[0].plot(nxs, medians, marker = ".", linestyle="", label = r"Median")
    # ax[0].plot(test_xs, y1, color=plot[0].get_color())

    plot = ax[0].plot(nxs, q16s, marker = ".", linestyle="", label = r"Median - $\sigma$")
    # ax[0].plot(test_xs, y2, color=plot[0].get_color())

    plot = ax[0].plot(nxs, q84s, marker = ".", linestyle="", label = r"Median + $\sigma$")
    # ax[0].plot(test_xs, y3, color=plot[0].get_color())

    
    ax[0].set_title(r"Median")
    ax[0].legend()
    
    ax[1].plot(nxs, np.abs(q25s - q75s))
    ax[1].set_title(r"InterQuartile Range")
    
    fig.supxlabel(r'$n_x$', fontsize=axes_size+4)
    fig.supylabel(r'$|\Delta_{rel}|(\%)$', fontsize=axes_size+4)
    fig.suptitle(r'$\Delta_{rel}$ Scaling at $n_y = %.2g$ ' % (ny), fontsize=title_size+4)
    
    fig.savefig(filename,dpi=200, format="pdf")
    ########################################################################################################
    filename = "../outputs/grid_x_test/test_mean2.pdf"
    
    fig = plt.figure(figsize=(12,8))
    ax = fig.subplots(1, 2).flatten()
    
    ax[0].plot(nxs, np.abs(means), marker = ".")
    ax[0].set_title(r"Mean")
    
    ax[1].plot(nxs, np.abs(std_devs))
    ax[1].set_title(r"Standard Deviation")
    
    fig.supxlabel(r'$n_x$', fontsize=axes_size+4)
    fig.supylabel(r'$|\Delta_{rel}|(\%)$', fontsize=axes_size+4)
    fig.suptitle(r'$\Delta_{rel}$ Scaling at $n_y = %.2g$ ' % (ny), fontsize=title_size+4)
    
    fig.savefig(filename,dpi=200, format="pdf")
    ###################################################################################################
    filename = "../outputs/grid_x_test/test_mean_cuts_2.pdf"
    
    fig = plt.figure(figsize=(12,8))
    ax = fig.subplots(1, 2).flatten()
    
    ax[0].plot(nxs[abs_means<10E8], means[abs_means<10E8], marker = ".")
    ax[0].set_title(r"Mean")
    
    ax[1].plot(nxs[abs_means<10E8], std_devs[abs_means<10E8])
    ax[1].set_title(r"Standard Deviation")
    
    fig.supxlabel(r'$n_x$', fontsize=axes_size+4)
    fig.supylabel(r'$|\Delta_{rel}|(\%)$', fontsize=axes_size+4)
    fig.suptitle(r'$\Delta_{rel}$ Scaling at $n_y = %.2g$ ' % (ny), fontsize=title_size+4)
    
    fig.savefig(filename,dpi=200, format="pdf")

    
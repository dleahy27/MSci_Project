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

def dataInitMedian(name):
    nx,ny,x,y,z,z_num = np.loadtxt(name, dtype='float64', delimiter=',', unpack=True, skiprows=1)
    nx = nx[0]
    ny = ny[0]
    
    ratio = (z_num/z - 1)*100
    
    median_error = np.median(ratio.flatten())

    # q10, q20, q30, q40 = np.percentile(ratio, [10 ,20, 30, 40])
    
    # return nx, ny, median_error, q10, q20, q30, q40
    q25, q75 = np.percentile(ratio, [25, 75])
    
    return nx, ny, median_error, q25, q75

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
    
    xs = 100
    
    medians = np.empty(xs)
    q10s = np.empty(xs)
    q20s = np.empty(xs)
    q30s = np.empty(xs)
    q40s = np.empty(xs)
    q25s = np.empty(xs)
    q75s = np.empty(xs)
    nxs = np.empty(xs)

    for i in range(xs):
            # nx, ny, median_error, q10, q20, q30, q40 = dataInitMedian(names[i])
            nx, ny, median_error, q25, q75 = dataInitMedian(names[i])
            medians[i] = median_error
            # q10s[i] = q10
            # q20s[i] = q20
            # q30s[i] = q30
            # q40s[i] = q40
            q25s[i] = q25
            q75s[i] = q75
            
            nxs[i] = nx

    print(nxs)
    # print(std_devs)
    filename = "../outputs/grid_x_test/test.pdf"
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1, 2).flatten()
    
    ax[0].plot(nxs, np.abs(medians), marker = ".", label = r"Median")
    # ax[0].plot(nxs, np.abs(q10s), marker = ".", label = r"10\% Quartile")
    # ax[0].plot(nxs, np.abs(q20s), marker = ".", label = r"20\% Quartile")
    # ax[0].plot(nxs, np.abs(q30s), marker = ".", label = r"30\% Quartile")
    # ax[0].plot(nxs, np.abs(q40s), marker = ".", label = r"40\% Quartile")
    ax[0].plot(nxs, np.abs(q25s), marker = ".", label = r"$25\%$ Quartile")
    ax[0].plot(nxs, np.abs(q75s), marker = ".", label = r"$75\%$ Quartile")
    ax[0].legend()
    #ax.set_ylim([-0.05,0.4])
    
    ax[1].plot(nxs, np.abs(q25s - q75s))
    ax[1].set_title(r"Inter Quartile Range")
    
    fig.supxlabel(r'$n_x$', fontsize=axes_size+4)
    fig.supylabel(r'$|\epsilon_{rel}|(\%)$', fontsize=axes_size+4)
    fig.suptitle(r'$\epsilon_{rel}$ Scaling at $n_y = %.2g$ ' % (ny), fontsize=title_size+4)
    
    fig.savefig(filename,dpi=200, format="pdf")

    
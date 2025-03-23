import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats.mstats import mquantiles
from ct10nlo import CT10nlo

plt.rcParams["text.usetex"] = True

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
        order = int(np.floor(np.log10(np.abs(minmax[1]))))
        norm = TwoSlopeNorm(vmin=-minmax[1], vcenter=0, vmax=minmax[1])
        levels = np.linspace(-minmax[1], minmax[1], 100, endpoint=True)
        ticks = np.linspace(-minmax[1], minmax[1], 9, endpoint=True)
        ticks = [np.round(tick, decimals=-order+2) for tick in ticks]      
    else:
        order = 1
        minmax[0] = -1
        minmax[1] = 1
        norm = TwoSlopeNorm(vmin=minmax[0], vcenter=0, vmax=minmax[1])
        levels = np.linspace(minmax[0], minmax[1], 100, endpoint=True)
        ticks = np.linspace(minmax[0], minmax[1], 9, endpoint=True)
        ticks = [np.round(tick, decimals=-order+2) for tick in ticks]
        
    return norm, levels, ticks

if __name__=='__main__':
    # read in data
    xs,qs,xfs_new = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/CT10nlo.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)
    xs,qs,xfs_legacy = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/gitlhapdf/share/LHAPDF/legacy_up.dat", dtype='float', delimiter=' ', unpack=True)
    
    grid_x = CT10nlo.logx
    grid_y = CT10nlo.logQ2
    
    nx = 801
    ny = 901
    
    # Plot data
    title_size = 18
    cbar_size = 18
    axes_size = 18
    msize = 4
    
    xfs_legacy = xfs_legacy[xs>=1E-8]
    xs = xs[xs>=1E-8]
    
    xs = np.log10(xs[::ny])
    qs = np.log10(qs[:ny])
    
    grid_y = grid_y[grid_y>=np.min(qs)]
    
    xfs_legacy = xfs_legacy.reshape((nx,ny)).T
    xfs_new = xfs_new.reshape((ny,nx))
    diff = 100*(-xfs_legacy+xfs_new)/xfs_legacy
    
    XX,YY = np.meshgrid(grid_x,grid_y)
    
    # 2D plots
    fig = plt.figure(figsize=(12,6))
    axs = fig.subplots(1,2)

    plot1 = axs[0].pcolormesh(xs,qs,xfs_new, rasterized=True)
    axs[0].scatter(XX, YY, color="black", s=msize)
    axs[0].set_title(r"Bicubic Spline")
    
    plot2 = axs[1].pcolormesh(xs,qs,xfs_legacy, rasterized=True)
    axs[1].scatter(XX, YY, color="black", s=msize)
    axs[1].set_title(r"Cubic Spline + Local")
    
    
    cb1 = plt.colorbar(plot1, ax=axs[0], label=r"$xf(x,Q)$")
    cb2 = plt.colorbar(plot2, ax=axs[1], label=r"$xf(x,Q)$")
    cb1.set_label(r"$xf(x,Q)$", rotation=270, labelpad=20)
    cb2.set_label(r"$xf(x,Q)$", rotation=270, labelpad=20)

    fig.supxlabel(r'$\log_{10}x$')
    fig.supylabel(r'$\log_{10}Q^2$')
    
    fig.tight_layout()
    fig.savefig(r"../outputs/LHAPDF/up/LHAPDF_test.pdf", format="pdf")
    
    fig = plt.figure(figsize=(6,6))
    ax = fig.subplots(1,1)
    norm,levels,ticks = dynamicColorbar(diff) 
    plot = ax.pcolormesh(xs,qs,diff, norm=norm, cmap="seismic", rasterized=True)
    ax.scatter(XX, YY, color="black", s=msize)
    cb = plt.colorbar(plot, ticks=ticks, ax=ax)
    cb.set_label(r"$\Delta xf(x,Q^2) (\%)$", rotation=270,fontsize=cbar_size+4, labelpad=20)
    fig.supxlabel(r'$\log_{10}x$',fontsize=cbar_size+4)
    fig.supylabel(r'$\log_{10}Q^2$',fontsize=cbar_size+4)
    fig.tight_layout()
    fig.savefig(r"../outputs/LHAPDF/up/LHAPDF_diff.pdf", format="pdf")
    
    # 2D plots
    fig = plt.figure(figsize=(6,6))
    ax = fig.subplots(1,1)

    plot1 = ax.pcolormesh(xs,qs,xfs_new, rasterized=True)
    ax.scatter(XX, YY, color="black", s=msize)
    ax.set_xlabel(r'$\log_{10}x$',fontsize=cbar_size,)
    ax.set_ylabel(r'$\log_{10}Q^2$',fontsize=cbar_size,)
    cb1 = plt.colorbar(plot1, ax=ax)
    cb1.set_label(r"$xf(x,Q^2)$", rotation=270, fontsize=cbar_size, labelpad=20)

    fig.tight_layout()
    fig.savefig(r"../outputs/LHAPDF/up/bicubic_presentation.pdf", format="pdf")
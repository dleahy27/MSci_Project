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
    
if __name__ == '__main__':
    # read in data
    an_d1x, an_d1y, an_d2x, an_d2y = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/ct10_derivs.csv", dtype='float', delimiter=',', unpack=True, skiprows=1, usecols=[2,3,4,5])
    x,y,xfs_legacy = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/gitlhapdf/share/LHAPDF/legacy_up.dat", dtype='float', delimiter=' ', unpack=True)
    
    grid_x = CT10nlo.logx
    grid_y = CT10nlo.logQ2
    
    xfs_legacy = xfs_legacy[x>=1E-8]
    x = x[x>=1E-8]
    # imported data is 1D columns
    sizex = 801
    sizey = 901
    x = x[::sizey]
    y = y[:sizey]
    
    xfs_legacy = xfs_legacy.reshape((sizex,sizey)).T

    num_d1y, num_d1x = np.gradient(xfs_legacy, y, x)
    num_d2y = np.gradient(num_d1y, y, axis=0) 
    num_d2x = np.gradient(num_d1x, x, axis=1) 
    
    an_d1x = an_d1x.reshape((sizey,sizex))
    an_d1y = an_d1y.reshape((sizey,sizex))
    an_d2x = an_d2x.reshape((sizey,sizex))
    an_d2y = an_d2y.reshape((sizey,sizex))
    
    # calculate percentage errors
    ratio1 = (num_d1x/an_d1x - 1)*100
    ratio2 = (num_d1y/an_d1y - 1)*100
    ratio3 = (num_d2x/an_d2x - 1)*100
    ratio4 = (num_d2y/an_d2y - 1)*100
    
    x = np.log10(x)
    y = np.log10(y)
    
    grid_y = grid_y[grid_y>=np.min(y)]
    XX,YY = np.meshgrid(grid_x,grid_y)
    
    # Plot data
    title_size = 18
    cbar_size = 18
    axes_size = 18
    msize = 2
    alph = 0.5

    # plot 2x2 errors with actual function at top left
    fig = plt.figure(figsize=(8,8))
    axs = fig.subplots(2,2).flatten()
    
    norm,levels,ticks1 = dynamicColorbar(ratio1) # Each has seprate colorbar as they model different things
    plot1 = axs[0].pcolormesh(x,y,ratio1, norm=norm, cmap="seismic", rasterized=True)
    axs[0].scatter(XX, YY, color="black", s=msize, alpha=alph)
    #axs[0].set_title(r"$\Delta \partial_x xf(x,Q)$", fontsize=title_size)
    axs[0].set_xlabel(r"(a)", fontsize=axes_size)
    
    norm,levels,ticks2 = dynamicColorbar(ratio2) # Each has seprate colorbar as they model different things
    plot2 = axs[1].pcolormesh(x,y,ratio2, norm=norm, cmap="seismic", rasterized=True)
    axs[1].scatter(XX, YY, color="black", s=msize, alpha=alph)
    #axs[1].set_title(r"$\Delta\partial_Q xf(x,Q)$", fontsize=title_size)
    axs[1].set_xlabel(r"(b)", fontsize=axes_size)
    
    norm,levels,ticks3 = dynamicColorbar(ratio3) # Each has seprate colorbar as they model different things
    plot3 = axs[2].pcolormesh(x,y,ratio3, norm=norm, cmap="seismic", rasterized=True)
    axs[2].scatter(XX, YY, color="black", s=msize, alpha=alph)
    #axs[2].set_title(r"$\Delta \partial^2_x xf(x,Q)$", fontsize=title_size)
    axs[2].set_xlabel(r"(c)", fontsize=axes_size)
    
    norm,levels,ticks4 = dynamicColorbar(ratio4) # Each has seprate colorbar as they model different things
    plot4 = axs[3].pcolormesh(x,y,ratio4, norm=norm, cmap="seismic", rasterized=True)
    axs[3].scatter(XX, YY, color="black", s=msize, alpha=alph)
    #axs[3].set_title(r"$\Delta\partial^2_Q xf(x,Q)$", fontsize=title_size)
    axs[3].set_xlabel(r"(d)", fontsize=axes_size)
    
    cb1 = plt.colorbar(plot1, ax=axs[0], ticks=ticks1)
    cb2 = plt.colorbar(plot2, ax=axs[1], ticks=ticks2)
    cb3 = plt.colorbar(plot3, ax=axs[2], ticks=ticks3)
    cb4 = plt.colorbar(plot4, ax=axs[3], ticks=ticks4)
    
    cb1.set_label(r"$\Delta \partial_x xf(x,Q^2)$ ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb2.set_label(r"$\Delta\partial_{Q^2} xf(x,Q^2)$ ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb3.set_label(r"$\Delta \partial^2_x xf(x,Q^2)$ ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb4.set_label(r"$\Delta\partial^2_{Q^2} xf(x,Q^2)$ ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    
    fig.supxlabel(r'$\log_{10}x$', fontsize=axes_size)
    fig.supylabel(r'$\log_{10}Q^2$', x=0.01, fontsize=axes_size) 
    fig.tight_layout()
    
    fig.savefig("../outputs/LHAPDF/up/LHAPDF_deriv_err.pdf", dpi=400, format="pdf")
    
    # plot 2x2 errors new
    fig = plt.figure(figsize=(8,8))
    axs = fig.subplots(2,2).flatten()
    
    plot1 = axs[0].pcolormesh(x,y,an_d1x, rasterized=True)
    axs[0].scatter(XX, YY, color="black", s=msize, alpha=alph)
    #axs[0].set_title(r"$\Delta \partial_x xf(x,Q)$", fontsize=title_size)
    axs[0].set_xlabel(r"(a)", fontsize=axes_size)
    
    plot2 = axs[1].pcolormesh(x,y,an_d1y, rasterized=True)
    axs[1].scatter(XX, YY, color="black", s=msize, alpha=alph)
    #axs[1].set_title(r"$\Delta\partial_Q xf(x,Q)$", fontsize=title_size)
    axs[1].set_xlabel(r"(b)", fontsize=axes_size)
    
    plot3 = axs[2].pcolormesh(x,y,an_d2x, rasterized=True)
    axs[2].scatter(XX, YY, color="black", s=msize, alpha=alph)
    #axs[2].set_title(r"$\Delta \partial^2_x xf(x,Q)$", fontsize=title_size)
    axs[2].set_xlabel(r"(c)", fontsize=axes_size)
    
    plot4 = axs[3].pcolormesh(x,y,an_d2y, rasterized=True)
    axs[3].scatter(XX, YY, color="black", s=msize, alpha=alph)
    #axs[3].set_title(r"$\Delta\partial^2_Q xf(x,Q)$", fontsize=title_size)
    axs[3].set_xlabel(r"(d)", fontsize=axes_size)
    
    cb1 = plt.colorbar(plot1, ax=axs[0])
    cb2 = plt.colorbar(plot2, ax=axs[1])
    cb3 = plt.colorbar(plot3, ax=axs[2])
    cb4 = plt.colorbar(plot4, ax=axs[3])
    
    cb1.set_label(r"$\partial_x xf(x,Q)$ ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb2.set_label(r"$\partial_Q xf(x,Q)$ ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb3.set_label(r"$\partial^2_x xf(x,Q)$ ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb4.set_label(r"$\partial^2_Q xf(x,Q)$ ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    
    fig.supxlabel(r'$\log_{10}x$', fontsize=axes_size)
    fig.supylabel(r'$\log_{10}Q$', x=0.01, fontsize=axes_size) 
    fig.tight_layout()
    
    fig.savefig("../outputs/LHAPDF/up/LHAPDF_deriv_new.pdf", dpi=400, format="pdf")
    
     # plot 2x2 errors legacy
    fig = plt.figure(figsize=(8,8))
    axs = fig.subplots(2,2).flatten()
    
    plot1 = axs[0].pcolormesh(x,y,num_d1x, rasterized=True)
    axs[0].scatter(XX, YY, color="black", s=msize, alpha=alph)
    #axs[0].set_title(r"$\Delta \partial_x xf(x,Q)$", fontsize=title_size)
    axs[0].set_xlabel(r"(a)", fontsize=axes_size)
    
    plot2 = axs[1].pcolormesh(x,y,num_d1y, rasterized=True)
    axs[1].scatter(XX, YY, color="black", s=msize, alpha=alph)
    #axs[1].set_title(r"$\Delta\partial_Q xf(x,Q)$", fontsize=title_size)
    axs[1].set_xlabel(r"(b)", fontsize=axes_size)
    
    plot3 = axs[2].pcolormesh(x,y,num_d2x, rasterized=True)
    axs[2].scatter(XX, YY, color="black", s=msize, alpha=alph)
    #axs[2].set_title(r"$\Delta \partial^2_x xf(x,Q)$", fontsize=title_size)
    axs[2].set_xlabel(r"(c)", fontsize=axes_size)
    
    plot4 = axs[3].pcolormesh(x,y,num_d2y, rasterized=True)
    axs[3].scatter(XX, YY, color="black", s=msize, alpha=alph)
    #axs[3].set_title(r"$\Delta\partial^2_Q xf(x,Q)$", fontsize=title_size)
    axs[3].set_xlabel(r"(d)", fontsize=axes_size)
    
    cb1 = plt.colorbar(plot1, ax=axs[0])
    cb2 = plt.colorbar(plot2, ax=axs[1])
    cb3 = plt.colorbar(plot3, ax=axs[2])
    cb4 = plt.colorbar(plot4, ax=axs[3])
    
    cb1.set_label(r"$\partial_x xf(x,Q)$ ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb2.set_label(r"$\partial_Q xf(x,Q)$ ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb3.set_label(r"$\partial^2_x xf(x,Q)$ ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb4.set_label(r"$\partial^2_Q xf(x,Q)$ ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    
    fig.supxlabel(r'$\log_{10}x$', fontsize=axes_size)
    fig.supylabel(r'$\log_{10}Q$', x=0.01, fontsize=axes_size) 
    fig.tight_layout()
    
    fig.savefig("../outputs/LHAPDF/up/LHAPDF_deriv_legacy.pdf", dpi=400, format="pdf")
    
    
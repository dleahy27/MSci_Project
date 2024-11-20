import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats.mstats import mquantiles

plt.rcParams["text.usetex"] = True

def dynamicColorbar(z):
    # 5% and 95% quantiles
    minmax = mquantiles(z, prob=[0.05,0.95])
    # colorbar norm, levels and ticks set to each case
    # Bigger negative value
    if ( np.abs(minmax[0]) > np.abs(minmax[1]) ):
        order = int(np.floor(np.log10(np.abs(minmax[0]))))
        
        norm = TwoSlopeNorm(vmin=minmax[0], vcenter=0, vmax=-minmax[0])
        levels = np.linspace(minmax[0], -minmax[0], 100, endpoint=True)
        ticks = np.linspace(minmax[0], -minmax[0], 9, endpoint=True)
        ticks = [np.round(tick, decimals=-order+2) for tick in ticks]
        
    # Bigger positive value    
    else:
        # colorbar norm, levels and ticks set to these
        order = int(np.floor(np.log(np.abs(minmax[1]))))
        norm = TwoSlopeNorm(vmin=-minmax[1], vcenter=0, vmax=minmax[1])
        levels = np.linspace(-minmax[1], minmax[1], 100, endpoint=True)
        ticks = np.linspace(-minmax[1], minmax[1], 9, endpoint=True)
        ticks = [np.round(tick, decimals=-order+2) for tick in ticks]
        
    return norm, levels, ticks
    
if __name__ == '__main__':
    # read in data
    an_d1x, an_d1y, an_d2x, an_d2y = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/an_ds_trig.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)
    x,y,num_d1x, num_d1y, num_d2x, num_d2y = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/num_ds_trig.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)
    
    # imported data is 1D columns
    sizex = 101
    sizey = 126
    
    x = x.reshape((sizey,sizex))
    y = y.reshape((sizey,sizex))
    
    an_d1x = an_d1x.reshape((sizey,sizex))
    an_d1y = an_d1y.reshape((sizey,sizex))
    an_d2x = an_d2x.reshape((sizey,sizex))
    an_d2y = an_d2y.reshape((sizey,sizex))
    
    num_d1x = num_d1x.reshape((sizey,sizex))
    num_d1y = num_d1y.reshape((sizey,sizex))
    num_d2x = num_d2x.reshape((sizey,sizex))
    num_d2y = num_d2y.reshape((sizey,sizex))
    
    # calculate percentage errors
    ratio1 = (num_d1x/an_d1x - 1)*100
    ratio2 = (num_d1y/an_d1y - 1)*100
    ratio3 = (num_d2x/an_d2x - 1)*100
    ratio4 = (num_d2y/an_d2y - 1)*100
    
    
    # Plot data
    title_size = 14
    cbar_size = 14
    axes_size = 14

    # plot 2x2 errors with actual function at top left
    fig = plt.figure(figsize=(10,8))
    axs = fig.subplots(2,2).flatten()
    
    norm,levels,ticks1 = dynamicColorbar(ratio1) # Each has seprate colorbar as they model different things
    plot1 = axs[0].contourf(x,y,ratio1, levels=levels, norm=norm, cmap="seismic")
    axs[0].set_title(r"$\Delta \partial z / \partial x$", fontsize=title_size)
    axs[0].label_outer()
    # axs[0].set_xscale("log")
    # axs[0].set_yscale("log")
    
    norm,levels,ticks2 = dynamicColorbar(ratio2) # Each has seprate colorbar as they model different things
    plot2 = axs[1].contourf(x,y,ratio2, levels=levels, norm=norm, cmap="seismic")
    axs[1].set_title(r"$\Delta \partial z / \partial y$", fontsize=title_size)
    axs[1].label_outer()
    # axs[1].set_xscale("log")
    # axs[1].set_yscale("log")
    
    norm,levels,ticks3 = dynamicColorbar(ratio3) # Each has seprate colorbar as they model different things
    plot3 = axs[2].contourf(x,y,ratio3, levels=levels, norm=norm, cmap="seismic")
    axs[2].set_title(r"$\Delta \partial^2 z / \partial x^2$", fontsize=title_size)
    axs[2].label_outer()
    # axs[2].set_xscale("log")
    # axs[2].set_yscale("log")
    
    norm,levels,ticks4 = dynamicColorbar(ratio4) # Each has seprate colorbar as they model different things
    plot4 = axs[3].contourf(x,y,ratio4, levels=levels, norm=norm, cmap="seismic")
    axs[3].set_title(r"$\Delta \partial^2 z / \partial y^2$", fontsize=title_size)
    axs[3].label_outer()
    # axs[3].set_xscale("log")
    # axs[3].set_yscale("log")
    
    cb1 = plt.colorbar(plot1, ax=axs[0], ticks=ticks1)
    cb2 = plt.colorbar(plot2, ax=axs[1], ticks=ticks2)
    cb3 = plt.colorbar(plot3, ax=axs[2], ticks=ticks3)
    cb4 = plt.colorbar(plot4, ax=axs[3], ticks=ticks4)
    
    cb1.set_label(r"relative error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb2.set_label(r"relative error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb3.set_label(r"relative error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb4.set_label(r"relative error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    
    fig.supxlabel(r'x', fontsize=axes_size)
    fig.supylabel(r'y', x=0.01, fontsize=axes_size) 
    fig.tight_layout()
    
    fig.savefig("../outputs/trig/deriv_err.pdf", dpi=800, format="pdf")
    
    ixs = np.random.uniform(1,sizex,9)
    iys = np.random.uniform(1,sizey,9)
    
    # fig = plt.figure(figsize=(12,12))
    # axs = fig.subplots(2,2).flatten()
    # for i in np.arange(9):
    #     axs[0].plot(y[1][-3:], an_d2x[i][-3:], ls="dashed", label=r"analytic at x={:.2E}".format(x[i][1]))
    #     axs[0].plot(y[1][-3:], num_d2x[i][-3:], label=r"numerical at x={:.2E}".format(x[i][1]))
    #     axs[0].set_title(r"Top boundary")
    #     #axs[0].legend()
        
    #     axs[1].plot(y[1][:4], an_d2x[i][:4], ls="dashed", label=r"analytic at x={:.2E}".format(x[i][1]))
    #     axs[1].plot(y[1][:4], num_d2x[i][:4], label=r"numerical at x={:.2E}".format(x[i][1]))
    #     axs[1].set_title(r"Bottom boundary")
    #     # axs[1].legend()
        
        
    #     axs[2].plot(np.array([x[-3:][0][0],x[-3:][1][0],x[-3:][2][0]]), an_d2y[i][-3:], ls="dashed", label=r"analytic at y={:.2E}".format(y[1][i]))
    #     axs[2].plot(np.array([x[-3:][0][0],x[-3:][1][0],x[-3:][2][0]]), num_d2y[i][-3:], label=r"numerical at y={:.2E}".format(y[1][i]))
    #     axs[2].set_title(r"Right boundary")
    #     # axs[2].legend()
        
    #     axs[3].plot(np.array([x[:4][0][0],x[:4][1][0],x[:4][2][0]]), an_d2y[i][:4], ls="dashed", label=r"analytic at y={:.2E}".format(y[1][i]))
    #     axs[3].plot(np.array([x[:4][0][0],x[:4][1][0],x[:4][2][0]]), num_d2y[i][:4], label=r"numerical at y={:.2E}".format(y[1][i]))
    #     axs[3].set_title(r"Left boundary")
    #     # axs[3].legend()
    
    
    # fig.savefig("../outputs/trig/extrapolation_compare.pdf", format="pdf")
        
    
    
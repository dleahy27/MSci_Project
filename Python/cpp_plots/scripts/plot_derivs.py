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
    an_d1x, an_d1y, an_d2x, an_d2y = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/analytical_ds_pdf.csv", dtype='float', delimiter=',', unpack=True, skiprows=1, usecols=[2,3,4,5])
    x,y,num_d1x, num_d1y, num_d2x, num_d2y, num_d3x, num_d3y = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/num_ds_pdf.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)
    x,y,an_num_d1x, an_num_d1y, an_num_d2x, an_num_d2y, an_num_d3x, an_num_d3y = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/an_ds_pdf.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)
    # x,y,zero_d1x, zero_d1y, zero_d2x, zero_d2y, zero_d3x, zero_d3y = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/zero_ds_pdf.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)
    
    
    # imported data is 1D columns
    sizex = 53
    sizey = 53    
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
    num_d3x = num_d3x.reshape((sizey,sizex))
    num_d3y = num_d3y.reshape((sizey,sizex))
    
    an_num_d1x = an_num_d1x.reshape((sizey,sizex))
    an_num_d1y = an_num_d1y.reshape((sizey,sizex))
    an_num_d2x = an_num_d2x.reshape((sizey,sizex))
    an_num_d2y = an_num_d2y.reshape((sizey,sizex))
    an_num_d3x = an_num_d3x.reshape((sizey,sizex))
    an_num_d3y = an_num_d3y.reshape((sizey,sizex))
    
    # calculate percentage errors
    ratio1 = (num_d1x/an_d1x - 1)*100
    ratio2 = (num_d1y/an_d1y - 1)*100
    ratio3 = (num_d2x/an_d2x - 1)*100
    ratio4 = (num_d2y/an_d2y - 1)*100
    
    an_ratio1 = (an_num_d1x/an_d1x - 1)*100
    an_ratio2 = (an_num_d1y/an_d1y - 1)*100
    an_ratio3 = (an_num_d2x/an_d2x - 1)*100
    an_ratio4 = (an_num_d2y/an_d2y - 1)*100
    
    # Plot data
    title_size = 14
    cbar_size = 14
    axes_size = 14

    # plot 2x2 errors with actual function at top left
    fig = plt.figure(figsize=(8,8))
    axs = fig.subplots(2,2).flatten()
    
    norm,levels,ticks1 = dynamicColorbar(ratio1) # Each has seprate colorbar as they model different things
    plot1 = axs[0].pcolormesh(x,y,ratio1, norm=norm, cmap="seismic")
    axs[0].set_title(r"$\Delta \partial_x z$", fontsize=title_size)
    axs[0].label_outer()
    
    norm,levels,ticks2 = dynamicColorbar(ratio2) # Each has seprate colorbar as they model different things
    plot2 = axs[1].pcolormesh(x,y,ratio2, norm=norm, cmap="seismic")
    axs[1].set_title(r"$\Delta\partial_y z$", fontsize=title_size)
    axs[1].label_outer()
    
    norm,levels,ticks3 = dynamicColorbar(ratio3) # Each has seprate colorbar as they model different things
    plot3 = axs[2].pcolormesh(x,y,ratio3, norm=norm, cmap="seismic")
    axs[2].set_title(r"$\Delta \partial^2_x z$", fontsize=title_size)
    axs[2].label_outer()
    
    norm,levels,ticks4 = dynamicColorbar(ratio4) # Each has seprate colorbar as they model different things
    plot4 = axs[3].pcolormesh(x,y,ratio4, norm=norm, cmap="seismic")
    axs[3].set_title(r"$\Delta\partial^2_y z$", fontsize=title_size)
    axs[3].label_outer()
    
    cb1 = plt.colorbar(plot1, ax=axs[0], ticks=ticks1)
    cb2 = plt.colorbar(plot2, ax=axs[1], ticks=ticks2)
    cb3 = plt.colorbar(plot3, ax=axs[2], ticks=ticks3)
    cb4 = plt.colorbar(plot4, ax=axs[3], ticks=ticks4)
    
    cb2.set_label(r"Relative Error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb4.set_label(r"Relative Error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    
    fig.supxlabel(r'$\log_{10}x$', fontsize=axes_size)
    fig.supylabel(r'$\log_{10}Q$', x=0.01, fontsize=axes_size) 
    fig.tight_layout()
    
    fig.savefig("../outputs/pdflike/deriv_err.pdf", dpi=800, format="pdf")
    
    fig = plt.figure(figsize=(12,6))
    axs = fig.subplots(1,2).flatten()
    
    plot5 = axs[0].pcolormesh(x,y,num_d3x)
    axs[0].set_title(r"$\partial^3_x z$", fontsize=title_size)
    axs[0].label_outer()
    
    plot6 = axs[1].pcolormesh(x,y,num_d3y)
    axs[1].set_title(r"$\partial^3_y z$", fontsize=title_size)
    axs[1].label_outer()
    
    cb5 = plt.colorbar(plot5)
    cb6 = plt.colorbar(plot6)
    
    fig.supxlabel(r'$\log_{10}x$', fontsize=axes_size)
    fig.supylabel(r'$\log_{10}Q$', x=0.01, fontsize=axes_size) 
    fig.tight_layout()
    
    fig.savefig("../outputs/pdflike/third_derivs.pdf", dpi=800, format="pdf")
    
    # plot 2x2 errors with actual function at top left
    fig = plt.figure(figsize=(8,8))
    axs = fig.subplots(2,2).flatten()
    
    norm,levels,ticks1 = dynamicColorbar(an_ratio1) # Each has seprate colorbar as they model different things
    plot1 = axs[0].pcolormesh(x,y,an_ratio1, norm=norm, cmap="seismic")
    axs[0].set_title(r"$\Delta \partial_x z$", fontsize=title_size)
    axs[0].label_outer()
    
    norm,levels,ticks2 = dynamicColorbar(an_ratio2) # Each has seprate colorbar as they model different things
    plot2 = axs[1].pcolormesh(x,y,an_ratio2, norm=norm, cmap="seismic")
    axs[1].set_title(r"$\Delta\partial_y z$", fontsize=title_size)
    axs[1].label_outer()
    
    norm,levels,ticks3 = dynamicColorbar(an_ratio3) # Each has seprate colorbar as they model different things
    plot3 = axs[2].pcolormesh(x,y,an_ratio3, norm=norm, cmap="seismic")
    axs[2].set_title(r"$\Delta \partial^2_x z$", fontsize=title_size)
    axs[2].label_outer()
    
    norm,levels,ticks4 = dynamicColorbar(an_ratio4) # Each has seprate colorbar as they model different things
    plot4 = axs[3].pcolormesh(x,y,an_ratio4, norm=norm, cmap="seismic")
    axs[3].set_title(r"$\Delta\partial^2_y z$", fontsize=title_size)
    axs[3].label_outer()
    
    cb1 = plt.colorbar(plot1, ax=axs[0], ticks=ticks1)
    cb2 = plt.colorbar(plot2, ax=axs[1], ticks=ticks2)
    cb3 = plt.colorbar(plot3, ax=axs[2], ticks=ticks3)
    cb4 = plt.colorbar(plot4, ax=axs[3], ticks=ticks4)
    
    cb2.set_label(r"Relative Error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb4.set_label(r"Relative Error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    
    fig.supxlabel(r'$\log_{10}x$', fontsize=axes_size)
    fig.supylabel(r'$\log_{10}Q$', x=0.01, fontsize=axes_size) 
    fig.tight_layout()
    
    fig.savefig("../outputs/pdflike/an_deriv_err.pdf", dpi=800, format="pdf")    
        
    
    
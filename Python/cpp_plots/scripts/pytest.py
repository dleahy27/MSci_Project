import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats.mstats import mquantiles

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
    
if __name__ == '__main__':
    # read in data
    x_c,y_c,z_c,z_cds,z_cds_num,z_cspline = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/pdfLike_bicubic.csv", dtype='float64', delimiter=',', unpack=True, skiprows=1)
    
    # imported data is 1D columns
    sizex = 153
    sizey = 153
    x_c = x_c.reshape((sizey,sizex))
    # x_c = x_c[:,:]
    y_c = y_c.reshape((sizey,sizex))
    # y_c = y_c[:, 2:]
    z_c = z_c.reshape((sizey,sizex))
    # z_c = z_c[:, 2:]
    z_cds = z_cds.reshape((sizey,sizex))
    # z_cds = z_cds[:, 2:]
    z_cds_num = z_cds_num.reshape((sizey,sizex))
    # z_cds_num = z_cds_num[:, 2:]
    z_cspline = z_cspline.reshape((sizey,sizex))
    # z_cspline = z_cspline[:, 2:]
    
    # calculate percentage errors
    ratio = (z_c/z_cspline - 1)*100
    ratio_ds = (z_c/z_cds - 1)*100
    ratio_ds_num = (z_c/z_cds_num - 1)*100
    
    
    # Plot data
    title_size = 14
    cbar_size = 14
    axes_size = 14

    # plot the grid data for each case
    fig = plt.figure(figsize=(10,8))
    axs = fig.subplots(2,2).flatten()
    
    plot1 = axs[0].pcolormesh(x_c,y_c,z_c)
    axs[0].set_title(r"Analytic")
    axs[0].label_outer()
    
    plot2 = axs[1].pcolormesh(x_c,y_c,z_cds)
    axs[1].set_title(r"Spline with analytic second derivatives")
    axs[1].label_outer()
    
    plot3 = axs[2].pcolormesh(x_c,y_c,z_cds_num)
    axs[2].set_title(r"Spline with numerical second derivatives")
    axs[2].label_outer()
    
    plot4 = axs[3].pcolormesh(x_c,y_c,z_cspline)
    axs[3].set_title(r"Spline with zeroed second derivatives")
    axs[3].label_outer()
    
    cb1 = plt.colorbar(plot1, ax=axs[0])
    cb2 = plt.colorbar(plot2, ax=axs[1])
    cb3 = plt.colorbar(plot3, ax=axs[2])
    cb4 = plt.colorbar(plot4, ax=axs[3])
    
    cb1.set_label(r"z(x,y)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb2.set_label(r"z(x,y)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb3.set_label(r"z(x,y)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb4.set_label(r"z(x,y)", rotation=270, labelpad=20, fontsize=cbar_size)
    
    fig.supxlabel(r'$\log(x)$', fontsize=axes_size)
    fig.supylabel(r'$\log(Q)$', x=0.01, fontsize=axes_size) 
    fig.tight_layout()
    fig.savefig("../outputs/pdflike/function_plots.pdf",dpi=800, format="pdf")
    
    # plot 2x2 errors with actual function at top left
    fig = plt.figure(figsize=(10,8))
    axs = fig.subplots(2,2).flatten()
    
    plot1 = axs[0].pcolormesh(x_c,y_c,z_c)
    axs[0].set_title(r"Analytic function", fontsize=title_size)
    axs[0].label_outer()
    
    norm,levels,ticks1 = dynamicColorbar(ratio)
    plot2 = axs[1].pcolormesh(x_c,y_c,ratio, norm=norm, cmap="seismic")
    axs[1].set_title(r"Error with second derivatives = 0", fontsize=title_size)
    axs[1].label_outer()
    
    norm,levels,ticks2 = dynamicColorbar(ratio_ds_num)
    plot3 = axs[2].pcolormesh(x_c,y_c,ratio_ds, norm=norm, cmap="seismic")
    axs[2].set_title(r"Error with analytic second derivatives", fontsize=title_size)
    axs[2].label_outer()
    
    norm,levels,ticks3 = dynamicColorbar(ratio_ds_num)
    plot4 = axs[3].pcolormesh(x_c,y_c,ratio_ds_num, norm=norm, cmap="seismic")
    axs[3].set_title(r"Error with numerical second derivatives", fontsize=title_size)
    axs[3].label_outer()
    
    cb1 = plt.colorbar(plot1, ax=axs[0])
    cb2 = plt.colorbar(plot2, ax=axs[1], ticks=ticks1)
    cb3 = plt.colorbar(plot3, ax=axs[2], ticks=ticks2)
    cb4 = plt.colorbar(plot4, ax=axs[3], ticks=ticks3)
    
    cb1.set_label(r"z(x,y)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb2.set_label(r"relative error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb3.set_label(r"relative error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    cb4.set_label(r"relative error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    
    fig.supxlabel(r'$\log(x)$', fontsize=axes_size)
    fig.supylabel(r'$\log(Q)$', x=0.01, fontsize=axes_size) 
    fig.tight_layout()
    fig.savefig("../outputs/pdflike/func_and_err_milli.pdf", dpi=400, format="pdf")
    
    # plot percentage error for all
    # fig = plt.figure(figsize=(12,6))
    # axs = fig.subplots(1,3)
    
    # # Choose ratio as it will have the greatest error
    # norm,levels,ticks = dynamicColorbar(ratio)
    
    # plot1 = axs[0].pcolormesh(x_c,y_c,ratio, norm=norm, cmap="seismic")
    # axs[0].set_title(r"second derivatives = 0", fontsize=title_size)
    # axs[0].label_outer()
    
    # plot2 = axs[1].pcolormesh(x_c,y_c,ratio_ds, norm=norm, cmap="seismic")
    # axs[1].set_title(r"analytic second derivatives", fontsize=title_size)
    # axs[1].label_outer()
    
    # plot3 = axs[2].pcolormesh(x_c,y_c,ratio_ds_num, norm=norm, cmap="seismic")
    # axs[2].set_title(r"numerical second derivatives", fontsize=title_size)
    # axs[2].label_outer()
    
    # cbar3 = fig.colorbar(plot3, ticks=ticks, ax=axs[2])
    # cbar3.set_label(r"relative error ($\%$)", rotation=270, labelpad=20, fontsize=cbar_size)
    
    # fig.supxlabel(r'x', fontsize=axes_size)
    # fig.supylabel(r'y', x=0.01, fontsize=axes_size) 
    # fig.tight_layout()
    # fig.savefig("../outputs/ratios_bicubic_all.pdf", dpi=800, format="pdf")
    
    # # plot percentage error comparing numerical and analytical
    # fig = plt.figure(figsize=(12,6))
    # axs = fig.subplots(1,2)
    
    # # Choose ratio_ds_num as it will have the greatest error
    # norm,levels,ticks = dynamicColorbar(ratio_ds_num)
    
    # plot1 = axs[0].pcolormesh(x_c,y_c,ratio_ds, norm=norm, cmap="seismic")
    # axs[0].set_title(r"analytic second derivatives", fontsize=title_size)
    
    # plot2 = axs[1].pcolormesh(x_c,y_c,ratio_ds_num, norm=norm, cmap="seismic")
    # axs[1].set_title(r"numerical second derivatives", fontsize=title_size)
    
    # cbar2 = fig.colorbar(plot2, ticks=ticks, ax=axs[1])
    # cbar2.set_label(r"relative error ($\%$)", rotation=270, labelpad=20, fontsize=title_size)
    
    # fig.supxlabel(r'x', fontsize=axes_size)
    # fig.supylabel(r'y', x=0.01, fontsize=axes_size) 
    # fig.tight_layout()
    # fig.savefig("../outputs/ratios_num_vs_analytic.pdf", dpi=800, format="pdf")
    
    
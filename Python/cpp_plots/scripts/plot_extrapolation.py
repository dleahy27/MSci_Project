import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats.mstats import mquantiles

plt.rcParams["text.usetex"] = True

# NEED TO WRITE THIS
if __name__ == '__main__':
    # Plot params
    title_size = 14
    cbar_size = 14
    axes_size = 14
    legend_size = 14
    
    sizex = 51
    sizey = 51
    
    z = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/derivtest_grid.csv", dtype='float', delimiter=',', unpack=True, skiprows=1, usecols=[2])
    z_top,z_bottom = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/tb_boundary_zs.csv", dtype='float', delimiter=',', unpack=True, skiprows=1, usecols=[1,2])
    z_left,z_right = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/lr_boundary_zs.csv", dtype='float', delimiter=',', unpack=True, skiprows=1, usecols=[1,2])
    y,d2x_left, d2x_right = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/x_boundary_derivs.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)
    x,d2y_top, d2y_bottom = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/y_boundary_derivs.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)
    d4d2xd2y = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/corner_boundary_derivs.csv", dtype='float', delimiter=',', unpack=True, skiprows=1, usecols=[2])
    
    y,an_d2x_left, an_d2x_right = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/an_bound_d2x.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)
    x,an_d2y_top, an_d2y_bottom = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/an_bound_d2y.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)
    an_d4d2xd2y = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/an_bound_corners.csv", dtype='float', delimiter=',', unpack=True, skiprows=1, usecols=[2])
    
    z = z.reshape(sizex,sizey)
    
    # Boundaries
    fig = plt.figure(figsize=(8,8))
    axs = fig.subplots(2,2).flatten()

    plot = axs[0].plot(x, an_d2y_top, marker="x", linestyle='None', label=r"Analytic")
    axs[0].plot(x, d2y_top,  marker='.', linestyle='None', label=r"Extrapolation")
    axs[0].set_title(r"Top boundary",fontsize=title_size)
    axs[0].set_ylabel(r"$\partial^2_v xf(u;v)$",fontsize=axes_size)
    axs[0].label_outer()
    axs[0].legend(fontsize = legend_size, loc="upper left")
    
    plot = axs[2].plot(x, an_d2y_bottom, marker="x", linestyle='None', label=r"nalytic")
    axs[2].plot(x, d2y_bottom,  marker='.', linestyle='None', label=r"numerical extrapolation")
    axs[2].set_title(r"Bottom boundary",fontsize=title_size)
    axs[2].set_ylabel(r"$\partial^2_v xf(u;v)$",fontsize=axes_size)
    axs[2].set_xlabel(r"$u$",fontsize=axes_size)
    
    plot = axs[1].plot(y, an_d2x_right, marker="x", linestyle='None', label=r"analytic")
    axs[1].plot(y, d2x_right,  marker='.', linestyle='None', label=r"numerical extrapolation")
    axs[1].set_title(r"Right boundary",fontsize=title_size)
    axs[1].set_ylabel(r"$\partial^2_u xf(v;u)$",fontsize=axes_size)
    axs[1].tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=False)
    
    
    plot = axs[3].plot(y, an_d2x_left, marker="x", linestyle='None', label=r"analytic")
    axs[3].plot(y, d2x_left,  marker='.', linestyle='None', label=r"numerical extrapolation")
    axs[3].set_title(r"Left boundary",fontsize=title_size)
    axs[3].set_ylabel(r"$\partial^2_u xf(v;u)$",fontsize=axes_size)
    axs[3].set_xlabel(r"$v$",fontsize=axes_size)
    
    fig.tight_layout()
    fig.savefig("../outputs/pdflike/x_power_law/boundary_extrapolation.pdf", format="pdf")

    # Boundaries
    fig = plt.figure(figsize=(12,6))
    axs = fig.subplots(1,2).flatten()

    for i in range(10):
        if i==0:
            plot1, = axs[1].plot(y[-1], z_top[i],  marker='x', linestyle='None', label=r"Extrapolation")
            color1 = plot1.get_color()
            axs[1].plot(y[-4:], z.T[i,-4:],  marker='.', linestyle='None', color=color1, label=r"Analytic")
        else:
            plot1, = axs[1].plot(y[-1], z_top[i],  marker='x', linestyle='None')
            color1 = plot1.get_color()
            axs[1].plot(y[-4:], z.T[i,-4:],  marker='.', linestyle='None', color=color1)
        
        plot2, = axs[0].plot(y[0], z_bottom[i],  marker='x', linestyle='None')
        color2 = plot2.get_color()
        axs[0].plot(y[:4], z.T[i,:4],  marker='.', color=color2, linestyle='None')
     
    axs[1].legend(fontsize = title_size,loc = "upper left")
    axs[0].set_title(r"Bottom boundary",fontsize=title_size)
    
    axs[1].set_title(r"Top boundary",fontsize=title_size)
    
    fig.supxlabel(r"$v$",fontsize=title_size)
    fig.supylabel(r"$xf(v;u)$",fontsize=title_size)
    fig.tight_layout()
    fig.savefig("../outputs/pdflike/x_power_law/top_bottom_divergences.pdf", format="pdf")

    fig = plt.figure(figsize=(12,6))
    axs = fig.subplots(1,2).flatten()

    for i in range(10):
        if i==0:
            plot1, = axs[1].plot(x[-1], z_right[i],  marker='x', linestyle='None', label=r"Extrapolation")
            color1 = plot1.get_color()
            axs[1].plot(x[-4:], z[i,-4:],  marker='.', linestyle='None', color=color1, label=r"Analytic")
        else:
            plot1, = axs[1].plot(x[-1], z_right[i],  marker='x', linestyle='None')
            color1 = plot1.get_color()
            axs[1].plot(x[-4:], z[i,-4:],  marker='.', linestyle='None', color=color1)
        
        plot2, = axs[0].plot(x[0], z_left[i],  marker='x', linestyle='None')
        color2 = plot2.get_color()
        axs[0].plot(x[:4], z[i,:4],  marker='.', color=color2, linestyle='None')
     
    axs[1].legend()
    axs[0].set_title(r"Left boundary",fontsize=title_size)
           
    axs[1].set_title(r"Right boundary",fontsize=title_size)
    
    fig.supxlabel(r"$u$",fontsize=title_size)
    fig.supylabel(r"$xf(u;v)$",fontsize=title_size)
    fig.tight_layout()
    fig.savefig("../outputs/pdflike/x_power_law/left_right_divergences.pdf", format="pdf")
    
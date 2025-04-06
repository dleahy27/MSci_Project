import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline as csh
from ct10nlo import CT10nlo 

plt.rcParams["text.usetex"] = True
    
if __name__ == '__main__':
    # read in data
    an_d1x, an_d1y, an_d2x, an_d2y = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/ct10_derivs.csv", dtype='float', delimiter=',', unpack=True, skiprows=1, usecols=[2,3,4,5])
    x_bigger,y_bigger,xfs_legacy = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/gitlhapdf/share/LHAPDF/legacy_up.dat", dtype='float', delimiter=' ', unpack=True)
    xs,qs,xfs_new = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/CT10nlo.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)

    grid_xfx = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/test_down_xfx.dat", delimiter=",")
    
    grid_x = CT10nlo.logx
    grid_y = CT10nlo.logQ2
    
    y = y_bigger[x_bigger>=1E-8]
    xfs_legacy = xfs_legacy[x_bigger >= 1E-8]
    x = x_bigger[x_bigger>=1E-8]
    
    # imported data is 1D columns
    sizex = 801
    sizey = 901
    
    x = x[::sizey]
    y = y[:sizey]
    
    x_bigger = x_bigger[::sizey]
    
    xfs_legacy = xfs_legacy.reshape((sizex,sizey))
    xfs_legacy = xfs_legacy.reshape((sizex,sizey))
    xfs_new = xfs_new.reshape((sizey,sizex))
    grid_xfx = grid_xfx.reshape((grid_y.shape[0],grid_x.shape[0]))
    print(grid_xfx.shape)
    num_d1x, num_d1y = np.gradient(xfs_legacy, x, y)
    num_d2y = np.gradient(num_d1y, y, axis=1) 
    num_d2x = np.gradient(num_d1x, x, axis=0) 
    
    an_d1x = an_d1x.reshape((sizey,sizex))
    an_d1y = an_d1y.reshape((sizey,sizex))
    an_d2x = an_d2x.reshape((sizey,sizex))
    an_d2y = an_d2y.reshape((sizey,sizex))
    
    x = np.log10(x)
    y = np.log10(y)
    grid_y = grid_y[grid_y>=np.min(y)]
    x_bigger = np.log10(x_bigger)
    
    # Plot data
    title_size = 14
    cbar_size = 14
    axes_size = 14
    alph = 0.5
    
    vals_x = np.linspace(0, sizex, 10, dtype=int)
    vals_y = np.linspace(0, sizex, 10, dtype=int)
###################################################################################################################
    # Horizontal Derivatives in the bulk
    fig1 = plt.figure(figsize=(14,8))
    axs1 = fig1.subplots(2,4).flatten()
    fig2 = plt.figure(figsize=(14,8))
    axs2 = fig2.subplots(2,4).flatten()
    
    for i in range(10):
        if i == 0 or i == 9:
            continue
        plot = axs1[i-1].plot(x, an_d2y[vals_y[i]][:], marker=".",  label=r"Bicubic")
        axs1[i-1].plot(x, num_d2y.T[vals_y[i]][:],  marker='.',  label=r"Legacy")
        axs1[i-1].set_title(r"$\log_{{10}}Q={0:.2}$".format(y[vals_y[i]]),fontsize=title_size)
        [axs1[i-1].axvline(_x, color="black", alpha=alph, zorder=-1) for _x in grid_x]
        axs1[i-1].set_xscale("symlog")
        #axs1[i-1].set_ylabel(r"$\partial^2_y xf\(x,Q\)$",fontsize=axes_size)
        #axs1[i-1].set_xlabel(r"$\log_{10}x$",fontsize=axes_size)
    
        plot = axs2[i-1].plot(x, an_d1y[vals_y[i]][:], marker=".",  label=r"Bicubic")
        axs2[i-1].plot(x, num_d1y.T[vals_y[i]][:],  marker='.',  label=r"Legacy")
        axs2[i-1].set_title(r"$\log_{{10}}Q={0:.2}$".format(y[vals_y[i]]),fontsize=title_size)
        axs2[i-1].set_xscale("symlog")
        [axs2[i-1].axvline(_x, color="black", alpha=alph, zorder=-1) for _x in grid_x]
        #axs2[i-1].set_ylabel(r"$\partial_y xf\(x,Q\)$",fontsize=axes_size)
        #axs2[i-1].set_xlabel(r"$\log_{10}x$",fontsize=axes_size)
    
        if i == 1:
            axs1[i-1].legend(fontsize = title_size)
            axs2[i-1].legend(fontsize = title_size)
        
    fig1.supxlabel(r'$\log_{10}x$', fontsize=title_size+4)
    fig2.supxlabel(r'$\log_{10}x$', fontsize=title_size+4)
    
    fig1.supylabel(r'$\partial^2_Q xf(x;Q)$', fontsize=title_size+4)
    fig2.supylabel(r'$\partial_Q xf(x;Q)$', fontsize=title_size+4)
    
    fig1.suptitle(r"Horizontal Slices",fontsize=title_size+4)
    fig2.suptitle(r"Horizontal Slices",fontsize=title_size+4)
    
    fig1.tight_layout()
    fig2.tight_layout()
    
    fig1.savefig("../outputs/LHAPDF/up/d2y_horizontal_slices.pdf", format="pdf")
    fig2.savefig("../outputs/LHAPDF/up/d1y_horizontal_slices.pdf", format="pdf") 
    
    # Vertical Derivatives in the bulk
    fig1 = plt.figure(figsize=(14,8))
    axs1 = fig1.subplots(2,4).flatten()
    fig2 = plt.figure(figsize=(14,8))
    axs2 = fig2.subplots(2,4).flatten()
    
    for i in range(10):
        if i == 0 or i == 9:
            continue
        plot = axs1[i-1].plot(y, an_d2y.T[vals_x[i]][:], marker=".",  label=r"Bicubic")
        axs1[i-1].plot(y, num_d2y[vals_x[i]][:],  marker='.',  label=r"Legacy")
        axs1[i-1].set_title(r"$\log_{{10}}x={0:.2}$".format(x[vals_x[i]]),fontsize=title_size)
        axs1[i-1].set_xscale("log")
        [axs1[i-1].axvline(_y, color="black", alpha=alph, zorder=-1) for _y in grid_y]
        #axs1[i-1].set_ylabel(r"$\partial^2_y xf\(x,Q\)$",fontsize=axes_size)
        #axs1[i-1].set_xlabel(r"$\log_{10}x$",fontsize=axes_size)
    
        plot = axs2[i-1].plot(y, an_d1y.T[vals_x[i]][:], marker=".",  label=r"Bicubic")
        axs2[i-1].plot(y, num_d1y[vals_x[i]][:],  marker='.',  label=r"Legacy")
        axs2[i-1].set_title(r"$\log_{{10}}x={0:.2}$".format(x[vals_x[i]]),fontsize=title_size)
        axs2[i-1].set_xscale("log")
        [axs2[i-1].axvline(_y, color="black", alpha=alph, zorder=-1) for _y in grid_y]
        #axs2[i-1].set_ylabel(r"$\partial_y xf\(x,Q\)$",fontsize=axes_size)
        #axs2[i-1].set_xlabel(r"$\log_{10}x$",fontsize=axes_size)
    
        if i == 1:
            axs1[i-1].legend(fontsize = title_size)
            axs2[i-1].legend(fontsize = title_size)
        if i == 1:
            fig = plt.figure(figsize=(6,6))
            ax = fig.subplots(1,1)
            plot = ax.plot(y, an_d2y.T[vals_x[i]][:], marker=".",  label=r"Bicubic")
            ax.plot(y, num_d2y[vals_x[i]][:],  marker='.',  label=r"Legacy")
            ax.set_title(r"$\log_{{10}}x={0:.2}$".format(x[vals_x[i]]),fontsize=title_size)
            ax.set_xscale("log")
            [ax.axvline(_y, color="black", alpha=alph, zorder=-1) for _y in grid_y]
            ax.set_xlabel(r'$\log_{10}Q^2$', fontsize=title_size+4)
            ax.set_ylabel(r'$\partial^2_{Q^2} xf(Q^2;x)$', fontsize=title_size+4)
            ax.set_xlim(0.9,2.5)
            fig.tight_layout()
            fig.savefig("../outputs/LHAPDF/up/deriv_slice.pdf", format="pdf")
            
        
    fig1.supxlabel(r'$\log_{10}Q$', fontsize=title_size+4)
    fig2.supxlabel(r'$\log_{10}Q$', fontsize=title_size+4)
    
    fig1.supylabel(r'$\partial^2_{Q^2} xf(Q;x)$', fontsize=title_size+4)
    fig2.supylabel(r'$\partial_Q^2 xf(Q;x)$', fontsize=title_size+4)
    
    fig1.suptitle(r"Vertical Slices",fontsize=title_size+4)
    fig2.suptitle(r"Vertical Slices",fontsize=title_size+4)
    
    fig1.tight_layout()
    fig2.tight_layout()
    
    fig1.savefig("../outputs/LHAPDF/up/d2y_vertical_slices.pdf", format="pdf")
    fig2.savefig("../outputs/LHAPDF/up/d1y_vertical_slices.pdf", format="pdf") 
    
###################################################################################################################################
    # Vertical Derivatives in the bulk
    fig1 = plt.figure(figsize=(14,8))
    axs1 = fig1.subplots(2,4).flatten()
    fig2 = plt.figure(figsize=(14,8))
    axs2 = fig2.subplots(2,4).flatten()

    fig_cubic = plt.figure(figsize=(14,8))
    axs_cubic = fig_cubic.subplots(2,4).flatten()
    pp = np.linspace(1,10,1000)
    for i in range(10):
        if i == 0 or i == 9:
            continue
        
        mask = (y>=1) & (y<=10)
        print(grid_y.shape)
        #print(grid_xfx.T[vals_x[i]][:].shape)
        #spline = csh(grid_y,grid_xfx.T[vals_x[i]][:])
        
        #val = spline(pp, 2)
        #axs_cubic[i-1].plot(pp, val, marker=".",  label=r"Cubic")
        #[axs_cubic[i-1].axvline(_y, color="black", alpha=alph, zorder=-1) for _y in grid_y]
        plot = axs1[i-1].plot(y[mask], an_d2x.T[vals_x[i]][mask], marker=".",  label=r"Bicubic")
        axs1[i-1].plot(y[mask], num_d2x[vals_x[i]][mask],  marker='.',  label=r"Legacy")
        axs1[i-1].set_title(r"$\log_{{10}}x={0:.2}$".format(x[vals_x[i]]),fontsize=title_size)
        #axs1[i-1].set_xscale("log")
        #axs1[i-1].set_xlim(1.5,2)
        [axs1[i-1].axvline(_y, color="black", alpha=alph, zorder=-1) for _y in grid_y]
        # axs1[i-1].set_ylabel(r"$\partial^2_x xf\(x,Q\)$",fontsize=axes_size)
        # axs1[i-1].set_xlabel(r"$\log_{10}Q$",fontsize=axes_size)
        
        plot = axs2[i-1].plot(y[:], an_d1x.T[vals_x[i]][:], marker=".",  label=r"Bicubic")
        axs2[i-1].plot(y[:], num_d1x[vals_x[i]][:],  marker='.',  label=r"Legacy")
        axs2[i-1].set_title(r"$\log_{{10}}x={0:.2}$".format(x[vals_x[i]]),fontsize=title_size)
        axs2[i-1].set_xscale("log")
        [axs2[i-1].axvline(_y, color="black", alpha=alph, zorder=-1) for _y in grid_y]
        # axs2[i-1].set_ylabel(r"$\partial_x xf\(x,Q\)$",fontsize=axes_size)
        # axs2[i-1].set_xlabel(r"$\log_{10}Q$",fontsize=axes_size)
    
        if i == 1:
            axs1[i-1].legend(fontsize = title_size)
            axs2[i-1].legend(fontsize = title_size)
###############################################################################################        
        if i == 5:
            mask = (y>=1.5) & (y<=2)
            fig = plt.figure(figsize=(4,4))
            ax = fig.subplots(1,1)
            ax.plot(y[mask], an_d2x.T[vals_x[i]][mask], marker=".",  label=r"Bicubic")
            ax.plot(y[mask], num_d2x[vals_x[i]][mask], marker=".",  label=r"Legacy")
            #ax.set_title(r"$\log_{{10}}x={0:.2}$".format(x[vals_x[i]]),fontsize=18)
            [ax.axvline(_y, color="black", alpha=alph, zorder=-1) for _y in grid_y]
            ax.set_xlim(1.5,2)
            ax.set_xlabel(r'$\log_{10}Q^2$', fontsize=14)
            ax.set_ylabel(r'$\partial^2_{xx} xf(Q^2;x)$', fontsize=16)
            fig.tight_layout()
            fig.savefig("../outputs/LHAPDF/up/deriv_slice_vert.pdf", format="pdf")
        
    fig1.supxlabel(r'$\log_{10}Q$', fontsize=title_size+4)
    fig2.supxlabel(r'$\log_{10}Q$', fontsize=title_size+4)
    
    fig1.supylabel(r'$\partial^2_x xf(Q;x)$', fontsize=title_size+4)
    fig2.supylabel(r'$\partial_x xf(Q;x)$', fontsize=title_size+4)
    
    fig1.suptitle(r"Vertical Slices",fontsize=title_size+4)
    fig2.suptitle(r"Vertical Slices",fontsize=title_size+4)
    
    fig1.tight_layout()
    fig2.tight_layout()
    
    fig_cubic.savefig("../outputs/LHAPDF/up/CUBIC.pdf", format="pdf")
    fig1.savefig("../outputs/LHAPDF/up/d2x_vertical_slices.pdf", format="pdf")
    fig2.savefig("../outputs/LHAPDF/up/d1x_vertical_slices.pdf", format="pdf")
    
    # Horizontal Derivatives in the bulk
    fig1 = plt.figure(figsize=(14,8))
    axs1 = fig1.subplots(2,4).flatten()
    fig2 = plt.figure(figsize=(14,8))
    axs2 = fig2.subplots(2,4).flatten()
    
    for i in range(10):
        if i == 0 or i == 9:
            continue
        plot = axs1[i-1].plot(x, an_d2x[vals_y[i]][:], marker=".",  label=r"Bicubic")
        axs1[i-1].plot(x, num_d2x.T[vals_y[i]][:],  marker='.',  label=r"Legacy")
        axs1[i-1].set_title(r"$\log_{{10}}Q={0:.2}$".format(y[vals_y[i]]),fontsize=title_size)
        #axs1[i-1].set_xscale("symlog")
        [axs1[i-1].axvline(_x, color="black", alpha=alph, zorder=-1) for _x in grid_x]
        #axs1[i-1].set_xlim(-8,-7.5)
        #print(an_d2x[vals_y[i]][:][x>=-2][np.logical_not(np.isnan(an_d2x[vals_y[i]][:][x>=-2]))])
        #axs1[i-1].set_ylim(0,np.max(an_d2x[vals_y[i]][:][x>=-2][np.logical_not(np.isnan(an_d2x[vals_y[i]][:][x>=-2]))]))
        #axs1[i-1].set_ylabel(r"$\partial^2_y xf\(x,Q\)$",fontsize=axes_size)
        #axs1[i-1].set_xlabel(r"$\log_{10}x$",fontsize=axes_size)
    
        plot = axs2[i-1].plot(x, an_d1x[vals_y[i]][:], marker=".",  label=r"Bicubic")
        axs2[i-1].plot(x, num_d1x.T[vals_y[i]][:],  marker='.',  label=r"Legacy")
        axs2[i-1].set_title(r"$\log_{{10}}Q={0:.2}$".format(y[vals_y[i]]),fontsize=title_size)
        axs2[i-1].set_xscale("symlog")
        [axs2[i-1].axvline(_x, color="black", alpha=alph, zorder=-1) for _x in grid_x]
        #axs2[i-1].set_ylabel(r"$\partial_y xf\(x,Q\)$",fontsize=axes_size)
        #axs2[i-1].set_xlabel(r"$\log_{10}x$",fontsize=axes_size)
    
        if i == 1:
            axs1[i-1].legend(fontsize = title_size)
            axs2[i-1].legend(fontsize = title_size)
        
    fig1.supxlabel(r'$\log_{10}x$', fontsize=title_size+4)
    fig2.supxlabel(r'$\log_{10}x$', fontsize=title_size+4)
    
    fig1.supylabel(r'$\partial^2_x xf(x;Q)$', fontsize=title_size+4)
    fig2.supylabel(r'$\partial_x xf(x;Q)$', fontsize=title_size+4)
    
    fig1.suptitle(r"Horizontal Slices",fontsize=title_size+4)
    fig2.suptitle(r"Horizontal Slices",fontsize=title_size+4)
    
    fig1.tight_layout()
    fig2.tight_layout()
    
    fig1.savefig("../outputs/LHAPDF/up/d2x_horizontal_slices.pdf", format="pdf")
    fig2.savefig("../outputs/LHAPDF/up/d1x_horizontal_slices.pdf", format="pdf")
    
###################################################################################################################################
    vals_x = np.linspace(0, sizex, 10, dtype=int)
    # Vertical Derivatives in the bulk
    fig1 = plt.figure(figsize=(14,8))
    axs1 = fig1.subplots(2,4).flatten()
    fig2 = plt.figure(figsize=(14,8))
    axs2 = fig2.subplots(2,4).flatten()

    for i in range(10):
        if i == 0 or i == 9:
            continue
        
        plot = axs1[i-1].plot(y[:], xfs_new.T[vals_x[i]][:], marker=".",  label=r"Bicubic")
        axs1[i-1].plot(y[:], xfs_legacy[vals_x[i]][:],  marker='.',  label=r"Legacy")
        axs1[i-1].set_title(r"$\log_{{10}}x={0:.2}$".format(x_bigger[vals_x[i]]),fontsize=title_size)
        [axs1[i-1].axvline(_y, color="black", alpha=alph, zorder=-1) for _y in grid_y]
        # axs1[i-1].set_ylabel(r"$\partial^2_x xf\(x,Q\)$",fontsize=axes_size)
        # axs1[i-1].set_xlabel(r"$\log_{10}Q$",fontsize=axes_size)
        
        plot = axs2[i-1].plot(x[:], xfs_new[vals_y[i]][:], marker=".",  label=r"Bicubic")
        axs2[i-1].plot(x[:], xfs_legacy.T[vals_y[i]][:],  marker='.',  label=r"Legacy")
        axs2[i-1].set_title(r"$\log_{{10}}Q={0:.2}$".format(y[vals_y[i]]),fontsize=title_size)
        [axs2[i-1].axvline(_x, color="black", alpha=alph, zorder=-1) for _x in grid_x]
        # axs2[i-1].set_ylabel(r"$\partial_x xf\(x,Q\)$",fontsize=axes_size)
        # axs2[i-1].set_xlabel(r"$\log_{10}Q$",fontsize=axes_size)
    
        if i == 1:
            axs1[i-1].legend(fontsize = title_size)
            axs2[i-1].legend(fontsize = title_size)
        
    fig1.supxlabel(r'$\log_{10}x$', fontsize=title_size+4)
    fig2.supxlabel(r'$\log_{10}Q$', fontsize=title_size+4)
    
    fig1.supylabel(r'$\partial^2_x xf(x;Q)$', fontsize=title_size+4)
    fig2.supylabel(r'$\partial_x xf(Q;x)$', fontsize=title_size+4)
    
    fig1.suptitle(r"Horizontal Slices",fontsize=title_size+4)
    fig2.suptitle(r"Vertical Slices",fontsize=title_size+4)
    
    fig1.tight_layout()
    fig2.tight_layout()
    
    fig1.savefig("../outputs/LHAPDF/up/xfx_horizontal_slices.pdf", format="pdf")
    fig2.savefig("../outputs/LHAPDF/up/xfx_vertical_slices.pdf", format="pdf")
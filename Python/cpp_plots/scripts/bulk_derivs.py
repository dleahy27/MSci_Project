import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["text.usetex"] = True

if __name__ == '__main__':
    # read in data
    an_d1x, an_d1y, an_d2x, an_d2y = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/analytical_ds_pdf.csv", dtype='float', delimiter=',', unpack=True, skiprows=1, usecols=[2,3,4,5])
    xs,ys,num_d1x, num_d1y, num_d2x, num_d2y, num_d3x, num_d3y = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/num_ds_pdf.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)
    
    # imported data is 1D columns
    sizex = 153
    sizey = 153
    
    xs = xs.reshape((sizey,sizex))
    ys = ys.reshape((sizey,sizex))       
    
    an_d1x = an_d1x.reshape((sizey,sizex))
    an_d1y = an_d1y.reshape((sizey,sizex))
    an_d2x = an_d2x.reshape((sizey,sizex))
    an_d2y = an_d2y.reshape((sizey,sizex))
    
    num_d1x = num_d1x.reshape((sizey,sizex))
    num_d1y = num_d1y.reshape((sizey,sizex))
    num_d2x = num_d2x.reshape((sizey,sizex))
    num_d2y = num_d2y.reshape((sizey,sizex))
    
    x = xs[0][:]
    y = ys.T[0][:]
    
    # Plot data
    title_size = 14
    cbar_size = 14
    axes_size = 14
    
    # Horizontal Derivatives in the bulk
    fig1 = plt.figure(figsize=(12,8))
    axs1 = fig1.subplots(2,4).flatten()
    fig2 = plt.figure(figsize=(12,8))
    axs2 = fig2.subplots(2,4).flatten()
    
    vals = np.linspace(0, sizex, 10, dtype=int)
    for i in range(10):
        if i == 0 or i == 9:
            continue
        plot = axs1[i-1].plot(y, an_d2y.T[vals[i]][:], marker="x", linestyle='None', label=r"analytic")
        axs1[i-1].plot(y, num_d2y.T[vals[i]][:],  marker='.', linestyle='None', label=r"numerical extrapolation")
        axs1[i-1].set_title(r"$\log_{{10}}Q={0:.2}$".format(y[vals[i]]),fontsize=title_size)
        #axs1[i-1].set_ylabel(r"$\partial^2_y z(x,y)$",fontsize=axes_size)
        #axs1[i-1].set_xlabel(r"$\log_{10}x$",fontsize=axes_size)
        
        plot = axs2[i-1].plot(y, an_d1y.T[vals[i]][:], marker="x", linestyle='None', label=r"analytic")
        axs2[i-1].plot(y, num_d1y.T[vals[i]][:],  marker='.', linestyle='None', label=r"numerical extrapolation")
        axs2[i-1].set_title(r"$\log_{{10}}Q={0:.2}$".format(y[vals[i]]),fontsize=title_size)
        #axs2[i-1].set_ylabel(r"$\partial_y z(x,y)$",fontsize=axes_size)
        #axs2[i-1].set_xlabel(r"$\log_{10}x$",fontsize=axes_size)
    
    fig1.supxlabel(r'$\log_{10}x$', fontsize=title_size+4)
    fig2.supxlabel(r'$\log_{10}x$', fontsize=title_size+4)
    
    fig1.supylabel(r'$\partial^2_y z(x,y)$', fontsize=title_size+4)
    fig2.supylabel(r'$\partial_y z(x,y)$', fontsize=title_size+4)
    
    fig1.suptitle(r"Horizontal Slices",fontsize=title_size+4)
    fig2.suptitle(r"Horizontal Slices",fontsize=title_size+4)
    
    fig1.savefig("../outputs/pdflike/bulk_2derivs_horizontal.pdf", format="pdf")
    fig2.savefig("../outputs/pdflike/bulk_1derivs_horizontal.pdf", format="pdf") 
    
    # Vertical Derivatives in the bulk
    fig1 = plt.figure(figsize=(12,8))
    axs1 = fig1.subplots(2,4).flatten()
    fig2 = plt.figure(figsize=(12,8))
    axs2 = fig2.subplots(2,4).flatten()
    
    vals = np.linspace(0, sizey, 10, dtype=int)
    for i in range(10):
        if i == 0 or i == 9:
            continue
        
        plot = axs1[i-1].plot(y[:], an_d2x[vals[i]][:], marker="x", linestyle='None', label=r"analytic")
        axs1[i-1].plot(y[:], num_d2x[vals[i]][:],  marker='.', linestyle='None', label=r"numerical extrapolation")
        axs1[i-1].set_title(r"$\log_{{10}}x={0:.2}$".format(x[vals[i]]),fontsize=title_size)
        # axs1[i-1].set_ylabel(r"$\partial^2_x z(x,y)$",fontsize=axes_size)
        # axs1[i-1].set_xlabel(r"$\log_{10}Q$",fontsize=axes_size)
        
        plot = axs2[i-1].plot(y[:], an_d1x[vals[i]][:], marker="x", linestyle='None', label=r"analytic")
        axs2[i-1].plot(y[:], num_d1x[vals[i]][:],  marker='.', linestyle='None', label=r"numerical extrapolation")
        axs2[i-1].set_title(r"$\log_{{10}}x={0:.2}$".format(x[vals[i]]),fontsize=title_size)
        # axs2[i-1].set_ylabel(r"$\partial_x z(x,y)$",fontsize=axes_size)
        # axs2[i-1].set_xlabel(r"$\log_{10}Q$",fontsize=axes_size)
    
    fig1.supxlabel(r'$\log_{10}Q$', fontsize=title_size+4)
    fig2.supxlabel(r'$\log_{10}Q$', fontsize=title_size+4)
    
    fig1.supylabel(r'$\partial^2_x z(x,y)$', fontsize=title_size+4)
    fig2.supylabel(r'$\partial_x z(x,y)$', fontsize=title_size+4)
    
    fig1.suptitle(r"Vertical Slices",fontsize=title_size+4)
    fig2.suptitle(r"Vertical Slices",fontsize=title_size+4)
    
    fig1.savefig("../outputs/pdflike/bulk_2derivs_vertical.pdf", format="pdf")
    fig2.savefig("../outputs/pdflike/bulk_1derivs_vertical.pdf", format="pdf")
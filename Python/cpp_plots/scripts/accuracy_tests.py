import numpy as np
import matplotlib.pyplot as plt


plt.rcParams["text.usetex"] = True

if __name__=='__main__':
    # read in data
    ns,ms,ts,errors = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/accuracy_test.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)
    
    ts = ts.reshape((101,101))
    errors = errors.reshape(((101,101)))
    
    # 2D plots
    fig = plt.figure()
    axs = fig.subplots(1,2)

    plot1 = axs[0].contourf(ms,ns,ts)
    axs[0].set_title(r"Processor Time")
    axs[0].label_outer()
    
    plot2 = axs[1].contourf(ms,ns,errors)
    axs[1].set_title(r"Total Error")
    axs[1].label_outer()
    
    
    plt.colorbar(plot1, ax=axs[0], label=r"$\sum_{m,n}\Delta z$")
    plt.colorbar(plot2, ax=axs[1], label=r"t(s)")

    fig.supxlabel(r'm')
    fig.supylabel(r'n')
    
    fig.savefig(r"../outputs/accuracy_plots2D")
    
    # 1D plots
    fig = plt.figure()
    axs = fig.subplots(2,2).flatten

    plot1 = axs[0].plot(ms[50,:],ts[:,50])
    axs[0].set_title(r"Processor Time")
    axs[0].label_outer()
    
    plot2 = axs[1].plot(ms[50,:],errors[:,50])
    axs[1].set_title(r"Total Error")
    axs[1].label_outer()
    
    plot3 = axs[2].contourf(ns[50,:],ts[50,:])
    axs[2].set_title(r"Processor Time")
    axs[2].label_outer()
    
    plot4 = axs[3].contourf(ns[50,:],errors[50,:])
    axs[3].set_title(r"Total Error")
    axs[3].label_outer()
    
    fig.savefig(r"../outputs/accuracy_plots1D")        
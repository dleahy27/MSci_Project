import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import lognorm


plt.rcParams["text.usetex"] = True

if __name__=='__main__':
    
    ts_init, ts_eval = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/outputs/time_test_1000.csv", dtype='float', delimiter=',', unpack=True, skiprows=1)
    
    fig = plt.figure()
    axs = fig.subplots(1,2)
    
    # Fit a normal distribution to the data:
    # shape, loc, scale = lognorm.fit(ts)
    
    n, bins, patches = axs[0].hist(ts_init,25, density=True)
    mean = np.mean(ts_init)
    axs[0].axvline(mean, ls="dashed", color="black", label=r"mean time={0:.4f}".format(mean))
    axs[0].set_title(r"Bicubic Initialisation")
    axs[0].label_outer()
    axs[0].legend()
    
    
    n, bins, patches = axs[1].hist(ts_eval,25, density=True)
    mean = np.mean(ts_eval)
    axs[1].axvline(mean, ls="dashed", color="black", label=r"mean time={0:.4f}".format(mean))
    axs[1].set_title(r"10,000 Grid Evaluations")
    axs[1].label_outer()
    axs[1].legend()

    # x = np.linspace(np.min(ts), np.max(ts), 1000)
    # p = lognorm.pdf(x,shape, loc, scale)
    
    # plt.plot(x,p,label=r'Log-normal Fit')
    # plt.xlim(0.5,0.7)

    fig.supxlabel(r'Time(s)')
    fig.supylabel(r'Density')
    
    fig.savefig(r"../outputs/time_histo_1000")    
    
    
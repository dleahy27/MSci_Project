import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats.mstats import mquantiles

plt.rcParams["text.usetex"] = True
    
if __name__ == '__main__':
    n_x1, n_q1, new = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/gitlhapdf/share/LHAPDF/new_times.dat", dtype='float', delimiter=' ', unpack=True)
    n_x2, n_q2, legacy = np.loadtxt("/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/gitlhapdf/share/LHAPDF/legacy_times.dat", dtype='float', delimiter=' ', unpack=True)
    
    fsize = 12
    
    total1 = n_x1*n_q1
    total2 = n_x2*n_q2
    
    fig = plt.figure(figsize=(4,4))
    ax = fig.subplots(1,1)
    
    ax.plot(total1, new, label=r"Bicubic")
    ax.plot(total2, legacy, label=r"Legacy")
    ax.set_ylabel(r"Time (s)", fontsize = fsize)
    ax.set_xlabel(r"Number of evaluations", fontsize = fsize)
    ax.legend(fontsize = fsize)
    
    fig.tight_layout()
    fig.savefig("../outputs/LHAPDF/times.pdf", format="pdf")
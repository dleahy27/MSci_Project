import numpy as np 
import matplotlib.pyplot as plt
import os
import sys
import glob
import re

plt.rcParams["text.usetex"] = True

def PdgNameConverter(pdg_name):
    match pdg_name:
        case -1:
            return r"d$\mathbf{^+}$"
        case -2:
            return r"u$\mathbf{^-}$"
        case -3:
            return r"s$\mathbf{^+}$"
        case -4:
            return r"c$\mathbf{^-}$"
        case -5:
            return r"b$\mathbf{^+}$"
        case 1:
            return r"d$\mathbf{^-}$"
        case 2:
            return r"u$\mathbf{^+}$"
        case 3:
            return r"s$\mathbf{^-}$"
        case 4:
            return r"c$\mathbf{^+}$"
        case 5:
            return r"b$\mathbf{^-}$"
        case 21:
            return r"g"
        case _:
            print(f"PDG number {pdg_name} is not a quark/anti-quark or a gluon.")
            return
    
    return

if __name__ == '__main__':
    input_path = "/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/C++/test"
    filename = "CT10nlo_0_*.dat"
    filename_glob = os.path.join(input_path,filename)
    files = np.array(glob.glob(filename_glob))
    
    Q2_plot = np.logspace(1, 8, num=8, base=10 )
    
    fig = plt.figure(figsize=(16,8))
    axs = fig.subplots(2,4)
    for file in files:
        x,Q2,fx = np.loadtxt(file, dtype='float', usecols=(0,1,2), unpack=True)
        pname = PdgNameConverter([int(d) for d in re.findall('-?\d+', file)][-1])
        
        j = 0
        for ax in axs.flat:
            ax.plot(np.log(x[Q2 == Q2_plot[j]]), fx[Q2 == Q2_plot[j]]*x[Q2 == Q2_plot[j]], label=pname)
            ax.set_title(r"Q$^2$={:.1E}".format(Q2_plot[j]))
            ax.set_ylim(0,0.18)
            ax.set_xlim(-6,0)
            if j == 0:
                ax.legend(fontsize=8)
            j+=1
    
    fig.supxlabel(r'$\log$(x)')
    fig.supylabel(r'xf(x)') 
    fig.tight_layout()
    fig.savefig("plot.pdf", format="pdf")
        
        
        
    
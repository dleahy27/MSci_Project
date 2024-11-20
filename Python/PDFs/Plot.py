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
            return r"$\bar{d}$"
        case -2:
            return r"$\bar{u}$"
        case -3:
            return r"$\bar{s}$"
        case -4:
            return r"$\bar{c}$"
        case -5:
            return r"$\bar{b}$"
        case 1:
            return r"$d$"
        case 2:
            return r"$u$"
        case 3:
            return r"$s$"
        case 4:
            return r"$c$"
        case 5:
            return r"$b$"
        case 21:
            return r"$g$"
        case _:
            print(f"PDG number {pdg_name} is not a quark/anti-quark or a gluon.")
            return
    
    return

if __name__ == '__main__':
    input_path = "/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/MSci_Project/C++/test"
    filename = "CT10nlo_0_*.dat"
    filename_glob = os.path.join(input_path,filename)
    files = np.array(glob.glob(filename_glob))
    
    title_sze = 20
    lbl_sze = 20
    legend_sze = 13
    
    Q2_plot = np.logspace(1, 8, num=8, base=10 )
    
    fig = plt.figure(figsize=(16,8))
    axs = fig.subplots(2,4)
    for file in files:
        x,Q2,fx = np.loadtxt(file, dtype='float', usecols=(0,1,2), unpack=True)
        pname = PdgNameConverter([int(d) for d in re.findall('-?\d+', file)][-1])
        
        j = 0
        for ax in axs.flat:
            ax.plot(np.log(x[Q2 == Q2_plot[j]]), fx[Q2 == Q2_plot[j]], label=pname)
            ax.set_title(r"$Q^2={:.1E}\left(\mathrm{{eV}}\right)^2$".format(Q2_plot[j]), size=title_sze)
            ax.set_ylim(0,1)
            ax.set_xlim(-8,0)
            ax.label_outer()
            if j == 0:
                ax.legend(fontsize=legend_sze, loc='upper left')
            j+=1
    
    fig.supxlabel(r'$\log$(x)', fontsize=lbl_sze)
    fig.supylabel(r'$xf(x;Q^2)$', x=0.01, fontsize=lbl_sze) 
    fig.tight_layout()
    fig.savefig("intro_pdfs.pdf", format="pdf")
        
        
        
    
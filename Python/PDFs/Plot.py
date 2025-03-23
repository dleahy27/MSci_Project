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
    legend_sze = 20
    
    Q2_plot = np.logspace(1, 8, num=8, base=10 )
    
    fig = plt.figure(figsize=(16,8))
    axs = fig.subplots(2,4)
    for file in files:
        x,Q2,fx = np.loadtxt(file, dtype='float', usecols=(0,1,2), unpack=True)
        pname = PdgNameConverter([int(d) for d in re.findall('-?\d+', file)][-1])
        
        j = 0
        for ax in axs.flat:
            if j ==0:
                ax.plot(np.log(x[Q2 == Q2_plot[j]]), fx[Q2 == Q2_plot[j]], label=pname)
                ax.set_title(r"$Q^2=10\left(\mathrm{{eV}}\right)^2$", size=title_sze+2)
                ax.set_ylim(0,1)
                ax.set_xlim(-8,0)
                ax.label_outer()
                j += 1
            else:
                ax.plot(np.log(x[Q2 == Q2_plot[j]]), fx[Q2 == Q2_plot[j]])
                ax.set_title(r"$Q^2=10^{:}\left(\mathrm{{eV}}\right)^2$".format(int(np.log10(Q2_plot[j]))), size=title_sze+2)
                ax.set_ylim(0,1)
                ax.set_xlim(-8,0)
                ax.label_outer()
                j += 1
            
    fig.legend(fontsize=legend_sze, loc='outside right center')
    fig.supxlabel(r'$\log_{10}(x)$', fontsize=lbl_sze+6)
    fig.supylabel(r'$xf(x;Q^2)$', x=0.08, fontsize=lbl_sze+6) 
    #fig.tight_layout()
    fig.savefig("intro_pdfs.pdf", format="pdf", bbox_inches='tight', pad_inches=0)
        
        
        
    
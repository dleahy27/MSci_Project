import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline,lagrange, CloughTocher2DInterpolator,griddata, LinearNDInterpolator
import lhapdf
import os
import sys
import glob
import re
import matplotlib

matplotlib.use("tkagg")
plt.rcParams["text.usetex"] = True

pi = np.pi

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Testing out using 1D spline (scipy) over x
# Local polynomial fit over Q^2
# Inputs: pdf name and set
# Outputs: Interpolated map
# Current Plan: Test for solely the up quark => extend to all quarks
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

def PdgNameConverter(pdg_name):
    match pdg_name:
        case -1:
            return r"\bar{d}"
        case -2:
            return r"\bar{u}"
        case -3:
            return r"\bar{s}"
        case -4:
            return r"\bar{c}"
        case -5:
            return r"\bar{b}"
        case 1:
            return r"d"
        case 2:
            return r"u"
        case 3:
            return r"s"
        case 4:
            return r"c"
        case 5:
            return r"b"
        case 21:
            return r"g"
        case _:
            print(f"PDG number {pdg_name} is not a quark/anti-quark or a gluon.")
            return
    
    return

if __name__ == '__main__':
    input_path = "/mnt/c/Users/dillo/Desktop/work/Uni/Year_5/Project/Code/C++/test"
    filename = "CT10nlo_0_*.dat"
    file = os.path.join(input_path,"CT10nlo_0_1.dat")
    # filename_glob = os.path.join(input_path,filename)
    # files = np.array(glob.glob(filename_glob))
    
    x,q2,fx = np.loadtxt(file, dtype='float', usecols=(0,1,2), unpack=True)
    
    # x_spline = CubicSpline(x[q2==q2[0]],fx[q2==q2[0]])
    # q2_polynomial = lagrange(q2[x==x[0]], fx[x==x[0]], )
    
    spline = LinearNDInterpolator(list(zip(x,q2)),fx)
    

    test_x = np.linspace(min(x),max(x),10000)
    test_fx = spline(test_x,np.ones(test_x.shape[0])*q2[1])    
    # def test(a,x):
    #     return a*np.sin(x)*np.cos(x)
    
    # def test2D(a,x,y):
    #     return a*np.sin(x)*np.cos(x)*np.sin(y)*np.cos(y)
    
    
    figure = plt.figure()
    plt.plot(np.log10(test_x),test_fx, label="Spline")
    plt.plot(np.log10(x[q2 == q2[1]]),fx[q2 == q2[1]], label="Raw Data")
    plt.title(f"q2={q2[1]}")
    plt.legend()
    plt.show()
    
    # figure = plt.figure()
    # plt.plot(np.log10(test_x),q2_polynomial(test_x), label="Local Polynomial")
    # plt.plot(np.log10(q2[q2 == q2[0]]),fx[q2 == q2[0]], label="Raw Data")
    # plt.title(f"q2={q2[0]}")
    # plt.legend()
    plt.show()
    
    
    
    
    
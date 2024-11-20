import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import matplotlib

matplotlib.use("tkagg")
plt.rcParams["text.usetex"] = True

pi = np.pi

def testFunc(x):
    
    return 5*np.exp(-((x)**2)/3)*(np.sin(x-3))**2

if __name__ == '__main__':
    lbl_sze = 18
    title_sze = 18
    test_x1= np.linspace(-5,5,4)
    test_x2 = np.linspace(-5,5,15)
    print(test_x1)


    test_y1 = testFunc(test_x1)
    test_y2 = testFunc(test_x2)
    print(test_y1)

    spline1 = CubicSpline(test_x1, test_y1)
    spline2 = CubicSpline(test_x2, test_y2)
    
    x= np.linspace(-5,5,1000)
    
    y1 = spline1(x)
    y2 = spline2(x)

    fig = plt.figure()
    axs = fig.subplots(1,2).flat
    
    axs[0].plot(x,y1, label=r"S(x)")
    axs[0].plot(x,testFunc(x), label=r"f(x)")
    axs[0].scatter(test_x1,test_y1, color="black")
    axs[0].set_title(r"Low Grid Density", fontsize=title_sze)
    axs[0].set_xlabel(r"$\left(a\right)$", fontsize=lbl_sze)
    axs[0].legend()
    
    axs[1].plot(x,y2, label=r"S(x)")
    axs[1].plot(x,testFunc(x),label=r"f(x)")
    axs[1].scatter(test_x2,test_y2, color="black")
    axs[1].set_xlabel(r"$\left(b\right)$", fontsize=lbl_sze)
    axs[1].set_title(r"High Grid Density", fontsize=title_sze)

    
    
    fig.supxlabel(r'$x$', fontsize=lbl_sze)
    fig.supylabel(r'$y$', fontsize=lbl_sze) 
    fig.tight_layout()
    fig.savefig("intro_cubic.pdf", format="pdf")
    fig.show()
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["text.usetex"] = True

def chebyshevNodes(N, kind=1):
    nodes = np.empty(N)
    
    if kind==1:
        for i in np.arange(N):
            x = np.cos(np.pi*(2*i + 1)/(2*N))
            nodes[i] = x
    elif kind==2:
        for i in np.arange(N):
            x = np.cos(np.pi*i/(N-1))
            nodes[i] = x
    else:
        raise ValueError("Only first and second kind nodes exist")
    
    return nodes

def semicircle(r, h, k):
    x0 = h - r  # determine x start
    x1 = h + r  # determine x finish
    x = np.linspace(x0, x1, 10000)  # many points to solve for y

    # use numpy for array solving of the semicircle equation
    y = k + np.sqrt(r**2 - (x - h)**2)  
    return x, y



if __name__ == '__main__':
    ticksize = 16
    axsize = 24
    titlesize = 24
    n = 10
    
    roots = chebyshevNodes(n)
    roots2 = chebyshevNodes(n,2)
    y = np.sqrt(np.ones(n)-(roots**2))
    y2 = np.sqrt(np.ones(n)-(roots2**2))
    
    x1,circle1 = semicircle(1,0,0)
    x2,circle2 = semicircle(1,0,0)
    
    fig, (ax1,ax2) = plt.subplots(1,2)
    fig.set_size_inches(18,6)
    ax1.plot(x1,circle1, color="black", zorder=0)
    ax1.hlines(0,-1,1,color="black",linewidth=2, zorder=0)
    ax1.scatter(roots, np.zeros(n), label=r"nodes", zorder=1.5, color="blue")
    ax1.scatter(roots, y, label=r"unit circle projection", zorder=1, color="red")
    ax1.set_title(r"Chebyshev Zeroes for n={0}".format(n), fontsize=titlesize)
    ax1.set_xlim(-1.05,1.05)
    ax1.set_xlabel(r"(a)", fontsize=axsize)
    ax1.tick_params("both", labelsize=ticksize)
    ax1.set_ylim(-0.05,1.1)
    for i in np.arange(roots.shape[0]):
        ax1.vlines(roots[i], 0, y[i], color="black", linestyle="dashed", linewidth=1, zorder=0)
    ax2.plot(x2,circle2, color="black", zorder=0)
    ax2.hlines(0,-1,1,color="black", linewidth=2, zorder=0)
    ax2.scatter(roots2, np.zeros(n), label=r"nodes", zorder=1.5, color="blue",)
    ax2.scatter(roots2, y2, label=r"unit circle projection", zorder=1, color="red")
    ax2.set_xlim(-1.05,1.05)
    ax2.set_ylim(-0.05,1.1)
    ax2.set_title(r"Chebyshev Extrema for n={0}".format(n), fontsize=titlesize)
    ax2.set_xlabel(r"(b)", fontsize=axsize)
    ax2.tick_params("both", labelsize=ticksize)
    for i in np.arange(roots.shape[0]):
        ax2.vlines(roots2[i], 0, y2[i], color="black", linestyle="dashed", linewidth=1, zorder=0)
    fig.tight_layout()
    fig.savefig("roots.pdf", format="pdf")
    
    
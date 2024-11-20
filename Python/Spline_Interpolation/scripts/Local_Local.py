import numpy as np 
import matplotlib.pyplot as plt
from bisect import bisect_left
from scipy.stats import lognorm
from time import time

plt.rcParams["text.usetex"] = True

pi = np.pi

def gauss(x, y):
    return (np.cos(y) * np.sin(x) + 1.5)

def constructor(xs, ys):

    zs = np.empty([ys.size, xs.size])

    for iy, y in enumerate(ys):
        for ix, x in enumerate(xs):
            zs[iy, ix] = gauss(x, y)

    return zs

def local_local(x,y,z,xs,ys):
    ## here !!!!!!
    dsx,_ = np.gradient(z)
    m = xs.shape[0]
    n = ys.shape[0]
    val = np.empty((m,n))
    
    ixs = np.searchsorted(x, xs, side='right') - 1
    if xs[-1] >= x[-1]:
            ixs[-1] -= 1
            
    iys = np.searchsorted(y, ys, side='right') - 1
    if ys[-1] >= y[-1]:
            iys[-1] -= 1
    
    for i in np.arange(m):
        for j in np.arange(n):
            if (iys[j] == 0):
                # Only 3 points at left side -- only use first order difference for the left derivative
                z0 = cubic_interpolation(xs[i],x[ixs[i]],x[ixs[i]+1], z[ixs[i],iys[j]], z[ixs[i]+1,iys[j]], dsx[ixs[i],iys[j]], dsx[ixs[i]+1,iys[j]])
                z1 = cubic_interpolation(xs[i],x[ixs[i]],x[ixs[i]+1], z[ixs[i],iys[j]+1], z[ixs[i]+1,iys[j]+1], dsx[ixs[i],iys[j]+1], dsx[ixs[i]+1,iys[j]+1])
                z2 = cubic_interpolation(xs[i],x[ixs[i]],x[ixs[i]+1], z[ixs[i],iys[j]+2], z[ixs[i]+1,iys[j]+2], dsx[ixs[i],iys[j]+2], dsx[ixs[i]+1,iys[j]+2])
                
                dsy = np.gradient([z0,z1,z2])
                
                val[i,j] = cubic_interpolation(ys[j], y[iys[j]], y[iys[j]+1], z0, z1, dsy[0], dsy[1])
                
                
            elif (iys[j] >= y.shape[0]-2):
                # Same here -- first order for the right derivative as there isnt enough points
                z0 = cubic_interpolation(xs[i],x[ixs[i]],x[ixs[i]+1], z[ixs[i],iys[j]-1], z[ixs[i]+1,iys[j]-1], dsx[ixs[i],iys[j]-1], dsx[ixs[i]+1,iys[j]-1])
                z1 = cubic_interpolation(xs[i],x[ixs[i]],x[ixs[i]+1], z[ixs[i],iys[j]], z[ixs[i]+1,iys[j]], dsx[ixs[i],iys[j]], dsx[ixs[i]+1,iys[j]])
                z2 = cubic_interpolation(xs[i],x[ixs[i]],x[ixs[i]+1], z[ixs[i],iys[j]+1], z[ixs[i]+1,iys[j]+1], dsx[ixs[i],iys[j]+1], dsx[ixs[i]+1,iys[j]+1])
                
                dsy = np.gradient([z0,z1,z2])
                
                val[i,j] = cubic_interpolation(ys[j], y[iys[j]], y[iys[j]+1], z0, z1, dsy[1], dsy[2])
                 
            else:
                z0 = cubic_interpolation(xs[i],x[ixs[i]],x[ixs[i]+1], z[ixs[i],iys[j]-1], z[ixs[i]+1,iys[j]-1], dsx[ixs[i],iys[j]-1], dsx[ixs[i]+1,iys[j]-1])
                z1 = cubic_interpolation(xs[i],x[ixs[i]],x[ixs[i]+1], z[ixs[i],iys[j]], z[ixs[i]+1,iys[j]], dsx[ixs[i],iys[j]], dsx[ixs[i]+1,iys[j]])
                z2 = cubic_interpolation(xs[i],x[ixs[i]],x[ixs[i]+1], z[ixs[i],iys[j]+1], z[ixs[i]+1,iys[j]+1], dsx[ixs[i],iys[j]+1], dsx[ixs[i]+1,iys[j]+1])
                z3 = cubic_interpolation(xs[i],x[ixs[i]],x[ixs[i]+1], z[ixs[i],iys[j]+2], z[ixs[i]+1,iys[j]+2], dsx[ixs[i],iys[j]+2], dsx[ixs[i]+1,iys[j]+2])
                
                dsy = np.gradient([z0,z1,z2,z3])
                
                val[i,j] = cubic_interpolation(ys[j], y[iys[j]], y[iys[j]+1], z0, z1, dsy[1], dsy[2])
    
    return val
    
def cubic_interpolation(x,x0,x1,y0,y1,m0,m1):

    a = ( (m0 + m1)*(x0 - x1) - 2*(y0 - y1)) / ((x0 - x1)**3 ) 
    b = ( -m0*(x0 - x1)*(x0 + 2*x1) + m1*(-2*x0**2 + x1*x0 + x1**2) + 3*(x0 + x1)*(y0 - y1) ) / ( (x0 - x1)**3 )
    c = ( m1*x0*(x0 - x1)*(x0 + 2*x1) - x1*( m0*(-2*x0**2 + x1*x0 + x1**2) + 6*x0*(y0 - y1) ) ) / ( (x0 - x1)**3 )
    d = ( y1*(x0**2)*(x0 - 3*x1) + x1*( x0*(x1 - x0)*(m1*x0 + m0*x1) - x1*y0*(x1 - 3*x0) ) ) / ((x0 - x1)**3)
        
    return a*(x**3) + b*(x**2) + c*(x) + d

if __name__ == '__main__':
    x = np.linspace(-5,5,51)
    y = np.linspace(-5,5,51)
    z = constructor(x,y)
    
    x_test = np.random.uniform(-5,5,100)
    y_test = np.random.uniform(-5,5,100)
    
    N = 1000
    
    times = np.empty(N)
    
    for i in range(N):
        start = time()
        z_local = local_local(x,y,z,x_test,y_test)
        end = time()
        times[i] = end-start
    
    z_gauss = constructor(x_test,y_test)
    
    fig = plt.figure()
    axs = fig.subplots(1,2)
    
    plot1 = axs[0].contourf(x_test,y_test,z_local)
    plot2 = axs[1].contourf(x_test,y_test,z_gauss)
    
    axs[0].set_title(r"local cubic")
    axs[1].set_title(r"analytic")
    
    cbar1 = fig.colorbar(plot1, ax=axs[0])
    cbar2 = fig.colorbar(plot2, ax=axs[1])
    
    plt.savefig("../outputs/local_local.pdf", format="pdf")
    
    plt.figure()
    
    # shape, loc, scale = lognorm.fit(times)
    
    n, bins, patches = plt.hist(times,25, density=True)
    mean = np.mean(times)
    plt.axvline(mean, ls="dashed", color="black", label=r"mean time={0:.4f}".format(mean))

    # x = np.linspace(np.min(times), np.max(times), 1000)
    # p = lognorm.pdf(x,shape, loc, scale)
    
    # plt.plot(x,p,label=r'Log-normal Fit')
    #plt.xlim(0.5,0.7)

    plt.xlabel(r'Time(s)')
    plt.ylabel(r'Density')
    plt.legend()
    plt.savefig(r"../outputs/local_time_histo_1000")
    
    
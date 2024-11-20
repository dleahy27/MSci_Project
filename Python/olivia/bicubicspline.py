#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from cubicspline import CubicSpline

# latex font configs
plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['text.usetex'] = True


class BicubicSpline:
    '''
    Defines a smooth bicubic spline across a 2D rectangular grid of data points.
    '''

    def __init__(self, xs, ys, zs, d2x2s_left=0, d2x2s_right=0, d2y2s_bottom=0, d2y2s_top=0, d4x2y2s_corners=0, fname='PPEplots/bicubictestplots/bicubictest', plot1d=False, index_x=0, index_y=0, verbose=True):
        '''
        Finds bicubic spline coefficients for given dataset.
        Each rectangular cell on the data grid requires 16
        coefficients.

        params:   xs: array of x-coordinates of spline knots (horizontal grid dimension)
                  ys: array of y-coordinates of spline knots (vertical grid dimension)
                  zs: 2D array of function values at spline knots (x, y)

        kwargs:   BOUNDARY CONDITIONS:
                  d2x2s_left:       array of second derivatives in x on the left
                                    edge of the grid (default 0)
                  d2x2s_right:      array of second derivatives in x on the right
                                    edge of the grid (default 0)
                  d2y2s_bottom:     array of second derivatives in y on the bottom
                                    edge of the grid (default 0)
                  d2y2s_top:        array of second derivatives in y on the top
                                    edge of the grid (default 0)
                  d4x2y2s_corners:  array of xxyy-derivatives at the 4 corners of
                                    the grid in order bottom-left; bottom-right; 
                                    top-right; top-left (default 0)
                  OTHER:
                  fname:            file path for plots, without the file type
                  plot1d:           additional 1D spline plots made along specified
                                    indices on grid
                  index_x:          plot 1D spline in x along y=index_x
                  index_y:          plot 1D spline in y along x=index_y
        '''

        if verbose:
            print('in constructor:')
        
        # define internal variables
        self.xs = xs
        self.ys = ys
        self.zs = zs

        # grid dimensions
        self.m = self.xs.size
        self.n = self.ys.size

        if verbose:
            print('Performing checks and initialising variables...')

        # check if list of function values has the appropriate length
        if self.zs.shape != (self.n, self.m):
            raise ValueError('Array of function values must have shape (m, n), where m and n are the grid dimensions along x and y, respectively!')
        
        
        # set up arrays of boundary conditions

        ## check if the user has given condition
        if type(d2x2s_left) == np.ndarray:

            ## check if given list has correct size
            if d2x2s_left.size != self.n:
                raise ValueError('Left edge boundary array must have the same length as array of y coordinates!')
            
            ## set internal variable
            else: self.d2x2s_left = d2x2s_left

        ## otherwise, set boundary list to array of zeros
        else: self.d2x2s_left = np.zeros(self.n)

        ## perform similar checks on all other boundary conditions
        if type(d2x2s_right) == np.ndarray:
            if d2x2s_right.size != self.n:
                raise ValueError('Right edge boundary array must have the same length as array of y coordinates!')
            else: self.d2x2s_right = d2x2s_right
        else: self.d2x2s_right = np.zeros(self.n)

        if type(d2y2s_top) == np.ndarray:
            if d2y2s_top.size != self.m:
                raise ValueError('Top edge boundary array must have the same length as array of x coordinates!')
            else: self.d2y2s_top = d2y2s_top
        else: self.d2y2s_top = np.zeros(self.m)

        if type(d2y2s_bottom) == np.ndarray:
            if d2y2s_bottom.size != self.m:
                raise ValueError('Bottom edge boundary array must have the same length as array of x coordinates!')
            else: self.d2y2s_bottom = d2y2s_bottom
        else: self.d2y2s_bottom = np.zeros(self.m)

        if type(d4x2y2s_corners) == np.ndarray:
            if d4x2y2s_corners.size != 4:
                raise ValueError('Corner boundary array must have length 4!')
            else: self.d4x2y2s_corners = d4x2y2s_corners
        else: self.d4x2y2s_corners = np.zeros(4)


        # initialise matrices

        ## matrix of spline coefficients
        S = np.empty([self.n-1, self.m-1, 16])

        ## matrix of x derivatives at knots
        zs_x = np.empty([self.n, self.m])

        ## matrix of y derivatives at knots
        zs_y = np.empty([self.n, self.m])

        ## matrix of xy derivatives at knots
        zs_xy = np.empty([self.n, self.m])

        ## boundary vectors of xyy derivatives at top and bottom knots
        zs_xyy_bottom = np.empty([self.m])
        zs_xyy_top = np.empty([self.m])

        
        if verbose:
            print('Constructing 1D splines along each grid line...')

        # construct 1D splines along x and y to find first 
        # and cross derivatives at spline knots

        ## construct splines along boundary rows
        spline_bottom = CubicSpline(self.xs, self.zs[0], d2t2_start=self.d4x2y2s_corners[0], d2t2_end=self.d4x2y2s_corners[1])
        spline_top = CubicSpline(self.xs, self.zs[-1], d2t2_start=self.d4x2y2s_corners[3], d2t2_end=self.d4x2y2s_corners[2])
        
        
        if plot1d:
            
            if verbose:
                print('Plotting 1D splines along top and bottom boundary conditions...')

            spline_bottom.plotSpline(fname=fname+'_1Dspline_bottom', scatter=True)
            spline_top.plotSpline(fname=fname+'_1Dspline_top', scatter=True)

        ## get spline coefficients
        S_bottom = spline_bottom.S
        S_top = spline_top.S

        ## calculate and store analytical x derivatives at boundary knots
        ## the rest of the terms of the derivative are zero at the knots
        ## due to (x - self.t) terms (see CubicSpline code)
        zs_xyy_bottom[:-1] = S_bottom[:, 1]
        zs_xyy_bottom[-1] = S_bottom[-1, 1] + 2 * S_bottom[-1, 2] * (self.xs[-1] - self.xs[-2]) + 3 * S_bottom[-1, 3] * (self.xs[-1] - self.xs[-2])**2
        zs_xyy_top[:-1] = S_top[:, 1]
        zs_xyy_top[-1] = S_top[-1, 1] + 2 * S_top[-1, 2] * (self.xs[-1] - self.xs[-2]) + 3 * S_top[-1, 3] * (self.xs[-1] - self.xs[-2])**2

        print(f"\n\n\n\n\n")
        print(zs_xyy_bottom)
        print(f"\n\n\n\n\n")
        print(zs_xyy_top)
        print(f"\n\n\n\n\n")
        ## loop over data rows
        for j in range(self.n):
            
            ## construct spline along ith row using the 
            ## user input data (excluding the boundary rows)
            spline_x = CubicSpline(self.xs, self.zs[j], d2t2_start=self.d2x2s_left[j], d2t2_end=self.d2x2s_right[j])

            ## get spline coefficients
            S_x = spline_x.S

            ## calculate and store analytical x derivatives 
            ## at knots along the ith row
            zs_x[j, :-1] = S_x[:, 1]
            zs_x[j, -1] = S_x[-1, 1] + 2 * S_x[-1, 2] * (self.xs[-1] - self.xs[-2]) + 3 * S_x[-1, 3] * (self.xs[-1] - self.xs[-2])**2

            if plot1d:
                
                if verbose:
                    print('Plotting 1D spline along specified row...')

                if j == index_x:
                    spline_x.plotSpline(fname=fname+'_1Dspline_x', scatter=True)


        ## loop over data columns
        for i in range(self.m):
            
            ## construct spline along jth column using the calculated 
            ## spline derivatives (excluding boundary columns)
            spline_xy = CubicSpline(self.ys, zs_x[:, i], d2t2_start=zs_xyy_bottom[i], d2t2_end=zs_xyy_top[i])

            ## get spline coefficients
            S_xy = spline_xy.S

            ## calculate and store analytical xy derivatives
            ## at knots along the jth column
            zs_xy[:-1, i] = S_xy[:, 1]
            zs_xy[-1, i] = S_xy[-1, 1] + 2 * S_xy[-1, 2] * (self.ys[-1] - self.ys[-2]) + 3 * S_xy[-1, 3] * (self.ys[-1] - self.ys[-2])**2
            
            ## construct spline along jth column using the user input data
            ## (excluding boundary columns since they have no further use)
            spline_y = CubicSpline(self.ys, self.zs[:, i], d2t2_start=self.d2y2s_bottom[i], d2t2_end=self.d2y2s_top[i])

            ## get spline coefficients
            S_y = spline_y.S

            ## calculate and store analytical xy derivatives
            ## at knots along the jth column
            zs_y[:-1, i] = S_y[:, 1]
            zs_y[-1, i] = S_y[-1, 1] + 2 * S_y[-1, 2] * (self.ys[-1] - self.ys[-2]) + 3 * S_y[-1, 3] * (self.ys[-1] - self.ys[-2])**2

            if plot1d:

                if verbose:
                    print('Plotting 1D splines along specified column...')

                if i == index_y:
                    spline_xy.plotSpline(fname=fname+'_1Dspline_xy', scatter=True)
                    spline_y.plotSpline(fname=fname+'_1Dspline_y', scatter=True)

        
        if verbose:
            print('Filling spline coefficient matrix...')

        # fill knot values and derivatives into matrix

        ## loop over the outer dimensions of the matrix
        for j in range(self.n-1):
            for i in range(self.m-1):

                ## array of knot values and knot derivatives
                ## at each iteration (grid cell)
                Z = np.empty([16])

                ## function values
                Z[0] = self.zs[j, i]
                Z[1] = self.zs[j, i+1]
                Z[2] = self.zs[j+1, i]
                Z[3] = self.zs[j+1, i+1]

                ## x derivatives
                Z[4] = zs_x[j, i]
                Z[5] = zs_x[j, i+1]
                Z[6] = zs_x[j+1, i]
                Z[7] = zs_x[j+1, i+1]

                ## y derivatives
                Z[8] = zs_y[j, i]
                Z[9] = zs_y[j, i+1]
                Z[10] = zs_y[j+1, i]
                Z[11] = zs_y[j+1, i+1]

                ## xy derivatives
                Z[12] = zs_xy[j, i]
                Z[13] = zs_xy[j, i+1]
                Z[14] = zs_xy[j+1, i]
                Z[15] = zs_xy[j+1, i+1]

                if verbose:
                    if i == 0 and j == 0:
                        print('Here is a matrix Z containing the function values and x, y, xy derivatives at the 4 knots of the first rectangle:')
                        print(Z)

                ## matrix A such that A * S[j, i] = Z
                A = np.empty([16, 16])

                ## loop over dimensions of A
                for l in range(4):
                    for k in range(4):

                        ## function value terms
                        A[0, 4*l+k] = self.xs[i]**k * self.ys[j]**l
                        A[1, 4*l+k] = self.xs[i+1]**k * self.ys[j]**l
                        A[2, 4*l+k] = self.xs[i]**k * self.ys[j+1]**l
                        A[3, 4*l+k] = self.xs[i+1]**k * self.ys[j+1]**l

                        # initialise remaining terms to zero
                        A[4, 4*l+k] = A[5, 4*l+k] = A[6, 4*l+k] = A[7, 4*l+k] = 0
                        A[8, 4*l+k] = A[9, 4*l+k] = A[10, 4*l+k] = A[11, 4*l+k] = 0
                        A[12, 4*l+k] = A[13, 4*l+k] = A[14, 4*l+k] = A[15, 4*l+k] = 0

                        if k > 0:

                            ## x derivative terms
                            A[4, 4*l+k] = k * self.xs[i]**(k-1) * self.ys[j]**l
                            A[5, 4*l+k] = k * self.xs[i+1]**(k-1) * self.ys[j]**l
                            A[6, 4*l+k] = k * self.xs[i]**(k-1) * self.ys[j+1]**l
                            A[7, 4*l+k] = k * self.xs[i+1]**(k-1) * self.ys[j+1]**l

                            if l > 0:

                                ## y derivative terms
                                A[8, 4*l+k] = l * self.xs[i]**k * self.ys[j]**(l-1)
                                A[9, 4*l+k] = l * self.xs[i+1]**k * self.ys[j]**(l-1)
                                A[10, 4*l+k] = l * self.xs[i]**k * self.ys[j+1]**(l-1)
                                A[11, 4*l+k] = l * self.xs[i+1]**k * self.ys[j+1]**(l-1)

                                ## xy derivative terms
                                A[12, 4*l+k] = k * l * self.xs[i]**(k-1) * self.ys[j]**(l-1)
                                A[13, 4*l+k] = k * l * self.xs[i+1]**(k-1) * self.ys[j]**(l-1)
                                A[14, 4*l+k] = k * l * self.xs[i]**(k-1) * self.ys[j+1]**(l-1)
                                A[15, 4*l+k] = k * l * self.xs[i+1]**(k-1) * self.ys[j+1]**(l-1)

                        elif l > 0:

                            ## remaining y derivative terms
                            A[8, 4*l+k] = l * self.xs[i]**k * self.ys[j]**(l-1)
                            A[9, 4*l+k] = l * self.xs[i+1]**k * self.ys[j]**(l-1)
                            A[10, 4*l+k] = l * self.xs[i]**k * self.ys[j+1]**(l-1)
                            A[11, 4*l+k] = l * self.xs[i+1]**k * self.ys[j+1]**(l-1)
                
                if verbose:
                    if i == 0 and j == 0:
                        print('Here is the corresponding matrix A:')
                        print(A)
                
                ## find determinant of A
                A_det = np.linalg.det(A)

                if verbose:
                    if i == 0 and j == 0:
                        print('and its determinant:')
                        print(A_det)

                ## invert A 
                A_inv = np.linalg.inv(A)

                if verbose:
                    if i == 0 and j == 0:
                        print('and its inverse:')
                        print(A_inv)

                ## calculate the coefficient array for this iteration (grid cell)
                ## meaning: S[j, i] = A^-1 * Z
                S[j, i] = np.linalg.matmul(A_inv, Z)

                if verbose:
                    if i == 0 and j == 0:
                        print('The corresponding spline coefficients then read:')
                        print(S[j, i])

        
        # define callable attribute
        self.S = S

        if verbose:
            print('Finally, this is the full spline coefficient matrix:')
            #print(self.S)
        
        # automatic print statement
        print('Spline coefficients have been successfully calculated and stored!')


    def evaluateSpline(self, x, y, verbose=True):
        '''
        Returns the value(s) of the interpolated function at x, y.

        params:  x:  x-coordinates(s) at which to evaluate function.
                     User can provide either a single input 
                     value or an array of inputs.
                 y:  y-coordinate(s) at which to evaluate function.
                     User can provide either a single input 
                     value or an array of inputs.
        '''
        
        if verbose:
            print('in evaluateSpline():')

        # store, then convert to array if user input is a scalar
        x_check = x
        y_check = y

        if type(x) != np.ndarray:
            x = np.array([x])
        if type(y) != np.ndarray:
            y = np.array([y])

        # set internal xy variables
        x_eval = x
        y_eval = y

        # initialise array storing function values
        splines = np.empty([y_eval.size, x_eval.size])

        # find indices before which each user input value should be
        # slotted into the arrays of knot coordinates such that the
        # arrays remain sorted; these indices (ix, iy) correspond to
        # the spline indices (i, j) at which to evaluate the spline
        # for each input point

        if verbose:
            print('Determining grid locations...')

        ## relevant indices
        ## np.searchsorted uses binary search
        ## side='right': self.xs[i] <= x_eval < self.xs[i+1] as required
        ## indices given are those of the knot values before which each
        ## input value should be slotted in so need to be shifted back by one
        ixs = np.searchsorted(self.xs, x_eval, side='right') - 1
        iys = np.searchsorted(self.ys, y_eval, side='right') - 1

        ## if final input value equal to (or greater than) the final knot 
        ## value, rightmost value evaluates to final index of knot array + 1 
        ## so it needs to be shifted back by one more
        if x_eval[-1] >= self.xs[-1]:
            ixs[-1] -= 1
        if y_eval[-1] >= self.ys[-1]:
            iys[-1] -= 1

        if verbose:
            print('x indices:')
            print(ixs)
            print('y indices:')
            print(iys)

            print('Evaluating bicubic spline at each input coordinate...')

        # loop over input coordinates
        for iy, y in enumerate(y_eval):
            for ix, x in enumerate(x_eval):
                print(iys[iy],ixs[ix])
                # evaluate spline at (ix, iy)-th coordinates and store in array
                splines[iy, ix] = (   self.S[iys[iy], ixs[ix], 0]                + self.S[iys[iy], ixs[ix], 1]  * x
                                    + self.S[iys[iy], ixs[ix], 2]  * x**2        + self.S[iys[iy], ixs[ix], 3]  * x**3
                                    + self.S[iys[iy], ixs[ix], 4]         * y    + self.S[iys[iy], ixs[ix], 5]  * x    * y
                                    + self.S[iys[iy], ixs[ix], 6]  * x**2 * y    + self.S[iys[iy], ixs[ix], 7]  * x**3 * y
                                    + self.S[iys[iy], ixs[ix], 8]         * y**2 + self.S[iys[iy], ixs[ix], 9]  * x    * y**2
                                    + self.S[iys[iy], ixs[ix], 10] * x**2 * y**2 + self.S[iys[iy], ixs[ix], 11] * x**3 * y**2
                                    + self.S[iys[iy], ixs[ix], 12]        * y**3 + self.S[iys[iy], ixs[ix], 13] * x    * y**3
                                    + self.S[iys[iy], ixs[ix], 14] * x**2 * y**3 + self.S[iys[iy], ixs[ix], 15] * x**3 * y**3 )

        if verbose:
            print('Here is the full array of function values:')
            print(splines)
        
        # return function values at each user input coordinate (x, y)
        if type(x_check) != np.ndarray and type(y_check) != np.ndarray:
            return splines[0, 0]
        
        return splines
    

    def plotSpline(self, num=50, scatter=False, func=None, diff_rel=False, pdf=False, fl=None, scale='linear', fname='PPEplots/bicubictestplots/bicubictest', axes=['$x$', '$y$', '$z$'], verbose=False):
        '''
        Creates plot of bicubic spline.

        kwargs: num:       number of points to sample within given range
                scatter:   overlays a scatter of the original data points 
                           if set to True
                func:      also plots (relative) difference between the spline 
                           and the true function, func, for visual error 
                           analysis
                diff_rel:  if True, the relative difference between spline
                           and true function will be plotted, otherwise just
                           the difference will be plotted
                pdf:       if True, the true function will evaluate a pdf
                           on the existing lhapdf patch algorithm for comparison
                fl:        pdf flavour (pid)
                scale:     scales treated as logarithmic if set to 'log';
                           default 'linear'
                fname:     file path for plots, without the file type
                axes:      name of the plot axes; default 'x', 'y', 'z'
        '''

        if verbose:
            print('in plotSpline():')

        # define array of locations to evaluate spline at
        xs_plot = np.array(np.linspace(self.xs[0], self.xs[-1], num=num))
        ys_plot = np.array(np.linspace(self.ys[0], self.ys[-1], num=num))

        if verbose:
            print('Test coordinates for plotting:')
            print(xs_plot)
            print(ys_plot)

            print('Evaluating spline at test coordinates...')

        # evaluate spline at each xy coordinate
        zs_plot = self.evaluateSpline(xs_plot, ys_plot, verbose=verbose)

        if func != None:

            if verbose:
                print('Evaluating true function at test coordinates...')

            # evaluate given function at each xy coordinate if specified
            zs_true = np.empty([ys_plot.size, xs_plot.size])

            if pdf:

                for iy, y in enumerate(10**ys_plot):
                    for ix, x in enumerate(10**xs_plot):
                        zs_true[iy, ix] = func(fl, x, y)

            else:

                for iy, y in enumerate(ys_plot):
                    for ix, x in enumerate(xs_plot):
                        zs_true[iy, ix] = func(x, y)

            if verbose:
                print('Finding (relative) difference between spline and true function...')

            # find (relative) difference between spline and true function
            if diff_rel:

                diff = (zs_true - zs_plot) / zs_true

                if verbose:
                    print('Relative difference:')
                    print(diff)
            
            else:

                diff = (zs_true - zs_plot)

                if verbose:
                    print('Difference:')
                    print(diff)

        # create a meshgrid for plotting
        Xs_plot, Ys_plot = np.meshgrid(xs_plot, ys_plot)

        if func != None:

            if verbose:
                print('Plotting (relative) difference between spline and true function...')

            # initialise figure
            plt.figure(figsize=(10, 8))

            # plot contour map of difference
            # can also try: pcolor, pcolormesh
            if scale == 'log':
                plt.contourf(10**Xs_plot, 10**Ys_plot, diff, levels=12, cmap='plasma')
            else:
                plt.contourf(Xs_plot, Ys_plot, diff, levels=12, cmap='plasma')

            # set colorbar
            cbar = plt.colorbar()
            cbar.ax.set_ylabel(axes[2], fontsize=20, labelpad=15)

            # add scatter of original data if specified
            if scatter:
                Xs_data, Ys_data = np.meshgrid(self.xs, self.ys)

                if scale == 'log':
                    plt.scatter(10**Xs_data, 10**Ys_data, color='white')
                else:
                    plt.scatter(Xs_data, Ys_data, color='white')

            # plot configs
            plt.xlabel(axes[0], fontsize=20, labelpad=15)
            plt.ylabel(axes[1], fontsize=20, labelpad=15)

            if scale == 'log':
                plt.xscale('log')
                plt.yscale('log')

            plt.savefig(fname=fname+'_errors.pdf')

        if verbose:
            print('Plotting spline...')

        # initialise main figure
        plt.figure(figsize=(10, 8))

        # plot contour map of spline
        if scale == 'log':
            plt.contourf(10**Xs_plot, 10**Ys_plot, zs_plot, levels=12, cmap='plasma')
        else:
            plt.contourf(Xs_plot, Ys_plot, zs_plot, levels=12, cmap='plasma')

        # set colorbar
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(axes[2], fontsize=20, labelpad=15)

        # add scatter of original data if specified
        if scatter:
            Xs_data, Ys_data = np.meshgrid(self.xs, self.ys)

            if scale == 'log':
                plt.scatter(10**Xs_data, 10**Ys_data, color='white')
            else:
                plt.scatter(Xs_data, Ys_data, color='white')

        # plot configs
        plt.xlabel(axes[0], fontsize=20, labelpad=15)
        plt.ylabel(axes[1], fontsize=20, labelpad=15)

        if scale == 'log':
            plt.xscale('log')
            plt.yscale('log')

        plt.savefig(fname=fname+'.pdf')


    def plotSplineDerivative(self, order=1, num=50, scatter=False, scale='linear', fname='PPEplots/bicubictestplots/bicubictest', axes=['$x$', '$y$', '$z$'], plot1d=False, index_x=0, index_y=0, verbose=False):
        '''
        Creates plot of bicubic spline derivative.

        kwargs: order:     order of derivative (1 or 2)
                num:       number of points to sample within given range
                scatter:   overlays a scatter of the original data points 
                           if set to True
                scale:     scales treated as logarithmic if set to 'log';
                           default 'linear'
                fname:     file path for plots, without the file type
                axes:      name of the plot axes; default 'x', 'y', 'z'
                plot1d:    create additional 1D plots along given indices 
                           in x and y
                index_x:  plot spline y derivative along x=index_x
                index_y:  plot spline x derivative along y=index_y
        '''

        if verbose:
            print('in plotSplineDerivative():')

        # define array of locations to evaluate spline at
        xs_dplot = np.array(np.linspace(self.xs[0], self.xs[-1], num=num))
        ys_dplot = np.array(np.linspace(self.ys[0], self.ys[-1], num=num))

        if verbose:
            print('Test coordinates for plotting:')
            print(xs_dplot)
            print(ys_dplot)

            print('Evaluating spline at test coordinates...')

        # evaluate spline at each xy coordinate
        zs_dplot = self.evaluateSpline(xs_dplot, ys_dplot, verbose=verbose)

        if scale == 'log':
            xs_dplot = np.array(np.logspace(self.xs[0], self.xs[-1], num=num))
            ys_dplot = np.array(np.logspace(self.ys[0], self.ys[-1], num=num))

        if verbose:
            print('Calculating derivatives of spline at each test coordinate...')

        # calculate derivatives of spline at each xy coordinate
        if order == 1:
            zs_dxplot = np.gradient(zs_dplot, xs_dplot, axis=1)
            zs_dyplot = np.gradient(zs_dplot, ys_dplot, axis=0)

            if verbose:
                print('dx:')
                print(zs_dxplot)
                print('dy:')
                print(zs_dyplot)
        
        elif order == 2:
            zs_dxplot = np.gradient(np.gradient(zs_dplot, xs_dplot, axis=1), xs_dplot, axis=1)
            zs_dyplot = np.gradient(np.gradient(zs_dplot, ys_dplot, axis=0), ys_dplot, axis=0)

            if verbose:
                print('d2x:')
                print(zs_dxplot)
                print('d2y:')

        else: raise ValueError('Derivatives available only up to 2nd order.')
        
        # create a meshgrid for plotting
        Xs_dplot, Ys_dplot = np.meshgrid(xs_dplot, ys_dplot)


        if verbose:
            print('Plotting x derivatives...')

        # heatmap of second derivatives along x

        # initialise figure
        plt.figure(figsize=(10, 8))

        # plot contour map of spline second derivative
        plt.contourf(Xs_dplot, Ys_dplot, zs_dxplot, levels=12, cmap='plasma')

        # set colorbar
        cbar = plt.colorbar()
        if order == 1:
            cbar.ax.set_ylabel('$\mathrm{d} \;$'+axes[2], fontsize=20, labelpad=15)
        elif order == 2:
            cbar.ax.set_ylabel('$\mathrm{d}^2 \;$'+axes[2], fontsize=20, labelpad=15)

        # add scatter of original data if specified
        if scatter:
            Xs_data, Ys_data = np.meshgrid(self.xs, self.ys)

            if scale == 'log':
                plt.scatter(10**Xs_data, 10**Ys_data, color='white')
            else:
                plt.scatter(Xs_data, Ys_data, color='white')

        # plot configs
        plt.xlabel(axes[0], fontsize=20, labelpad=15)
        plt.ylabel(axes[1], fontsize=20, labelpad=15)

        if scale == 'log':
            plt.xscale('log')
            plt.yscale('log')

        if order == 1: 
            plt.savefig(fname=fname+'_dx.pdf')
        elif order == 2: 
            plt.savefig(fname=fname+'_d2x.pdf')


        if verbose:
            print('Plotting y derivatives...')

        # heatmap of second derivatives along y

        # initialise figure
        plt.figure(figsize=(10, 8))

        # plot contour map of spline second derivative
        plt.contourf(Xs_dplot, Ys_dplot, zs_dyplot, levels=12, cmap='plasma')

        # set colorbar
        cbar = plt.colorbar()
        if order == 1:
            cbar.ax.set_ylabel('$\mathrm{d} \;$'+axes[2], fontsize=20, labelpad=15)
        elif order == 2:
            cbar.ax.set_ylabel('$\mathrm{d}^2 \;$'+axes[2], fontsize=20, labelpad=15)

        # add scatter of original data if specified
        if scatter:
            Xs_data, Ys_data = np.meshgrid(self.xs, self.ys)

            if scale == 'log':
                plt.scatter(10**Xs_data, 10**Ys_data, color='white')
            else:
                plt.scatter(Xs_data, Ys_data, color='white')

        # plot configs
        plt.xlabel(axes[0], fontsize=20, labelpad=15)
        plt.ylabel(axes[1], fontsize=20, labelpad=15)

        if scale == 'log':
            plt.xscale('log')
            plt.yscale('log')

        if order == 1: 
            plt.savefig(fname=fname+'_dy.pdf')
        elif order == 2: 
            plt.savefig(fname=fname+'_d2y.pdf')


        if plot1d:

            if verbose:
                print('Plotting x derivatives along a specific row...')

            # 1D plot of x derivative along given y grid line

            # relevant slice of x derivative values
            zsdx_1D = zs_dxplot[index_y]
            
            # initialise figure
            plt.figure(figsize=(10, 8))

            # plot spline derivative in x
            plt.plot(xs_dplot, zsdx_1D, color='r')

            # plot configs
            plt.xlabel(axes[0], fontsize=20)
            if order == 1: 
                plt.ylabel('$\mathrm{d} \;$'+axes[2], fontsize=20)
            elif order == 2: 
                plt.ylabel('$\mathrm{d}^2 \;$'+axes[2], fontsize=20)

            if scale == 'log':
                plt.xscale('log')

            if order == 1: 
                plt.savefig(fname=fname+'_dx_1D.pdf')
            elif order == 2: 
                plt.savefig(fname=fname+'_d2x_1D.pdf')


            if verbose:
                print('Plotting y derivatives along a specific column...')

            # 1D plot of y derivative along given x grid line

            # relevant slice of y derivative values
            zsdy_1D = zs_dyplot[:, index_x]
            
            # initialise figure
            plt.figure(figsize=(10, 8))

            # plot spline derivative in x
            plt.plot(ys_dplot, zsdy_1D, color='r')

            # plot configs
            plt.xlabel(axes[1], fontsize=20)
            if order == 1: 
                plt.ylabel('$\mathrm{d} \;$'+axes[2], fontsize=20)
            elif order == 2: 
                plt.ylabel('$\mathrm{d}^2 \;$'+axes[2], fontsize=20)

            if scale == 'log':
                plt.xscale('log')

            if order == 1: 
                plt.savefig(fname=fname+'_dy_1D.pdf')
            elif order == 2: 
                plt.savefig(fname=fname+'_d2y_1D.pdf')
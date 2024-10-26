#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# latex font configs
plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['text.usetex'] = True


class CubicSpline:
    '''
    Defines a smooth cubic spline across a 1D line of data points.
    '''

    def __init__(self, ts, ys, d2t2_start=0, d2t2_end=0):
        '''
        Finds cubic spline coefficients for given dataset.
        Each spline segment requires 4 coefficients.

        params:   ts: array of spline knots
                  ys: array of function values at spline knots

        kwargs:   d2t2_start: second derivative at first knot
                              (default 0)
                  d2t2_end:   second derivative at final knot
                              (default 0)
        '''

        # define internal variables
        self.ts = ts
        self.ys = ys

        # total number of knots
        self.n = self.ts.shape[0]

        # check if the two lists are the same size
        if self.ys.shape[0] != self.n:
            raise ValueError('Arrays must have the same length!')
        

        # initialise arrays

        ## useful values for matrix calculation
        hs = np.empty([self.n-1])
        bs = np.empty([self.n-1])
        us = np.empty([self.n-1])
        vs = np.empty([self.n-1])

        ## S''(ti) values
        zs = np.empty([self.n])

        ## coefficients of each sub-spline,
        ## where sub-splines Si have the form:
        ## Si(x) = Ai + Bi(x-ti) + Ci(x-ti)^2 + Di(x-ti)^3
        ## which is a Taylor expansion of Si about ti
        S = np.empty([self.n-1, 4])

        ## some initial array elements
        us[0] = vs[0] = 0
        zs[0] = d2t2_start
        zs[-1] = d2t2_end


        # loop over size of arrays
        for i in range(self.n-1):

            ## define useful values for matrix calculation
            hs[i] = self.ts[i+1] - self.ts[i]
            bs[i] = (self.ys[i+1] - self.ys[i]) / hs[i]

            ## forward elimination
            if i == 1:
                us[i] = 2 * (hs[i] + hs[i-1])
                vs[i] = 6 * (bs[i] - bs[i-1])

            if i > 1:
                us[i] = 2 * (hs[i] + hs[i-1]) - ( np.power(hs[i-1], 2) / us[i-1] )
                vs[i] = 6 * (bs[i] - bs[i-1]) - ( (hs[i-1] * vs[i-1]) / us[i-1] )


        # backwards loop over size of arrays
        for i in range(self.n-2, -1, -1):

            ## back substitution
            if (i > 0):

                ## the following works because zs[n-1] has been set to a constant
                zs[i] = ( vs[i] - (hs[i] * zs[i+1]) ) / us[i]

                ## find coefficients for each sub-spline
                S[i, 0] = self.ys[i]
                S[i, 1] = - (hs[i] * zs[i+1]) / 6 - (hs[i] * zs[i]) / 3 + bs[i]
                S[i, 2] = zs[i] / 2
                S[i, 3] = (zs[i+1] - zs[i]) / (6 * hs[i])

            ## find remaining coefficients
            else:
                S[i, 0] = self.ys[i]
                S[i, 1] = - (hs[i] * zs[i+1]) / 6 - (hs[i] * zs[i]) / 3 + bs[i]
                S[i, 2] = zs[i] / 2
                S[i, 3] = (zs[i+1] - zs[i]) / (6 * hs[i])

        
        # define callable coefficient variable
        self.S = S


    def evaluateSpline(self, x, verbose=False):
        '''
        Returns the value(s) of the interpolated function at x.

        params:  x:  value(s) at which to evaluate function.
                     User can provide either a single input 
                     value or an array of inputs.
        '''
        
        # initialise variables (define empty lists because np.empty()
        # needs to know the dtype in advance, but dtype could be
        # list of bools or bool for conditions and function for
        # subsplines)
        conditions = [] 
        subsplines = []

        # loop over spline coefficients
        for i in range(self.n-1):
            
            # define condition for ith sub-spline: each condition is a 
            # boolean list defining the location along the spline of 
            # each user input value

            ## at end of spline, value can also be equal to final knot
            if i == self.n-2:
                cond = np.logical_and(x >= self.ts[i], x <= self.ts[i+1])
                conditions.append(cond)

            ## otherwise
            else:
                cond = np.logical_and(x >= self.ts[i], x < self.ts[i+1])
                conditions.append(cond)

            # define ith sub-spline
            subspline = lambda y, ii=i: ( self.S[ii, 0] + self.S[ii, 1] * (y - self.ts[ii]) 
                                         + self.S[ii, 2] * np.power(y - self.ts[ii], 2) 
                                         + self.S[ii, 3] * np.power(y - self.ts[ii], 3) )
            subsplines.append(subspline)
            
        # piecewise define the full spline and evaluate 
        # at input values
        return np.piecewise(x, conditions, subsplines)
        
    
    def plotSpline(self, scatter=False, func=None, fname='PPEplots/cubictestplots/cubictest'):
        '''
        Creates plot of cubic spline.

        kwargs: scatter: overlays a scatter of the original
                         data points if set to True.
                func:    overlays a plot of the true function, 
                         func, to enable visual comparison between 
                         it and the spline function
                fname:   file path for plots, without the file type
        '''

        # define array of locations to evaluate spline at
        xs = np.array(np.linspace(self.ts[0], self.ts[-1], 500))

        # call function at array locations
        ys = self.evaluateSpline(xs)

        # initialise figure
        plt.figure(figsize=(10, 8))

        # plot spline
        plt.plot(xs, ys, color='r', label='cubic spline')

        # add scatter of original data if specified
        if scatter:
            plt.scatter(self.ts, self.ys, color='k', label='knots')

        # overlay true function if specified
        if func != None:
            true_ys = func(xs)
            plt.plot(xs, true_ys, color='b', linestyle='--', label='true function')

        # plot configs
        plt.xlabel('$x$', fontsize=20)
        plt.ylabel('$y$', fontsize=20)
        plt.legend(fontsize=15)
        plt.savefig(fname=fname+'.pdf')
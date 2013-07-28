"""A brute force approach to finding the error of a fitting procedure
by refitting the raw data with noise added.
Does not rely on the assumption of small errors, but does assume that
errors are uncorrelated.
Stephen Inglis, 2013"""

import numpy
from scipy.optimize import leastsq

def refit(x,y,dy,errf,c_guess,N_steps=100):
    """Main module to fit the data, assuming one-dimensional data.
x - dependant variable of some data
y - independant variable
dy - standard deviation of the mean of the independant data
erry(C,x,y,dy) - error function of the fit, typically of the form (fitf(...) - y)/dy
c_guess - list of intial guess of the solution C for fitf(...), [c1,c2,...]
N_steps - number of trials of the noise adding procedure, default 100"""
    # Ensure the data is in numpy arrays, for easier manipulation
    x = numpy.array(x)
    y = numpy.array(y)
    dy = numpy.array(dy)
    psum = None
    for _i in range(N_steps):
        # Generate data with noise, assuming each data point is uncorrelated
        ny = numpy.array([numpy.random.normal(loc=z[0],scale=z[1]) for z in zip(y,dy)])
        # Main workhorse from scipy, minimizing the sum of the squares from 'errf' by modifying the first parameters of the function
        # Insert your favourite minimization routine here
        res = leastsq(errf, c_guess, args=(x,ny,dy), full_output=1)
        (popt, pcov, infodict, errmsg, ier) = res
        if psum == None:
            psum = numpy.array(popt)
            psum2 = numpy.array(popt)**2
        else:
            psum += numpy.array(popt)
            psum2 += numpy.array(popt)**2
    psum /= N_steps
    psum2 /= N_steps
    perr = (psum2 - psum**2)
    # In the case of data without error, this prevents errors
    perr = (perr*(perr>0))**0.5
    return psum,perr

def main_test():
    """A simple experiment.
Imagine we have 3 independent methods for calculating a value X.
We assume no systematic bias, but we assume they have different accuracy
Let us see how using a flat average and weighted average compare for finding the true mean."""
    x = numpy.array([1.,2.,3.])
    error = numpy.array([1.,5.,25.])
    numtrials = 10
    y = [0,0,0]
    dy = [0,0,0]
    for i in range(3):
        temp = numpy.random.normal(loc=10.,scale=error[i],size=numtrials)
        y[i] = numpy.mean(temp)
        dy[i] = numpy.std(temp)/(numtrials**0.5)
    # 1-parameter fit, the average
    fitf = lambda fc,fx: fc[0]
    print main_test.__doc__
    print 'Real mean = 10.00'
    print 'Meaningful Error < %f' % dy[0]
    # Using a weighted average
    errf = lambda fc,fx,fy,fdy: (fitf(fc,fx) - fy)/fdy
    popt,perr = refit(x,y,dy,errf,[10])
    print 'Weighted Mean = %f +- %f' % (popt[0],perr[0])
    # Using an unweighted (flat) average
    errf2 = lambda fc,fx,fy,fdy: (fitf(fc,fx) - fy)
    popt,perr = refit(x,y,dy,errf2,[10])
    print 'Flat Mean = %f +- %f' % (popt[0],perr[0])

if __name__ == "__main__":
    main_test()

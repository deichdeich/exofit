import numpy as np
import matplotlib.pyplot as plt
import TransitMaker as tm
import math
import time


plt.ion()

### FAKE DATA ###
### the first option is just the exact curve
### the second option is the same curve with
### Gaussian random error.
#data = get_curve(np.array(real_y))
data = np.genfromtxt('fakedata.csv')


### FIT PARAMETERS ###

"""
actual values from which the fake data were derived
(i.e. the script should recover these):
e = .0	
aRs = 15.23								
i = 1.555													
p0 = 0.115
period = 0.58040464894       								
"""
real_y = np.matrix([0,15.23,1.55,0.115,0.58040464894])

### INITIAL GUESSES ###
e = 1e-5 #eccentricity
aRs = 14.5 #semi-major axis to stellar radius ratio							
i = 1.4 #inclination angle (radians)					
p0 = 0.1 #planet to star radius ratio									
period = 0.5 #planet orbital period (days)      



### OTHER PARAMETERS ###
			
t0 = 2454966.52507 #ephemeris (JD)
eps = 1.0e-7 #minimum eccentricity for solving Kepler's eqn
             #(otherwise circular orbit assumed)  					            
w0 = math.pi/2. #argument of periapse
u1 = 0.1 #quadratic limb darkening parameter 1
u2 = 0.3 #quadratic limb darkening parameter 2


### FIT STUFF ###
"""
These are given the same variable names as Algorithm 12.3
on p.312 in Joel Franklin's Methods for Computational Physics.
"""
t = np.linspace(t0-period/20., t0 + period/20., 1000)
statics = np.array([u1,u2,w0,t0,eps])[np.newaxis]



### Returns light curve.  'statics' array has values I'm not fitting for
def get_curve(x):
    curve = tm.get_lc(x[0][0],
                      x[0][1],
                      x[0][2],
                      x[0][3],
                      x[0][4],
                      statics[0][0],
                      statics[0][1],
                      statics[0][2],
                      statics[0][3],
                      statics[0][4],
                      t)
    return curve



### Returns the u(x) scalar
def u(x):
    curve_x = get_curve(x) 
    diffs = (data-curve_x)**2
    total = np.sum(diffs)
    return total


# Calculate the derivative at the point y
# with respect to the variable wrt.  wrt
# is the index of that component in y.
### SUGGESTION: steepest-descent.py is crapping out
### for inputs far from the actual fit values.
### This is probably because the gradient of the function
### near these input values is zero.  So put something in here
### to say that if the gradient is zero, increase
### the dx until it's not zero.
def gradient(a):
    dx = np.array([1e-7,0.1,0.01,0.001,0.001])
    df = np.zeros(len(np.array(a)[0]))
    for wrt in xrange(len(dx)):
        plus_a,min_a = np.array(a), np.array(a)
        plus_a[0][wrt] += dx[wrt]
        min_a[0][wrt] -= dx[wrt]
        plus_u = u(np.array(plus_a))
        min_u = u(np.array(min_a))
        difference = plus_u-min_u
        dfdwrt = difference/(2*dx[wrt])
        df[wrt] = dfdwrt
    return np.linalg.norm(df)

# make_plot takes:
# plotnum: (0,1)
# drawables is an array with ([data],[actual values],[fit values])
# y are the current input parameters
# i is the iteration
# nrm is the normalized 
def make_plot(plotnum, drawables, y, i, nrm):
    plt.figure(plotnum)
    plt.cla()
    if plotnum == 0:
        plt.plot(t,get_curve(np.array(x0)),color="green")
        plt.plot(t,get_curve(np.array(y)),color='red')
        plt.scatter(t,data)
        plt.ylim((0.975, 1.005))
        y = np.array(y)
        if np.isnan(np.sum(y)):
            nrm = nrm0
            y = y0
            print 'nan values'
        else:
            nrm0 = nrm
            y0 = y
        textstr = 'iteration %.0f e=%.6f aRs=%.6f i=%.6f\n p0 = %.6f period = %.6f' %(i,
                                                                        y[0][0],
                                                                        y[0][1],
                                                                        y[0][2],
                                                                        y[0][3],
                                                                        y[0][4])
        plt.title(textstr)
        plt.draw()
    
    elif plotnum == 1:
        plt.figure(1)
        plt.title('Norm of residuals: %.6f'%nrm)
        if nrm > ylimit-0.05:
            ylimit += ylimit-nrm
        plt.ylim(0,ylimit)
        if i == xlimit-50:
            xlimit += 300
        plt.xlim(0,xlimit)
        plt.scatter(i,nrm)
        plt.draw()


epsilon = 0.05

            
def steepest_descent(liveplot=False):
    if liveplot == True:
        plt.figure(0)
        plt.cla()
        plt.figure(1)
        plt.cla()
    x0  = real_y
    statics = np.array([u1,u2,w0,t0,eps])[np.newaxis]
    
    # These are the initial guesses.  x0 are the true values from which the 
    # fake data were generated
    y0 = x0 + 1*np.array([0.1,0.1,0.1,0.1,0.1])[np.newaxis]
    y = y0
    inp = y0
    z = np.array([0,0,0,0,0])[np.newaxis]
    
    # Step size
    eta = .00005
    
    # Initial gradient
    Du = gradient(y)

    
    # Dummy counting variable
    i = 0
    nrm0 = 10
    resid0 = np.linalg.norm(data-get_curve(np.array(y)))
    
    ylimit = 1.5*resid0
    xlimit = 700
    while np.linalg.norm(data-get_curve(np.array(y))) > epsilon:
        nrm = np.linalg.norm(data-get_curve(np.array(y)))
        
        resid = np.linalg.norm(data-get_curve(np.array(y)))
        
        if resid > resid0:
            eta = eta/2
        elif resid <= resid0:
            eta = eta*1.1
        
        if eta < 1e-5:
            eta = 1e-5
        
        
        resid0 = resid
        
        Du = gradient(y)
        z = y
        y = y - np.matrix(eta*Du)
        i+=1
        if liveplot == True:
            
            # figure(0) has the data and fit
            plt.figure(0)
            plt.cla()
            plt.plot(t,get_curve(np.array(x0)),color="white")
            plt.plot(t,get_curve(np.array(y)),color='red')
            plt.scatter(t,data)
            plt.ylim((0.975, 1.005))
            y = np.array(y)
            
            if np.isnan(np.sum(y)):
                nrm = nrm0
                y = y0
                print 'nan values'
            else:
                nrm0 = nrm
                y0 = y
            textstr = 'iteration %.0f e=%.6f aRs=%.6f i=%.6f\n p0 = %.6f period = %.6f' %(i,
                                                                        y[0][0],
                                                                        y[0][1],
                                                                        y[0][2],
                                                                        y[0][3],
                                                                        y[0][4])
            plt.title(textstr)
            plt.draw()
            
            
            # figure(1) has the norm of the residuals
            plt.figure(1)
            plt.title('Norm of residuals: %.6f'%nrm)
            if nrm > ylimit-0.05:
                ylimit += ylimit-nrm
            plt.ylim(0,ylimit)
            if i == xlimit-50:
                xlimit += 1
            plt.xlim(0,xlimit)
            plt.scatter(i,nrm)
            
            plt.draw()
            
    print 'input: ', inp.T, '\nfitted values: ', y.T, '\nreal values: ',real_y.T,'\n',i,'iterations'
    return y,np.linalg.norm(y-z)

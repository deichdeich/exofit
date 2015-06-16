import numpy as np
import matplotlib.pyplot as plt
import TransitMaker as tm
import math
import time


plt.ion()
plt.figure()

"""
An attempt at Levenberg-Marquadt fitting of transit data
pg 312 in Franklin
"""

### DATA ###
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
i = 1.5 #inclination angle (radians)					
p0 = 0.15 #planet to star radius ratio									
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
on p.312 in Franklin.
"""
t = np.linspace(t0-period/20., t0 + period/20., 1000)

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

def u(curve_x):
    diffs = np.zeros(1000)
    for i in xrange(1000):
        diffs[i] = (data[i]-curve_x[i])**2
    total = np.sum(diffs)
    return total

epsilon = 0.000005
delx = 0.000005

# Calculate the derivative at the point y
# with respect to the variable wrt.  wrt
# is the index of that component in y.
def derivative(a,wrt):
    curve0 = get_curve(np.array(a))
    newa = np.array(a)
    newa[0][wrt] = newa[0][wrt] + delx
    newcurve = get_curve(np.array(newa))
    difference = curve0-newcurve
    dfda = difference/delx
    return dfda

### LEVENBERG-MARQUARDT ###
x0 = real_y #np.array([e,aRs,i,p0,period])[np.newaxis] #make sure that get_lc takes the input in this order
statics = np.array([u1,u2,w0,t0,eps])[np.newaxis]
y0 = x0 + np.array([1,1,1,1,1])[np.newaxis]
y = y0
z = np.array([0,0,0,0,0])[np.newaxis]
lamb0 = 500 #I have no idea how to choose this
lamb = lamb0
lambx = 2
i=0
while np.linalg.norm(y-z) >= epsilon:
    Je = derivative(y,0)
    JaRs = derivative(y,1)
    Ji = derivative(y,2)
    Jp0 = derivative(y,3)
    Jperiod = derivative(y,4)
    Jt = np.array([Je,JaRs,Ji,Jp0,Jperiod])
    Jt = np.matrix(Jt)
        
    
    J = Jt.T
    B = Jt*J
    w = y
    curve_y = get_curve(np.array(y))
    curve_w = get_curve(np.array(w))
    u_y = u(curve_y)
    u_w = u(curve_w)
    if u_y > u_w:
        print 'wrong'
    n = 0
    while u_y <= u_w:
        A = B + lamb * np.diag(B)
        A = np.matrix(A)
        dmf = data-curve_y[np.newaxis]
        dmf = dmf.T
        dx = np.linalg.inv(A)*Jt*dmf
        w = y + dx.T
        curve_w = get_curve(np.array(w))
        u_w = u(curve_w)
        
        if u_y <= u_w:
            lamb = lamb*lambx
        else:
            lamb = lamb/lambx
        
        n += 1
        print 'n: ',n
        
        
    plt.cla()
    plt.plot(t,curve_w,color='red')
    plt.scatter(t,data)
    plt.ylim((0.985, 1.005))
    plt.title(n)
    plt.draw()
        
    i+=1
    print 'i: ',i   
    z = y
    y = w

    
    


print 'initial: ',np.matrix(y0).T
print '\nfitted: ',np.matrix(y).T
print '\nreal: ',np.matrix(real_y).T
plt.ioff()
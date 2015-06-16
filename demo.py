from ext_func.rsky import rsky
from ext_func.occultquad import occultquad
import numpy as np
import math
import matplotlib.pyplot as plt
import cartesian as cart

"""
Calculates a transit light curve assuming a
quadratic LD law:

I(mu) = 1 - u1*(1-mu) - u2*(1-mu)**2

Input:
e 	= eccentricity
i	= inclination angle (radians!)
u1 	= limb darkening parameter 1
u2 	= limb darkening parameter 2
p0 	= planet-to-star radius ratio
w	=
period	= planet orbital period (days)
t0	= ephemeris (BJD)
eps	= error tolerance for calculating the eccentric anomaly
t	= times (BJD) for which to calculate relative flux

Output:
mu_c 	= relative flux at times t
"""
def get_lc(e, aRs, i, u1, u2, p0, w, period, t0, eps,t):
        r_s = 1.0
        npoints = len(t)
        
        #calculates separation of centers between the planet and the star
        z0 = rsky(e, aRs, i, r_s, w, period, t0, eps, t)
        
        #returns limb darkened model lightcurve                
        mu_c = occultquad(z0, u1, u2, p0, npoints)   			
        return mu_c

#eccentricity
e = .0

#semi-major axis to stellar radius ratio	
aRs = 15.23								

#inclination angle (radians)
i = 1.555						

#quadratic limb darkening parameter 1		
u1 = 0.1								

#quadratic limb darkening parameter 2
u2 = 0.3								

#planet to star radius ratio
p0 = 0.115								

#argument of periapse
w = math.pi/2.							

#planet orbital period (days)
period = 0.58040464894       

#ephemeris (JD) 			
t0 = 2454966.52507

#minimum eccentricity for solving Kepler's eqn (otherwise circular orbit assumed)           			
eps = 1.0e-7						



	
t = np.linspace(t0-period/20., t0 + period/20., 1000)

lc = get_lc(e, aRs, i, u1, u2, p0, w, period, t0, eps, t)

"""
fakedata = lc + (.001*np.random.randn(len(lc)))


def get_rsquared(data_point,fit_point):
    datax,datay = data_point[0],data_point[1]
    fitx,fity = fit_point[0],fit_point[1]
    
    rsquared = (datax-fitx)**2+(datay-fity)**2
    
    return rsquared

erange = np.linspace(0, .5, num = 10)
aRsrange = np.linspace(5, 48, num = 10)
irange = np.linspace(80, 95, num = 10)
p0range = np.linspace(0, 4, num = 10)
period = np.linspace(0,55,num = 10)

params = [erange,aRsrange,irange,p0range,period]
perms = cart.cartesian(params)

print "generated permutations"

curves = []

for i in xrange(len(perms)):
    e,aRs,i,p0,period = perms[i]
    lc = get_lc(e,aRs,i,u1,u2,p0,w,period,t0,eps,t)
    curves.append(lc)
    if i%10000 == 0:
        print "generated curve ",i


list1,list2,list3 = [],[],[]


speed up tip:
you can remove the infinities in the loop below.


for curve in curves:
    list1 = []
    for i in xrange(0,1000):
        rsquared = get_rsquared([fakedata[i],i],[curve[i],i])
        list1.append(rsquared)
    list2.append(list1)
    if len(list2)%10000 == 0:
        print "generated r squared values of curve ", len(list2)

print "removing infinities..."
for item in list2:
    for thing in item:
        if thing > 1000000000000000000000000000000:
            item.remove(thing)


for item in list2:
   list3.append(np.nanmean(item))

print "looking for minimum r squared..."
floatMinRSquared = np.nanmin(list3)
indexBestFitCurve = list3.index(floatMinRSquared)

print "best fit curve: ",indexBestFitCurve

plt.plot(curves[indexBestFitCurve],color="red")
lc = get_lc(e, aRs, i, u1, u2, p0, w, period, t0, eps, t)
plt.plot(lc)
"""

plt.plot(t, lc, color='red')
#plt.scatter(t,fakedata)

plt.xlabel("Time (days)")
plt.ylabel("Relative flux")
plt.ylim((0.984, 1.005))
plt.show()
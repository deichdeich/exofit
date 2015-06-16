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


plt.plot(t, lc, color='red')

plt.xlabel("Time (days)")
plt.ylabel("Relative flux")
plt.ylim((0.984, 1.005))
plt.show()
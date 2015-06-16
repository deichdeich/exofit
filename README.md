# exoplanets
An implementation of a steepest descent algorithm for exoplanet transit data.

The only thing here that I wrote is 'ad-exofit.py', which implements a version of the algorithm mentioned above.
I wrote this mostly out of curiosity; it turns out that this method is a lot slower than the Markov Chain Monte
Carlo that's all the rage these days.  There's a neat little feature in steepest-descent.py that lets you watch the
fit and the norm of the residuals happen in real time.  Definition of exhilirating.

The rest of the code (i.e. most of it) is by Laura Kreidberg (http://astro.uchicago.edu/~kreidberg/).

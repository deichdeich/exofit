# exoplanets
An implementation of a steepest descent algorithm for fitting exoplanet transit data to Mandel and Agol light curves
(http://arxiv.org/abs/astro-ph/0210099).

The only thing here that I wrote is 'ad-exofit.py', which implements a version of the algorithm mentioned above.
I wrote this mostly out of curiosity; it turns out that this method is a lot slower than the Markov Chain Monte
Carlo that's all the rage these days.  It also gets stuck in local minima.  The solution to this is straight forward
-- if there has been no change in the fit for some number of steps, increase the dx until there is change -- I just
haven't gotten around to coding this.

There's a neat little feature that lets you watch the fit and the norm of the residuals happen in real time.
Definition of exhilirating.

The rest of the code (i.e. most of it) is by Laura Kreidberg (http://astro.uchicago.edu/~kreidberg/).

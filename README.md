# exoplanets
An implementation of a steepest descent algorithm for exoplanet transit data.

The only thing here that I wrote is 'steepest-descent.py', which implements a variety of the algorithm
that is its namesake.  I wrote this mostly out of curiosity; it turns out that this method is a lot
slower than the Markov Chain Monte Carlo that's all the rage these days.  There's a neat little feature
in steepest-descent.py that lets you watch the norm of the residual slowly drop as the fit hones in
on the answer.

The rest of the code (i.e. most of it) is by Laura Kreidberg (http://astro.uchicago.edu/~kreidberg/).

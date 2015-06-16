import numpy as np
from math import pi, log
"""
peakfinder.py:  finds peaks and troughs (max and min) in spectra data

Author: Alex Deich

inputs-
y_axis: list-like record of the spectra data
x_axis (optional): list-like record of the x-axis.  Could be wavenumber, whatever.
If not provided, x_axis is taken to be the index of y_axis.
howfar: parameter which indicates how far ahead the script should look for the
next peak.  Too low is inefficient, too high will miss peaks.
delta: parameter which indicates an average peak width.
returns-
[spect_peaks,spect_mins], 2 2x2 lists of x,y points corresponding to
the peaks and troughs.
"""
#checks sanity of input data
def datasanitize_peakfinder(x_axis, y_axis):
    if x_axis is None:
        x_axis = range(len(y_axis))
    if len(y_axis) != len(x_axis):
        raise (ValueError,
                ’Input vectors y_axis and x_axis must have same length’)
    y_axis = np.array(y_axis)
    x_axis = np.array(x_axis)
    return x_axis, y_axis
def peakfinder(y_axis, x_axis = None, howfar = 2, delta=0):
    spect_peaks = []
    spect_mins = []
    toss = []
    x_axis, y_axis = _datasanitize_peakfinder(x_axis, y_axis)
    # store data length for later use

length = len(y_axis)
if howfar < 1:
    raise ValueError, "howfar must be ’1’ or above in value"
if not (np.isscalar(delta) and delta >= 0):
    raise ValueError, "delta must be a positive number"
  mn, mx = np.Inf, -np.Inf
for index, (x, y) in enumerate(zip(x_axis[:-howfar],
                                    y_axis[:-howfar])):
if y > mx: mx = y
        mxpos = x
    if y < mn:
mn = y mnpos = x
    if y < mx-delta and mx != np.Inf:
        if y_axis[index:index+howfar].max() < mx:
            spect_peaks.append([mxpos, mx])
            toss.append(True)
            mx = np.Inf
            mn = np.Inf
            if index+howfar >= length:
                break
            continue
    if y > mn+delta and mn != -np.Inf:
try:
if y_axis[index:index+howfar].min() > mn:
    spect_mins.append([mnpos, mn])
    toss.append(False)
    mn = -np.Inf
    mx = -np.Inf
    if index+howfar >= length:
break

if toss[0]:
        spect_peaks.pop(0)
    else:
        spect_mins.pop(0)
    del toss
except IndexError:
    pass
return [spect_peaks, spect_mins]
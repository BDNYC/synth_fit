#!/usr/bin/python
# Utilities
import numpy as np


def smooth(x, beta):
    """
    Smooths a spectrum *x* using a Kaiser-Bessel smoothing window of narrowness *beta* (~1 => very smooth, ~100 => not smooth)
    """
    window_len = 11
    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    w = np.kaiser(window_len, beta)
    y = np.convolve(w / w.sum(), s, mode='valid')
    return y[5:len(y) - 5] * (x.unit if hasattr(x, 'unit') else 1)


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


def rebin_spec(spec, wavnew, waveunits='um'):
    from pysynphot import spectrum, observation
    # Gives same error answer: Err = np.array([np.sqrt(sum(spec[2].value[idx_include(wavnew,[((wavnew[0] if n==0 else
    #  wavnew[n-1]+wavnew[n])/2,wavnew[-1] if n==len(wavnew) else (wavnew[n]+wavnew[n+1])/2)])]**2)) for n in range(
    # len(wavnew)-1)])*spec[2].unit if spec[2] is not '' else ''
    if len(spec) == 2:
        spec += ['']
    try:
        Flx, Err, filt = spectrum.ArraySourceSpectrum(wave=spec[0].value,
                                                      flux=spec[1].value), spectrum.ArraySourceSpectrum(
            wave=spec[0].value, flux=spec[2].value) if spec[2] else '', spectrum.ArraySpectralElement(spec[0].value,
                                                                                                      np.ones(
                                                                                                          len(spec[0])),
                                                                                                      waveunits=waveunits)
    except:
        spec, wavnew = [i * q.Unit('') for i in spec], wavnew * q.Unit('')
        Flx, Err, filt = spectrum.ArraySourceSpectrum(wave=spec[0].value,
                                                      flux=spec[1].value), spectrum.ArraySourceSpectrum(
            wave=spec[0].value, flux=spec[2].value) if spec[2] else '', spectrum.ArraySpectralElement(spec[0].value,
                                                                                                      np.ones(
                                                                                                          len(spec[0])),
                                                                                                      waveunits=waveunits)
    return [wavnew, observation.Observation(Flx, filt, binset=wavnew.value, force='taper').binflux * spec[1].unit,
            observation.Observation(Err, filt, binset=wavnew.value, force='taper').binflux * spec[2].unit if spec[
                2] else np.ones(len(wavnew)) * spec[1].unit]
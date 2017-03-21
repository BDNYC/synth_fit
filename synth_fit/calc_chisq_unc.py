# Calculate Chi-Squared for all models in a grid to determine the starting
# point for the emcee walkers
# Stephanie Douglas
################################################################################

import logging

import numpy as np
from astropy import units as u
import pickle
from smooth import *
import matplotlib.colors as colors

from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
from pylab import *
cmap = cm.get_cmap('RdPu')

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

cmap = plt.get_cmap('RdPu')
sub_cmap = LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=0.2, b=1),cmap(np.linspace(0.2, 1, 6)))
new_cmap = discrete_cmap(6, base_cmap=sub_cmap)

def calc_chisq(data_flux,data_unc,model_flux):
    a = (data_flux-model_flux)**2
    b = (data_unc**2)
    cs=[]
    for i in range(len(a)):
    	if float(b[i].value)==0.0:
    		bb==0.000000000000000001
    	else: bb=float(b[i].value)	
    	c = float(a[i].value)/float(bb)
    	cs.append(c)
    d = np.sum(cs)
    return c

def test_all(data_wave, data_flux, data_unc, model_dict, params,
    smooth=False,resolution=None,shortname=''):
    """
    Calculates chi-squared for all models in a grid to determine
    the starting point for the emcee walkers (or just to find the
    model that minimizes chi-squared...)

    NOTE: at this point I have not accounted for model parameters
    that are NOT being used for the fit

    Parameters
    ----------
    data_wave: array; astropy.units Quantity

    data_flux: array; astropy.units Quantity

    data_unc: array; astropy.units Quantity

    model_dict: dictionary
        keys 'wavelength' and 'flux' should correspond to model wavelength and 
        flux arrays, and those should be astropy.units Quantities
        other keys should correspond to params

    params: array of strings
        the model parameters to be interpolated over.  These should 
        correspond to keys of model_dict

    smooth: boolean (default=True)
        whether or not to smooth the model spectra before interpolation 
        onto the data wavelength grid 
        (a check will be performed before interpolation to see if it's
        it's necessary)

    resolution: astropy.units Quantity (optional)
        Resolution of the input DATA, to be used in smoothing the model.
        Only relevant if smooth=True

    """

    #logging.debug(data_wave.unit.to_string('fits'))
    #logging.debug(data_unc.unit.to_string('fits'))
    #logging.debug(model_dict['wavelength'].unit.to_string('fits'))
    if data_wave.unit!=model_dict['wavelength'].unit:
         logging.debug('changing units')
         data_wave = data_wave.to(model_dict['wavelength'].unit)
    #logging.debug(data_wave.unit.to_string('fits'))

    if ((len(model_dict['wavelength'])==len(data_wave)) and 
        (np.sum(model_dict['wavelength']-data_wave)<(model_dict['wavelength'].unit*1e-12))):
        interp = False
        logging.info('calc_chisq.test_all: NO INTERPOLATION')
    else:
        interp = True
        logging.info('calc_chisq.test_all: INTERPOLATION NEEDED')


    ndim = len(params)

    num_models = len(model_dict['flux'])

    chisq = np.ones(num_models)*(99e15)
		
    save_chisq = []
    
    for i in range(num_models):
#        logging.debug('%d %d %f %f',i, num_models, model_dict['logg'][i], 
#            model_dict['teff'][i])
        if smooth:
            mod_flux = falt2(model_dict['wavelength'],model_dict['flux'][i],resolution)
        else:
            mod_flux = model_dict['flux'][i]
            #logging.debug('shape flux {} mf {}'.format(np.shape(model_dict['flux']), np.shape(mod_flux)))

        #logging.debug('lengths dw {} modw {} modf {}'.format(
        #    len(data_wave),len(model_dict['wavelength']),len(mod_flux)))
        if interp:
            mod_flux = np.interp(data_wave,model_dict['wavelength'],mod_flux)
#         mod_flux=mod_flux*u.erg/u.AA/u.cm**2/u.s
#        logging.debug(str(mod_flux[100:110]))
#        logging.debug('stdev %f', np.std(mod_flux))
#     
#         mult1 = data_flux*mod_flux
#         bad = np.isnan(mult1)
#         mult = np.sum(mult1[~bad])
#         sq1 = mod_flux**2
#         square = np.sum(sq1[~bad])
#         ck = mult/square
#         mod_flux = mod_flux*ck
        mult1 = data_flux*mod_flux/(data_unc**2)
        bad = np.isnan(mult1)
        mult = mult1[~bad]
        bad2 = np.isinf(mult)
        mult = np.sum(mult[~bad2])
        sq1 = mod_flux*mod_flux/(data_unc**2)
        mult2 = sq1[~bad]
        
        mult2 = float(sum(mult2[~bad2]))
        ck = mult/mult2
        mod_flux=mod_flux*ck
        
        chisq[i] = calc_chisq(data_flux, data_unc, mod_flux)
        foo = plt.scatter(model_dict['teff'][i],chisq[i],c=model_dict['logg'][i],cmap=new_cmap,edgecolor='None',vmin=2.75,vmax=5.75)
        params_list = []
        for p in params:
         params_list.append(model_dict[p][i])
        save_chisq.append([params_list,chisq[i]])    
    
    min_loc = np.argmin(chisq)
#    logging.debug('min_loc %d', min_loc)
    best_params = np.array([])
    for p in params:
        best_params = np.append(best_params,model_dict[p][min_loc])
        
    plt.xlabel('$T_{eff}$',fontsize='x-large')
    plt.ylabel('Goodness of Fit',fontsize='x-large')
    plt.annotate('{}, {}'.format(shortname,best_params),xy=(0.5,0.95),xycoords='axes fraction')
    fig = plt.gcf()
    labels = np.arange(3.0,6.0,0.5)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    cbar = fig.colorbar(foo, cax=cbar_ax, ticks=labels)
    cbar.set_label(label='$log(g)$',size=18)
    plt.savefig('/Users/paigegiorla/Desktop/Teff_goodnesses_{}'.format(shortname)+'.pdf')
    plt.clf()
   
        
    fb = open('/Users/paigegiorla/Desktop/chisquares.pkl','wb')
    pickle.dump(save_chisq,fb)
    fb.close()
    return best_params,min(chisq)

import astropy.units as q
import pickle
import mcmc_fit as mc
import numpy as np


# Read in SED of 1256-0224
SED_path = '/Users/EileenGonzales/Dropbox/BDNYC/BDNYC_Research/Python/Modules/Atmospheres_paper/Redone/best1256-0224 (L3.5sd) SED.txt'
w, f, e = np.loadtxt(SED_path, delimiter= " ", unpack=True)


# Adding units to the arrays
w=w*q.um
f=f* q.erg/q.AA/q.cm**2/q.s
e=e* q.erg/q.AA/q.cm**2/q.s


# Load the pickle file of the model grid
pickle_path = '/Users/EileenGonzales/Dropbox/BDNYC/BDNYCdb_copy/BTSettl_mcmc.pkl'
file = open(pickle_path,'rb’)  #rb means read binary
models = pickle.load(file)
file.close()


# This makes the model grid
mg=mc.make_model_db('btsettl', 'model_atmosphere_db', model_grid=models, grid_data='spec', param_lims=[('teff',400,1600,50),('logg',3.5,5.5,0.5)], fill_holes=False, bands=[], rebin_models=w,use_pandas=False)

# This does the actual fitting to the spectra
bdsamp=mc.fit_spectrum([w,f,e], mg, 'btsettl', 'GJ758B',25, 50, mask=[], db='', extents=None,object_name='Test', log=False, plot=True, prnt=True, generate=True,outfile=None)
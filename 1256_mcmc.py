import astropy.units as q
import pickle
import mcmc_fit as mc
import numpy as np


# Read in SED of 1256-0224
SED_path = '/Users/EileenGonzales/Dropbox/BDNYC/BDNYC_Research/Python/Modules/Atmospheres_paper/Redone/best1256-0224 (L3.5sd) SED.txt'
w, f, e = np.loadtxt(SED_path, delimiter=" ", unpack=True)


# Adding units to the arrays
w = w*q.um
f = f* q.erg/q.AA/q.cm**2/q.s
e = e* q.erg/q.AA/q.cm**2/q.s


# Load the pickle file of the model grid
pickle_path = '/Users/EileenGonzales/Dropbox/BDNYC/BDNYCdb_copy/BTSettl_mcmc.pkl'
file = open(pickle_path, 'rb')  # rb means read binary
models = pickle.load(file)
file.close()


# This makes the model grid
mg=mc.make_model_db('btsettl', 'model_atmosphere_db', model_grid=models, grid_data='spec', param_lims=[('teff',1400,2000,50),('logg',3.5,5.5,0.5)], fill_holes=False, bands=[], rebin_models=w,use_pandas=False)

# I can't get the model grid to be made error:

# Traceback (most recent call last):
#   File "/Users/EileenGonzales/anaconda/lib/python2.7/site-packages/IPython/core/interactiveshell.py", line 3066, in run_code
#     exec(code_obj, self.user_global_ns, self.user_ns)
#   File "<ipython-input-55-e03a4670e4f1>", line 1, in <module>
#     mg=mc.make_model_db('btsettl', 'model_atmosphere_db', model_grid=models, grid_data='spec', param_lims=[('teff',400,1600,50),('logg',3.5,5.5,0.5)], fill_holes=False, bands=[], rebin_models=w,use_pandas=False)
#   File "/Users/EileenGonzales/Dropbox/BDNYC/BDNYC_Research/Python/Modules/synth_fit/mcmc_fit.py", line 116, in make_model_db
#     models = pd.DataFrame(model_grid)
#   File "/Users/EileenGonzales/anaconda/lib/python2.7/site-packages/pandas/core/frame.py", line 226, in __init__
#     mgr = self._init_dict(data, index, columns, dtype=dtype)
#   File "/Users/EileenGonzales/anaconda/lib/python2.7/site-packages/pandas/core/frame.py", line 363, in _init_dict
#     dtype=dtype)
#   File "/Users/EileenGonzales/anaconda/lib/python2.7/site-packages/pandas/core/frame.py", line 5158, in _arrays_to_mgr
#     index = extract_index(arrays)
#   File "/Users/EileenGonzales/anaconda/lib/python2.7/site-packages/pandas/core/frame.py", line 5206, in extract_index
#     raise ValueError('arrays must all be same length')
# ValueError: arrays must all be same length



# This does the actual fitting to the spectra
bdsamp=mc.fit_spectrum([w,f,e], mg, 'btsettl', 'GJ758B',25, 50, mask=[], db='', extents=None,object_name='Test', log=False, plot=True, prnt=True, generate=True,outfile=None)

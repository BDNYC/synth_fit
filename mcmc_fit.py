import logging
import cPickle
import astropy.units as u
import BDdb
import synth_fit
import synth_fit.bdfit
import numpy as np

def make_model_db(model_grid_name):
	'''
	Given a **model_grid_name**, returns the grid from the model_atmospheres.db in the proper format to work with fit_spectrum()
	'''
	# Load the model_atmospheres database and pull all the data from the specified table
	db = BDdb.get_db('/Users/paigegiorla/Code/Python/BDNYC/model_atmospheres.db')
	mg = db.dict.execute("Select * from {}".format(model_grid_name)).fetchall()
	
	# Make a new dictionary with each parameter as an array of values
	mg = {'id':np.array([row['id'] for row in mg]),
	'teff':np.array([row['teff'] for row in mg]),
	'logg':np.array([row['logg'] for row in mg]),
	'wavelength':np.array([np.array(row['wavelength'],dtype=np.float64) for row in mg]),
	'flux':np.array([np.array(row['flux'],dtype=np.float64) for row in mg])}	
	
	# Apply astropy.units to wavelength and flux arrays
	mg['flux'] = (u.erg/u.AA/u.cm**2/u.s)*mg['flux']
	mg['wavelength'] = (u.um)*mg['wavelength']
	
	return mg

def fit_spectrum(raw_spectrum, model_grid_name, object_name='Test', log=False):
	'''
	Given **raw_spectrum** as an integer id from the SPECTRUM table or a [W,F,E] list with astropy units, 
	returns a marginalized distribution plot of best fit parameters from the specified **model_grid_name**.
	'''
	
	if log: logging.basicConfig(level=logging.DEBUG)
	
	# Turn the model_atmospheres.db grid into one that can talk to Steph's code
	mg = make_model_db(model_grid_name)
	
	# Input can be [W,F,E] or an id from the SPECTRUM table of the BDNYC Data Archive
	if isinstance(raw_spectrum,(float,int)):
		db = BDdb.get_db('/Users/paigegiorla/Desktop/PG_DB_2_16_15.db')
		query_spectrum = db.dict.execute("SELECT * from spectra where id={}".format(raw_spectrum)).fetchone()
	
		# turn into a dictionary with astropy units quantities
		wave_unit = u.Unit(query_spectrum['wavelength_units'])
		flux_unit = u.Unit(query_spectrum['flux_units'].replace("ergss","erg s"
		).replace('ergs','erg s'
		).replace('Wm','W m'
		).replace("A","AA"))
		
		spectrum = {'wavelength': query_spectrum['wavelength']*wave_unit,
		'flux': query_spectrum['flux']*flux_unit,
		'unc': query_spectrum['unc']*flux_unit}
	
	else: 
		spectrum={'wavelength':raw_spectrum[0],'flux':raw_spectrum[1],'unc':raw_spectrum[2]}
	
	if log: logging.debug(model_grid['wavelength'].shape); logging.debug(model_grid['wavelength'])
	
	# Specify the parameter space to be walked
	params = [i for i in model_grid.keys() if i in ['logg', 'teff', 'f_sed', 'k_zz']]
	
	# Set up the sampler object (it's a wrapper around emcee)
	bdsamp = synth_fit.bdfit.BDSampler(object_name, spectrum, model_grid,	params, smooth=False,	plot_title=plot_title, snap=True) # smooth=False if model already matches data, snap=True if no interpolation is needed on grid
																				        
	# Run the mcmc method
	bdsamp.mcmc_go(nwalk_mult=2,nstep_mult=10)
	
	# Plotting
	plot_title = "{}, {}".format(object_name,"BT-Settl 2013")
	bdsamp.plot_triangle()
	bdsamp.plot_chains()
	
	if log: logging.info("ran MCMC"); logging.info("all done!")

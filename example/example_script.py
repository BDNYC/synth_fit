## Example script for synth_fit

import logging
import cPickle
import astropy.units as u
import BDdb
import synth_fit
import synth_fit.bdfit
import numpy as np

mg=BDdb.get_db('/Users/paigegiorla/Code/Python/BDNYC/model_atmospheres.db').dict.execute("Select * from bt_settl_2013").fetchall()
mg={'id':np.array([i[0] for i in mg]),
	'teff':np.array([i[1] for i in mg]),
	'logg':np.array([i[2] for i in mg]),
	'wavelength':np.array([np.array(i[4],dtype=np.float64) for i in mg]),
	'flux':np.array([np.array(i[5],dtype=np.float64) for i in mg])}	

def fit_spectrum(raw_spectrum,model_grid=mg,object_name=''):
	'''raw_spectrum requires units! Must be astropy quantities.
	'''

	logging.basicConfig(level=logging.INFO)

	# load the database - replace with appropriate path
	
	# object_name = '0036+1821'
	# object_name = 86

	if isinstance(raw_spectrum,(float,int)):
		db=BDdb.get_db('/Users/paigegiorla/Desktop/PG_DB_2_16_15.db')
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
	

	# now open up the model file
	# infile = open("SpeX_marley_nolowg.pkl","rb")
	# model = cPickle.load(infile)
	# infile.close()

	# change into astropy units quantities
	model_grid['flux'] = (u.erg/u.AA/u.cm**2/u.s)*model_grid['flux']
	model_grid['wavelength'] = (u.um)*model_grid['wavelength']
	params = [i for i in model_grid.keys() if i in ['logg', 'teff', 'f_sed', 'k_zz']]

	# now set up the sampler object (it's a wrapper around emcee)

	plot_title = "TESTING {}, {}".format(object_name,"BT-Settl 2013")

	bdsamp = synth_fit.bdfit.BDSampler(object_name, 
									   spectrum,
									   model_grid,
									   params,
									   smooth=False, # model already matches data
									   plot_title=plot_title,
									   snap=True) # no interpolation on grid

	#this isn't enough for real results, just to make sure it runs
	bdsamp.mcmc_go(nwalk_mult=2,nstep_mult=10)

	bdsamp.plot_triangle()
# 	bdsamp.plot_chains()

	logging.info("ran MCMC")

	logging.info("all done!")

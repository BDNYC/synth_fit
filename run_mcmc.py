## START HERE FOR MCMC on one object

import logging
import cPickle
import astropy.units as u
import BDdb
import synth_fit
import synth_fit.bdfit
import numpy as np
from matplotlib import pyplot as plt

db=BDdb.get_db('/Users/paigegiorla/Code/Python/BDNYC/model_atmospheres.db')
mg=db.dict.execute("Select * from bt_settl_2013 where teff between 900 and 1700").fetchall()
mg={'id':np.array([row['id'] for row in mg]),
	'teff':np.array([row['teff'] for row in mg]),
	'logg':np.array([row['logg'] for row in mg]),
	'wavelength':np.array([np.array(row['wavelength'],dtype=np.float64) for row in mg]),
	'flux':np.array([np.array(row['flux'],dtype=np.float64) for row in mg])}	
# change into astropy units quantities
mg['flux'] = (u.erg/u.AA/u.cm**2/u.s)*mg['flux']
mg['wavelength'] = (u.um)*mg['wavelength']
print mg['flux']

def fit_spectrum(raw_spectrum,object_name,model_grid=mg):
	'''raw_spectrum requires units! Must be astropy quantities.
	'''

	logging.basicConfig(level=logging.DEBUG)

	# load the database - replace with appropriate path
	
	# object_name = '0036+1821'
# 	object_name = 1580

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
	
	logging.debug(model_grid['wavelength'].shape)
	logging.debug(model_grid['wavelength'])

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
 	plt.show()
 	bdsamp.plot_chains()
	plt.show()
	logging.info("ran MCMC")

	logging.info("all done!")

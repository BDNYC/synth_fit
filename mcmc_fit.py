import logging, BDdb, cPickle, synth_fit, astropy.units as q, utilities as u, numpy as np, matplotlib.pyplot as plt

def model_grid_interp_test(mg, teff, logg):
  idx1, idx2, idx3 = zip(mg['teff'],mg['logg']).index((teff[0]-teff[1],logg[0]-logg[1])), zip(mg['teff'],mg['logg']).index((teff[0],logg[0])), zip(mg['teff'],mg['logg']).index((teff[0]+teff[1],logg[0]+logg[1]))
  w1, f1 = mg['wavelength'][idx1], mg['flux'][idx1]
  w2, f2 = mg['wavelength'][idx2], mg['flux'][idx2]
  w3, f3 = mg['wavelength'][idx3], mg['flux'][idx3]
  
  W, F = interp_models(['teff','logg'], [teff[0],logg[0]], mg)
  
  plt.loglog(W, F, color='k', lw=2, label='Interpolated')
  
  for x in [[w1,f1,'r',(teff[0]-teff[1],logg[0]-logg[1])],[w2,f2,'g',(teff[0],logg[0])],[w3,f3,'b',(teff[0]+teff[1],logg[0]+logg[1])]]:
    f = u.smooth(x[1], 1)
    plt.loglog(x[0].value, f.value, color=x[2], label='{} {}'.format(*x[3]), alpha=0.5)
  plt.legend(loc=1)
  
def make_model_db(model_grid_name, model_atmosphere_db, param_lims=[('teff',400,2500)], rebin_models=True, complete_grid=False):
	'''
	Given a **model_grid_name**, returns the grid from the model_atmospheres.db in the proper format to work with fit_spectrum()
	'''
	# Load the model_atmospheres database and pull all the data from the specified table
	db = BDdb.get_db(model_atmosphere_db)
	if param_lims:
		limit_text = ' AND '.join(['{} between {} and {}'.format(i,j,k) for i,j,k in param_lims])
		model_grid = db.dict.execute("Select * from {} where {}".format(model_grid_name,limit_text)).fetchall()  
	else: model_grid = db.dict.execute("Select * from {}".format(model_grid_name)).fetchall()
	
	# Make a new dictionary with each parameter as an array of values
	mg, wavelength, flux = {k:np.array([row.get(k) for row in model_grid]) for k in ['teff','logg','k_zz','f_sed']}, [], []
	
  # Choose template wavelength array to rebin all other spectra
	W = np.array(model_grid[0]['wavelength'], dtype=np.float64)
	
  # Pull out model spectra and rebin if necessary
	for row in model_grid:
	  w, f, e = u.rebin_spec([np.array(row['wavelength'], dtype=np.float64)*q.um, np.array(row['flux'], dtype=np.float64)*q.erg/q.s/q.cm**2/q.AA], W*q.um) if rebin_models else [np.array(row['wavelength'], dtype=np.float64), np.array(row['flux'], dtype=np.float64), None]
	  wavelength.append(w), flux.append(f)
	
  # Add astropy units and put spectra in the model grid
	mg.update({'wavelength':(q.um)*np.array(wavelength), 'flux':(q.erg/q.AA/q.cm**2/q.s)*np.array(flux)})
	
	if complete_grid:
    # Put Paige's code in here!
	  
    # Pull out the missing values in each parameter space
    # Run interp_model() at those parameters to generate interpolated spectrum
    # Insert the new spectra into the model grid
	  pass
	
	return mg

def fit_spectrum(raw_spectrum, model_grid, walkers, steps, object_name='Test', log=False, plot=True, prnt=True, outfile=None):
	'''
	Given **raw_spectrum** as an integer id from the SPECTRUM table or a [W,F,E] list with astropy units, 
	returns a marginalized distribution plot of best fit parameters from the specified **model_grid_name**.
	'''
	
	if log: logging.basicConfig(level=logging.DEBUG)
		
	# Input can be [W,F,E] or an id from the SPECTRUM table of the BDNYC Data Archive
	if isinstance(raw_spectrum,(float,int)):
		db = BDdb.get_db('/Users/paigegiorla/Desktop/PG_DB_2_16_15.db')
		query_spectrum = db.dict.execute("SELECT * from spectra where id={}".format(raw_spectrum)).fetchone()
	
		# Turn into a dictionary with astropy units quantities
		wave_unit = u.Unit(query_spectrum['wavelength_units'])
		flux_unit = u.Unit(query_spectrum['flux_units'].replace("ergss","erg s").replace('ergs','erg s').replace('Wm','W m').replace("A","AA"))
		
		spectrum = {'wavelength':query_spectrum['wavelength']*wave_unit, 'flux':query_spectrum['flux']*flux_unit,	'unc':query_spectrum['unc']*flux_unit}
	
	else: 
		spectrum={'wavelength':raw_spectrum[0],'flux':raw_spectrum[1], 'unc':raw_spectrum[2]}
	
	if log: logging.debug(model_grid['wavelength'].shape); logging.debug(model_grid['wavelength'])
	
	# Specify the parameter space to be walked
	params = [i for i in model_grid.keys() if i in ['logg', 'teff', 'f_sed', 'k_zz']]
	
	# Set up the sampler object (it's a wrapper around emcee)
	bdsamp = synth_fit.bdfit.BDSampler(object_name, spectrum, model_grid,	params, smooth=False,	plot_title="{}, {}".format(object_name,"BT-Settl 2013"), snap=True) # smooth=False if model already matches data, snap=True if no interpolation is needed on grid
																				        
	# Run the mcmc method
	bdsamp.mcmc_go(nwalk_mult=walkers,nstep_mult=steps,outfile=outfile)
	
	# Plotting
	if plot:
	  bdsamp.plot_triangle()
	  bdsamp.plot_chains()
	
  # Printing
	if log: logging.info("ran MCMC"); logging.info("all done!")
	
	if prnt:
		print bdsamp.all_params
		print bdsamp.all_quantiles.T[1]
	
  # Generate best fit spectrum with a list of the parameter, 50th quantile value, and uncertainty.
	bdsamp.best_fit_spectrum = interp_models(bdsamp.all_params, bdsamp.all_quantiles.T[1], model_grid)
	
	return bdsamp
	
def interp_models(params, coordinates, model_grid, smoothing=1):
  """
  Interpolation code that accepts a model grid and a list of parameters/values to return an interpolated spectrum.

  Parameters
  ----------
  params: list
      A list of the model parameters, e.g. ['teff', 'logg', 'f_sed']
  coordinates: list
      A list of the coordinates in parameter space to evaluate, e.g. [1643, 5.1, 2.3]
  model_grid: object
      The output of make_model_db()
  
  Returns
  -------
  spectrum: list of arrays
       The wavelength and flux at the specified **values** in parameter space

  Notes
  -----
  You might have to update scipy and some other things to run this. Updating takes about 1.5 hours! Do:
  >>> ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
  >>> brew install gcc
  >>> pip install scipy --upgrade

  """
  from scipy.interpolate import LinearNDInterpolator
  
  # Transpose the flux arrays to interpolate the grid at each wavelength position
  flux_columns = [np.array([i[n] for i in model_grid['flux'].value]) for n in range(len(model_grid['flux'].value[0]))]
  
  # Build parameter space from arrays
  grid = [model_grid[p] for p in params]
  
  # Define the wavelength array and an empty flux array
  W, F = model_grid['wavelength'][0].value, np.zeros(len(flux_columns))
  
  # Interpolate to specified coordinates for each wavelength point in parameter space
  for n in range(len(flux_columns)):
    
    # Create grid interpolation function to pass coordinates to
    interpND = LinearNDInterpolator(zip(*grid), flux_columns[n], rescale=True)
  
    # Find flux value at desired coordinates in parameter space and insert into interpolated flux array
    F[n] = interpND(coordinates)
  
  return [W,u.smooth(F,smoothing) if smoothing else F]

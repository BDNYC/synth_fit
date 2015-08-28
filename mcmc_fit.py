import logging, BDdb, cPickle, synth_fit, synth_fit.bdfit, astropy.units as q, utilities as u, numpy as np, matplotlib.pyplot as plt, pandas as pd

def fit_spectrum(raw_spectrum, model_grid, walkers, steps, object_name='Test', log=False, plot=True, prnt=True, outfile=None):
	'''
	Given **raw_spectrum** as an integer id from the SPECTRUM table or a [W,F,E] list with astropy units, 
	returns a marginalized distribution plot of best fit parameters from the specified **model_grid** name.
	
	Parameters
	----------
	raw_spectrum: int, dict
	  An integer id for the desired spectrum from the SPECTRA table or a dictionary with the wavelength and flux arrays to be fit
	model_grid: str
	  The name of the model grid to be used in the fit, e.g. 'bt_settl_2013'
	walkers: int
	  The number of walkers to deploy
	steps: int
	  The number of steps for each walker to take
	
	Retiurns
	--------
	bdsamp: object
	  The MCMC result instance
	'''
	
	if log: logging.basicConfig(level=logging.DEBUG)
		
	# Input can be [W,F,E] or an id from the SPECTRUM table of the BDNYC Data Archive
	if isinstance(raw_spectrum,int):
		db = BDdb.get_db('/Users/paigegiorla/Desktop/PG_DB_2_16_15.db')
		query_spectrum = db.dict("SELECT * from spectra where id={}".format(raw_spectrum)).fetchone()
	
		# Turn into a dictionary with astropy units quantities
		wave_unit = u.Unit(query_spectrum['wavelength_units'])
		flux_unit = u.Unit(query_spectrum['flux_units'].replace("ergss","erg s").replace('ergs','erg s').replace('Wm','W m').replace("A","AA"))
		
		spectrum = {'wavelength':query_spectrum['wavelength']*wave_unit, 'flux':query_spectrum['flux']*flux_unit,	'unc':query_spectrum['unc']*flux_unit}
	
	else: 
		spectrum = {'wavelength':raw_spectrum[0], 'flux':raw_spectrum[1], 'unc':raw_spectrum[2]}
	
	if log: logging.debug(model_grid['wavelength'].shape); logging.debug(model_grid['wavelength'])
	
	# Specify the parameter space to be walked
	params = [i for i in model_grid.keys() if i in ['logg', 'teff', 'f_sed', 'k_zz']]
	
	# Set up the sampler object (it's a wrapper around emcee)
	bdsamp = synth_fit.bdfit.BDSampler(object_name, spectrum, model_grid,	params, smooth=False,	plot_title="{}, {}".format(object_name,"BT-Settl 2013"), snap=True) # smooth=False if model already matches data, snap=True if no interpolation is needed on grid
																				        
	# Run the mcmc method
	bdsamp.mcmc_go(nwalk_mult=walkers, nstep_mult=steps, outfile=outfile)
	
	# Plotting
	if plot:
	  bdsamp.plot_triangle()
	  bdsamp.plot_chains()
	
  # Printing
	if log: logging.info("ran MCMC"); logging.info("all done!")
	
	if prnt:
		print bdsamp.all_params
		print bdsamp.all_quantiles.T[1]
	
  # Generate best fit spectrum the 50th quantile value
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
  
  # Take out nusiance parameters and build parameter space from arrays
  params, coordinates = zip(*[i for i in zip(params,coordinates) if i[0] in ['teff','logg','k_zz','f_sed']])
  grid = [model_grid.get(p) for p in params]
  
  # Define the wavelength array and an empty flux array
  W, F = model_grid['wavelength'][0].value, np.zeros(len(flux_columns))
  
  # Interpolate to specified coordinates for each wavelength point in parameter space
  for n in range(len(flux_columns)):
    
    # Create grid interpolation function to pass coordinates to
    interpND = LinearNDInterpolator(zip(*grid), flux_columns[n], rescale=True)
  
    # Find flux value at desired coordinates in parameter space and insert into interpolated flux array
    F[n] = interpND(coordinates)
  
  return [W,u.smooth(F,smoothing) if smoothing else F]

def make_model_db(model_grid_name, model_atmosphere_db, param_lims=[('teff',400,1200)], rebin_models=True, grid_resolution=[('teff',50),('logg',0.5)]):
  '''
  Given a **model_grid_name**, returns the grid from the model_atmospheres.db in the proper format to work with fit_spectrum()
  '''
  # Load the model_atmospheres database and pull all the data from the specified table
  db = BDdb.get_db(model_atmosphere_db)
  if param_lims:
    limit_text = ' AND '.join(['{} between {} and {}'.format(i,j,k) for i,j,k in param_lims])
    model_grid = db.dict("Select * from {} where {}".format(model_grid_name,limit_text)).fetchall()  
  else: model_grid = db.dict("Select * from {}".format(model_grid_name)).fetchall()
  
  # Make a new dictionary with each parameter as an array of values
  mg, wavelength, flux = {k:np.array([row.get(k) for row in model_grid]) for k in ['teff','logg','k_zz','f_sed'] if k in model_grid[0]}, [], []
  
  # Choose template wavelength array to rebin all other spectra
  W = np.array(model_grid[0]['wavelength'], dtype=np.float64)
  
  # Pull out model spectra and rebin if necessary
  for row in model_grid:
    w, f, e = u.rebin_spec([np.array(row['wavelength'], dtype=np.float64)*q.um, np.array(row['flux'], dtype=np.float64)*q.erg/q.s/q.cm**2/q.AA], W*q.um) if rebin_models else [np.array(row['wavelength'], dtype=np.float64), np.array(row['flux'], dtype=np.float64), None]
    wavelength.append(w), flux.append(f)
    
  for k in grid_resolution:
    
    # Define new axis for this parameter at the given resolution using the min and max of the existing grid as end points
    # mg[k[0]+'_new'] = np.arange(min(mg[k[0]]),max(mg[k[0]]),k[1])
    pass
    
    
    # No no no. What I want to do is fill in a grid. For each missing point, I only want to fill it in if there is one above and below OR one left and right.
    
    #     o o o o o 
    # o o x o x o o
    # o x x o o
    #   o o o o
  
  # Add astropy units and put spectra in the model grid
  mg.update({'wavelength':(q.um)*np.array(wavelength), 'flux':(q.erg/q.AA/q.cm**2/q.s)*np.array(flux)})
  
  return mg
	
# def make_model_db(model_grid_name, model_atmosphere_db, param_lims=[('teff',400,2500)], rebin_models=True, grid_resolution=[('teff',50),('logg',0.5)]):
#   '''
#   Given a **model_grid_name**, returns the grid from the model_atmospheres.db in the proper format to work with fit_spectrum()
#   '''
#   # Load the model_atmospheres database and pull all the data from the specified table
#   db = BDdb.get_db(model_atmosphere_db)
#   if param_lims:
#     limit_text = ' AND '.join(['{} BETWEEN {} AND {}'.format(i,j,k) for i,j,k in param_lims])
#     model_grid = db.dict("SELECT * FROM {} WHERE {}".format(model_grid_name,limit_text)).fetchall()  
#   else: model_grid = db.dict("SELECT * FROM {}".format(model_grid_name)).fetchall()
# 
#   # Load the model atmospheres into a data frame
#   models = pd.DataFrame(model_grid)
# 
#   # Choose template wavelength array to rebin all other spectra
#   W = rebin_models if isinstance(rebin_models,(list,np.ndarray)) else models['wavelength'][0]
# 
#   # Rebin model spectra
#   def rebin(row): return [i.value for i in u.rebin_spec([row['wavelength']*q.um, row['flux']*q.erg/q.s/q.cm**2/q.AA], W*q.um)][:2]
#   models['wavelength'], models['flux'] = [[i[n] for i in models.apply(rebin, axis=1)] for n in [0,1]]
# 
#   return models

def model_grid_interp_test(model_grid, teff, logg):
  '''
  Perform a tet of the model grid interpolation by specifying the teff, logg and grid resolution of the given model grid.
  
  Parameters
  ----------
  model_grid: dict
    The model grid object that results from running make_model_db()
  teff: tuple, list
    A sequence of the teff value and increment over which to interpolate, e.g. (1200,100) tests the model at 1200K by interpolating between the 1100K and 1300K models
  logg: tuple, list
    A sequence of the logg value and increment over which to interpolate, e.g. (4.5,0.5) tests the model at 4.5dex by interpolating between the 4.0dex and 5.0dex models

  Returns
  -------
  None
  '''
  # Get the upper, lower and target teff and logg values
  (t1,g1), (t2,g2), (t3,g3) = [(teff[0]+(teff[1]*i),logg[0]+(logg[1]*i)) for i in [-1.,0.,1.]]
  
  # Find the indexes of the models with the appropriate parameters
  idx1, idx2, idx3 = [zip(model_grid['teff'],model_grid['logg']).index((t,g)) for t,g in [(t1,g1), (t2,g2), (t3,g3)]]
  
  # Pull out the wavelength and flux values of the appropriate models
  (w1,f1), (w2,f2), (w3,f3) = [(model_grid['wavelength'][i],model_grid['flux'][i]) for i in [idx1,idx2,idx3]]
  
  # Pop the true model and interpolate the model grid to see how they compare
  mg = model_grid.copy()
  mg['teff'], mg['logg'], mg['wavelength'], mg['flux'] = [np.delete(mg[p], idx2, 0) for p in ['teff','logg','wavelength','flux']]
  W, F = interp_models(['teff','logg'], [t2,g2], mg)
  
  # Plot the interpolated model and the upper, lower, anf true models
  plt.loglog(W, F, color='k', lw=2, label='Interpolated')
  for w,f,c,t,g in [[w1,f1,'r',t1,g1],[w2,f2,'g',t2,g2],[w3,f3,'b',t3,g3]]:
    f = u.smooth(f, 1)
    plt.loglog(w.value, f.value, color=c, label='{} {}'.format(t,g), alpha=0.5)
  plt.legend(loc=0, fontsize=14, frameon=False)

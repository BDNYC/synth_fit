## START HERE FOR MCMC LOOP OVER LIST OF OBJECTS

import run_mcmc
import numpy as np 
import astropy.units as u
from matplotlib import pylab as plt

def showme():
	output=np.load('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/model_fits/spex_objects.txt')
	'''output=[object,source,spectype,shortname]'''
	notworking, done=[],[]
	for i in output[:1]:	
	
## Why doesn't this work?
# 		flux = (u.erg/u.AA/u.cm**2/u.s)*i[0][1]
# 		print flux.value
# 		print flux.unit
# 		wavelength = (u.um)*i[0][0]
# 		raw_spectrum=[wavelength,flux,i[0][2]]

		try: 
			run_mcmc.fit_spectrum(i[1])
			done.append(i[1])
# 			plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/t_dwarfs/mcmc_fits/Marley_Saumon/triangle_plots/'+'{}_{}'.format(i[2],i[3])+'.pdf')

			# Turn the model_atmospheres.db grid into one that can talk to Steph's code
	model_grid = make_model_db(model_grid_name, model_atmosphere_db)
model_atmosphere_db='/Users/paigegiorla/Code/Python/BDNYC/model_atmospheres.db',
		except ValueError: notworking.append(i[1])
		except TypeError: notworking.append(i[1])
# 		except UnitsError: notworking.append(i[1])
## How do I make it print here?		
	
	return notworking, done	
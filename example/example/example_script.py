## Example script for synth_fit

import logging

import cPickle
import astropy.units as u

import BDdb
import synth_fit
import synth_fit.bdfit

logging.basicConfig(level=logging.INFO)

# load the database - replace with appropriate path
db = BDdb.get_db('/home/stephanie/Dropbox/BDNYC_new.db')
#db = BDdb.get_db('/home/stephanie/Dropbox/BDNYCdb/BDNYC.db')

object_name = '0036+1821'

query_spectrum = db.dict.execute(
        "SELECT spec.wavelength_units, spec.flux_units, "
        "spec.wavelength, spec.flux, spec.unc, spec.header FROM "
        "spectra AS spec JOIN sources as s ON spec.source_id=s.id "
        "WHERE s.shortname='0036+1821' AND spec.telescope_id=7 AND "
        "spec.instrument_id=6 AND mode_id=1").fetchone()

# turn into a dictionary with astropy units quantities
wave_unit = u.Unit(query_spectrum['wavelength_units'])
flux_unit = u.Unit(query_spectrum['flux_units'].replace("ergss","erg s"
                                                        ).replace('ergs','erg s'
                                                        ).replace('Wm','W m'
                                                        ).replace("A","AA"))

spectrum = {'wavelength': query_spectrum['wavelength']*wave_unit,
            'flux': query_spectrum['flux']*flux_unit,
            'unc': query_spectrum['unc']*flux_unit}

# now open up the model file
infile = open("SpeX_marley_nolowg.pkl","rb")
model = cPickle.load(infile)
infile.close()

# change into astropy units quantities
model['fsyn'] = model['fsyn']*(u.erg / (u.AA * u.cm**2 * u.s))
model['wsyn'] = model['wsyn']*(u.um)

params = ['logg', 'fsed', 'teff']

# now set up the sampler object (it's a wrapper around emcee)

plot_title = "TESTING {}, {}".format(object_name,"Saumon & Marley 08")

bdsamp = synth_fit.bdfit.BDSampler(object_name, 
                                   spectrum,
                                   model,
                                   params,
                                   smooth=False, # model already matches data
                                   plot_title=plot_title,
                                   snap=True) # no interpolation on grid

#this isn't enough for real results, just to make sure it runs
bdsamp.mcmc_go(nwalk_mult=2,nstep_mult=10)

logging.info("ran MCMC")

logging.info("all done!")

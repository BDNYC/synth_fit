# Module containing functions for working with emcee and running mcmc
# and plotting the output
# Stephanie Douglas
################################################################################

import datetime
import logging

## Third-party
import matplotlib

matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import emcee

from plotting.emcee_plot import emcee_plot
# https://github.com/adrn/streams/blob/master/streams/plot/emcee.py
# Want to update ^ so it shows the burn_in cut
from plotting import triangle
from make_model import *
from calc_chisq import *


class BDSampler(object):
    """
    Class to contain and run emcee on a spectrum and model grid
    and to plot/analyze outputs

    Call as:
       from bdmcmc import bdfit
       x = bdfit.BDSampler(obj_name,spectrum,model,params)

    Parameters for __init__
    -----------------------
    obj_name: string
        gives an identifier for the object

    spectrum: dictionary 
        contains 'wavelength','flux','unc' arrays
        (all much be astropy.units Quantities)

    model: dictionary 
        keys 'wsyn' and 'fsyn' should correspond to model wavelength and 
        flux arrays, and those should be astropy.units Quantities
        other keys should correspond to params

    params: list of strings
        parameters to vary in fit, must be keys of model

    smooth: boolean (default=True)
        whether or not to smooth the model spectra before interpolation 
        onto the data wavelength grid 

    plot_title (string,default='None')
        title for any plots created; also used as part of filenames for 
        output files

    Creates
    -------
    date (string)
    name (string)
    plot_title (string)
    model (bdmcmc.make_model.ModelGrid instance)
    model_ndim (integer) : number of parameters in model
    start_p (array_like) : starting parameters
    all_params (list) : all parameters (model params + ln(s))
    ndim (int) : total number of parameters

    mcmc_go(nwalk_mult=20, nstep_mult=50, outfile=None):
        chain (array_like)
        cropchain (array_like)

    plot_triangle():
        corner_fig (pyplot figure)
    plot_chains():
        chain_fig (pyplot figure)

    get_error_and_unc():
        all_quantiles (array_like) : 
        means (array_like) : 
        lower_lims (array_like) : 
        upper_lims (array_like) : 
        error_and_unc (array_like) : 

    Automatic Outputs
    -----------------
    pickled output file (filename is plot_title + '_chains.pkl')

    """

    def __init__(self, obj_name, spectrum, model, params, smooth=False,
                 plot_title='None', snap=False, wavelength_bins=[0.9, 1.4, 1.9, 2.5] * u.um):
        """
        Parameters 
        ----------
        obj_name: string
            gives an identifier for the object
    
        spectrum: dictionary 
            contains 'wavelength','flux','unc' arrays
            (all much be astropy.units Quantities)
    
        model: dictionary 
            keys 'wsyn' and 'fsyn' should correspond to model wavelength and 
            flux arrays, and those should be astropy.units Quantities
            other keys should correspond to params
        
        params: list of strings
            parameters to vary in fit, must be keys of model

        smooth: boolean (default=True)
            whether or not to smooth the model spectra before interpolation 
            onto the data wavelength grid 

        plot_title (string, default='None')
            title for any plots created; also used as part of filenames for 
            output files. If none is provided, object name and date are used


        """

        ## date string to version output files for a particular run
        self.date = datetime.date.isoformat(datetime.date.today())
        # Eventually - Add a timestamp?

        self.snap = snap
        self.name = obj_name
        logging.info('%s', self.name)

        ## If no plot_title is provided, create one
        if plot_title == 'None':
            self.plot_title = '{} {}'.format(self.name, self.date)
        else:
            self.plot_title = plot_title

        ## Set up the ModelGrid instance (this contains the data and 
        ## model dictionary. It is passed to emcee, and is used to 
        ## calculate the probabilities during the MCMC run)
        self.model = ModelGrid(spectrum, model, params, smooth=smooth,
                               snap=snap, wavelength_bins=wavelength_bins)
        # print spectrum.keys()
        logging.info('Set model')

        ## Calculate the number of parameters for the atmospheric model
        self.model_ndim = len(params)
        logging.info('{} params {}'.format(self.model_ndim,
                                           str(params)))

        ## Calculate starting parameters for the emcee walkers 
        ## by minimizing chi-squared just using the grid of synthetic spectra
        self.start_p, self.min_chi = test_all(spectrum['wavelength'], spectrum['flux'],
                                              spectrum['unc'], model, params, smooth=smooth, shortname=obj_name)
        for i in range(self.model_ndim):
            if (self.start_p[i] >= self.model.plims[params[i]]['max']):
                self.start_p[i] = self.start_p[i] * 0.95
            elif (self.start_p[i] <= self.model.plims[params[i]]['min']):
                self.start_p[i] = self.start_p[i] * 1.05
        print 'chisq:', self.start_p
        ## Add additional parameters beyond the atmospheric model parameters
        self.all_params = list(np.copy(params))

        if len(wavelength_bins) > 1:
            norm_number = len(wavelength_bins) - 1
        else:
            norm_number = 1
        for i in range(norm_number):
            self.all_params.append("N{}".format(i))

        # add normalization parameter
        self.start_p = np.append(self.start_p, np.ones(norm_number))

        # add (log of) tolerance parameter
        good_unc = np.where(np.isnan(self.model.unc) == False)[0]
        start_lns = np.log(2.0 * np.average(self.model.unc[good_unc]))
        logging.info('starting ln(s)={} s={}'.format(start_lns,
                                                     np.exp(start_lns)))
        self.start_p = np.append(self.start_p, start_lns)
        self.all_params.append("ln(s)".format(i))

        logging.info('All params %s', str(self.all_params))
        logging.info('Set starting params %s', str(self.start_p))

        ## The total number of dimensions for the fit is the number of
        ## parameters for the model plus any additional parameters added above
        self.ndim = len(self.all_params)

    def mcmc_go(self, nwalk_mult=20, nstep_mult=50, outfile=None):
        """
        Sets up and calls emcee to carry out the MCMC algorithm

        Parameters
        ----------
        nwalk_mult: integer (default=20)
            multiplied by ndim to get the number of walkers

        nstep_mult: integer (default=50)
            multiplied by ndim to get the number of steps

        outfile: string (default=None)
            filename for any output files; if none is provided, use plot_title

        Creates
        -------
        self.chain (output of all chains)
        self.cropchain (cuts out the first 10% of the steps, 
            then flattens the chain)
        """

        nwalkers, nsteps = self.ndim * nwalk_mult, self.ndim * nstep_mult
        logging.info('%d walkers, %d steps', nwalkers, nsteps)

        ## Initialize the walkers in a gaussian ball around start_p
        ## start_p was set in __init, with the minimum chi-squared model
        ## plus any additional parameters
        p0 = np.zeros((nwalkers, self.ndim))
        logging.debug('p0 shape %s', str(np.shape(p0)))
        for i in range(nwalkers):
            p0[i] = self.start_p + (1e-2 * np.random.randn(self.ndim) *
                                    self.start_p)
            logging.debug('p0[%s] shape %s', i, str(p0[i]))

        ## Set up the sampler
        sampler = emcee.EnsembleSampler(nwalkers, self.ndim, self.model)
        logging.info('sampler set')

        ## Burn in the walkers
        pos, prob, state = sampler.run_mcmc(p0, nsteps / 10)
        logging.debug('pos %s', str(pos))
        logging.debug('prob %s', str(prob))
        logging.debug('state %s', str(state))

        ## Reset the walkers, so the burn-in steps aren't included in analysis
        ## Now the walkers start at the position from the end of the burn-in
        ## Then run the actual MCMC run
        sampler.reset()
        logging.info('sampler reset')
        pos, prob, state = sampler.run_mcmc(pos, nsteps)
        logging.info('sampler completed')
        logging.info("avg accept {}".format(np.average(
            sampler.acceptance_fraction)))
        # logging.info("avg autocorrelation length {}".format(np.average(
        #    sampler.acor)))

        ## store chains for plotting/analysis
        ## Chains contains the positions for each parameter, for each walker
        self.chain = sampler.chain

        ## cut out the burn-in samples (first 10%, for now)
        burn_in = np.floor(nsteps * 0.1)
        self.cropchain = sampler.chain[:, int(burn_in):, :].reshape(
            (-1, self.ndim))

        if self.snap:
            chain_shape = np.shape(self.chain[:, burn_in:, :])
            logging.debug("starting to snap {}".format(chain_shape))
            self.cropchain = self.model.snap_full_run(self.cropchain)
            logging.debug("Snapped cropchains {} to {}".format(
                chain_shape, np.shape(self.cropchain)))
            self.chain = self.cropchain.reshape(chain_shape)
            logging.debug("Snapped chains")

        ## Save the chains to a pkl file for any diagnostics
        if outfile == None:
            outfile = '{}_chains.pkl'.format(self.plot_title)
        open_outfile = open(outfile, 'wb')
        cPickle.dump(self.chain, open_outfile)
        open_outfile.close()

        ## Reshape the chains (don't need to crop out burn-in b/c that's done)
        ## This makes one array with all the samples for each parameter
        self.cropchain = sampler.chain.reshape((-1, self.ndim))
        self.get_quantiles()

    def plot_triangle(self, extents=None):
        """
        Calls triangle module to create a corner-plot of the results
        """
        self.corner_fig = triangle.corner(self.cropchain, labels=self.all_params, quantiles=[.16, .5, .84],
                                          verbose=False, extents=extents)  # , truths=np.ones(3))
        plt.suptitle(self.plot_title)

    def plot_chains(self):
        """
        Calls Adrian's code to plot the development of the chains
        as well as 1D histograms of the results
        """
        self.chain_fig = emcee_plot(self.chain, labels=self.all_params)
        plt.suptitle(self.plot_title)

    def quantile(self, x, quantiles):
        """
        Calculate the quantiles given by quantiles for the array x

        From DFM's triangle code
        """
        xsorted = sorted(x)
        qvalues = [xsorted[int(q * len(xsorted))] for q in quantiles]
        return zip(quantiles, qvalues)

    def get_quantiles(self):
        """ calculates (16th, 50th, 84th) quantiles for all parameters """
        self.all_quantiles = np.ones((self.ndim, 3)) * -99.
        for i in range(self.ndim):
            quant_array = self.quantile(self.cropchain[:, i], [.16, .5, .84])
            self.all_quantiles[i] = [quant_array[j][1] for j in range(3)]

    def get_error_and_unc(self):
        """ Calculates 1-sigma uncertainties for all parameters """

        self.get_quantiles()

        ## The 50th quantile is the mean, the upper and lower "1-sigma" 
        ## uncertainties are calculated from the 16th- and 84th- quantiles
        ## in imitation of Gaussian uncertainties
        self.means = self.all_quantiles[:, 1]
        self.lower_lims = self.all_quantiles[:, 2] - self.all_quantiles[:, 1]
        self.upper_lims = self.all_quantiles[:, 1] - self.all_quantiles[:, 0]

        self.error_and_unc = np.ones((self.ndim, 3)) * -99.
        self.error_and_unc[:, 1] = self.all_quantiles[:, 1]
        self.error_and_unc[:, 0] = (self.all_quantiles[:, 2] -
                                    self.all_quantiles[:, 1])
        self.error_and_unc[:, 2] = (self.all_quantiles[:, 1]
                                    - self.all_quantiles[:, 0])

        return self.error_and_unc

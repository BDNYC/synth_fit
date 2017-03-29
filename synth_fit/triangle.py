#!/usr/bin/env python
# -*- coding: utf-8 -*-
#edited by Stephanie Douglas

from __future__ import print_function, absolute_import, unicode_literals

__all__ = ["corner", "hist2d", "error_ellipse"]
__version__ = "0.0.5"
__author__ = "Dan Foreman-Mackey (danfm@nyu.edu)"
__copyright__ = "Copyright 2013 Daniel Foreman-Mackey"
__contributors__ = [
    # Alphabetical by first name.
    "Adrian Price-Whelan @adrn",
    "Brendon Brewer @eggplantbren",
    "Stephanie Douglas @stephtdouglas",
    "Ekta Patel @ekta1224",
    "Emily Rice @emilurice",
    "Geoff Ryan @geoffryan",
    "Phil Marshall @drphilmarshall",
    "Pierre Gratier @pirg",]
import logging
import numpy as np
# import matplotlib
# matplotlib.use('WXAgg')
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap, colorConverter
from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

def corner(xs, labels=None, extents=None, truths=None, truth_color="#4682b4",
           scale_hist=False, quantiles=[], verbose=True, plot_contours=True,
           plot_datapoints=False, spec_grid=None, demand_fig=None, **kwargs):
    """
    Make a *sick* corner plot showing the projections of a data set in a
    multi-dimensional space. kwargs are passed to hist2d() or used for
    `matplotlib` styling.

    Parameters
    ----------
    xs : array_like (nsamples, ndim)
        The samples. This should be a 1- or 2-dimensional array. For a 1-D
        array this results in a simple histogram. For a 2-D array, the zeroth
        axis is the list of samples and the next axis are the dimensions of
        the space.

    labels : iterable (ndim,) (optional)
        A list of names for the dimensions.

    extents : iterable (ndim,) (optional)
        A list of length 2 tuples containing lower and upper bounds (extents)
        for each dimension, e.g., [(0.,10.), (1.,5), etc.]

    truths : iterable (ndim,) (optional)
        A list of reference values to indicate on the plots.

    truth_color : str (optional)
        A ``matplotlib`` style color for the ``truths`` makers.

    scale_hist : bool (optional)
        Should the 1-D histograms be scaled in such a way that the zero line
        is visible?

    quantiles : iterable (optional)
        A list of fractional quantiles to show on the 1-D histograms as
        vertical dashed lines.

    verbose : bool (optional)
        If true, print the values of the computed quantiles.

    plot_contours : bool (optional)
        Draw contours for dense regions of the plot.

    plot_datapoints : bool (optional)
        Draw the individual data points.

    spec_grids : matplotlib.gridspec.GridSpec (optional)
        array of 

    fig : matplotlib.Figure (optional)
        Overplot onto the provided figure object.
    """

    # Deal with 1D sample lists.
    xs = np.atleast_1d(xs)
    if len(xs.shape) == 1:
        xs = np.atleast_2d(xs)
    else:
        assert len(xs.shape) == 2, "The input sample array must be 1- or 2-D."
        xs = xs.T
    assert xs.shape[0] <= xs.shape[1], "I don't believe that you want more " \
                                       "dimensions than samples!"

    # backwards-compatibility
    plot_contours = kwargs.get("smooth", plot_contours)

    K = len(xs)
    factor = 2.0           # size of one side of one panel
    lbdim = 0.5 * factor   # size of left/bottom margin
    trdim = 0.05 * factor  # size of top/right margin
    whspace = 0.05         # w/hspace size
    plotdim = factor * K + factor * (K - 1.) * whspace
    dim = lbdim + plotdim + trdim

    if spec_grid is None:
        fig = pl.figure()
        spec_grid = gridspec.GridSpec(K,K)
    else:
        fig=pl.gcf()

    #set up a full grid for ease; unused spots will be set to invisible later

    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    #pl.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
    #                    wspace=whspace, hspace=whspace)

    if demand_fig is None:
        setup_axes = [[pl.subplot(spec_grid[i,j]) for j in np.arange(K)] 
                  for i in np.arange(K)]
        axes = np.array(setup_axes).reshape((K,K))
 
        fig, axes = pl.subplots(K, K, figsize=(dim, dim))
    else:
        try:
          fig,axes = demand_fig[0],demand_fig[1] 
          axes = np.array(axes).reshape((K, K))
          print(K)
        except:
            raise ValueError("Provided figure has {0} axes, but data has "
                             "other dimensions K={1}".format(len(demand_fig[1]), K))

    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
                        wspace=whspace, hspace=whspace)


    if extents is None:
        extents = [[x.min(), x.max()] for x in xs]

        # Check for parameters that never change.
        m = np.array([e[0] == e[1] for e in extents], dtype=bool)
        if np.any(m):
        	idxs = np.where(m==True)[0]
        	for i in idxs:
        		if min(xs[i])==0:      			
        			extents[i] = [min(xs[i])-1,max(xs[i])+1]
        		else:
        			extents[i] = [min(xs[i])*0.99,max(xs[i])*1.01]
#             raise ValueError(("It looks like the parameter(s) in column(s) "
#                               "{0} have no dynamic range. Please provide an "
#                               "`extent` argument.")
#                              .format(", ".join(map("{0}".format,
#                                                    np.arange(len(m))[m]))))

    for i, x in enumerate(xs):
        ax = axes[i, i]
#         # Plot the histograms. If it's not the first histogram, rotate 90 degrees with orientation='horizontal'.
#         orientation = 'vertical'
#         if i==1:
#             orientation = 'horizontal'
#         if i==2:
#             orientation = 'horizontal'    

        n, b, p = ax.hist(x, bins=kwargs.get("bins", 50), range=extents[i],
                          histtype="step", color=kwargs.get("color", "k"), lw=kwargs.get("linewidth"))
        if truths is not None:
            ax.axvline(truths[i], color=truth_color)

        # Plot quantiles if wanted.
        if len(quantiles) > 0:
            xsorted = sorted(x)
            qvalues = [xsorted[int(q * len(xsorted))] for q in quantiles]
            for q in qvalues:
                ax.axvline(q, ls="dashed", color=kwargs.get("color", "k"))

            if verbose:
                print("Quantiles:")
                print(zip(quantiles, qvalues))

        # Set up the axes.
        ax.set_xlim(extents[i])
        if scale_hist:
            maxn = np.max(n)
            ax.set_ylim(-0.1 * maxn, 1.1 * maxn)
        else:
            ax.set_ylim(0, 1.1 * np.max(n))
        ax.set_yticklabels([])
        ax.xaxis.set_major_locator(MaxNLocator(5))

        # Not so DRY.
        if i < K - 1:
            ax.set_xticklabels([])
        else:
            [l.set_rotation(30) for l in ax.get_xticklabels()]
            if labels is not None:
                ax.set_xlabel(labels[i],fontsize="x-large")
                ax.xaxis.set_label_coords(0.56, -0.19)

        for j, y in enumerate(xs):
            ax = axes[i, j]
            if j > i:
                ax.set_visible(False)
                ax.set_frame_on(False)
                continue
            elif j == i:
                continue

            hist2d(y, x, ax=ax, extent=[extents[j], extents[i]],
                   plot_contours=plot_contours,
                   plot_datapoints=plot_datapoints,
                   **kwargs)

            if truths is not None:
                ax.plot(truths[j], truths[i], "s", color=truth_color)
                ax.axvline(truths[j], color=truth_color)
                ax.axhline(truths[i], color=truth_color)

            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(5))

            if i < K - 1:
                ax.set_xticklabels([])
            else:
                [l.set_rotation(30) for l in ax.get_xticklabels()]
                if labels is not None:
                    ax.set_xlabel(labels[j],fontsize="x-large")
                    ax.xaxis.set_label_coords(0.43, -0.19)

            if j > 0:
                ax.set_yticklabels([])
            else:
                [l.set_rotation(30) for l in ax.get_yticklabels()]
                if labels is not None:
                    ax.set_ylabel(labels[i],fontsize="x-large")
                    ax.yaxis.set_label_coords(-0.19, 0.55)

    return fig, axes

def error_ellipse(mu, cov, ax=None, factor=1.0, **kwargs):
    """
    Plot the error ellipse at a point given its covariance matrix.

    """
    # some sane defaults
    facecolor = kwargs.pop('facecolor', 'none')
    edgecolor = kwargs.pop('edgecolor', 'k')

    x, y = mu
    U, S, V = np.linalg.svd(cov)
    theta = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
    ellipsePlot = Ellipse(xy=[x, y],
                          width=2 * np.sqrt(S[0]) * factor,
                          height=2 * np.sqrt(S[1]) * factor,
                          angle=theta,
                          facecolor=facecolor, edgecolor=edgecolor, **kwargs)

    if ax is None:
        ax = pl.gca()
    ax.add_patch(ellipsePlot)

    return ellipsePlot

def hist2d(x, y, bins=20, range=None, weights=None, levels=None, smooth=None,
           ax=None, color=None, plot_datapoints=False, plot_density=False,
           plot_contours=True, no_fill_contours=False, fill_contours=False,
           contour_kwargs=None, contourf_kwargs=None, data_kwargs=None,
           **kwargs):
    """
    Plot a 2-D histogram of samples.
    Parameters
    ----------
    x, y : array_like (nsamples,)
       The samples.
    levels : array_like
        The contour levels to draw.
    ax : matplotlib.Axes (optional)
        A axes instance on which to add the 2-D histogram.
    plot_datapoints : bool (optional)
        Draw the individual data points.
    plot_density : bool (optional)
        Draw the density colormap.
    plot_contours : bool (optional)
        Draw the contours.
    no_fill_contours : bool (optional)
        Add no filling at all to the contours (unlike setting
        ``fill_contours=False``, which still adds a white fill at the densest
        points).
    fill_contours : bool (optional)
        Fill the contours.
    contour_kwargs : dict (optional)
        Any additional keyword arguments to pass to the `contour` method.
    contourf_kwargs : dict (optional)
        Any additional keyword arguments to pass to the `contourf` method.
    data_kwargs : dict (optional)
        Any additional keyword arguments to pass to the `plot` method when
        adding the individual data points.
    """
    if ax is None:
        ax = pl.gca()

    # Set the default range based on the data range if not provided.
    if range is None:
        if "extent" in kwargs:
            logging.warn("Deprecated keyword argument 'extent'. "
                         "Use 'range' instead.")
            range = kwargs["extent"]
        else:
            range = [[x.min(), x.max()], [y.min(), y.max()]]

    # Set up the default plotting arguments.
    if color is None:
        color = "k"

    # Choose the default "sigma" contour levels.
    if levels is None:
        levels = 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5) ** 2)

    # This is the color map for the density plot, over-plotted to indicate the
    # density of the points near the center.
    density_cmap = LinearSegmentedColormap.from_list(
        "density_cmap", [color, (1, 1, 1, 0)])

    # This color map is used to hide the points at the high density areas.
    white_cmap = LinearSegmentedColormap.from_list(
        "white_cmap", [(1, 1, 1), (1, 1, 1)], N=2)

    # This "color map" is the list of colors for the contour levels if the
    # contours are filled.
    rgba_color = colorConverter.to_rgba(color)
    contour_cmap = [list(rgba_color) for l in levels] + [rgba_color]
    for i, l in enumerate(levels):
        contour_cmap[i][-1] *= float(i) / (len(levels)+1)

    # We'll make the 2D histogram to directly estimate the density.
    try:
        H, X, Y = np.histogram2d(x.flatten(), y.flatten(), bins=bins,
                                 range=range, weights=weights)
    except ValueError:
        raise ValueError("It looks like at least one of your sample columns "
                         "have no dynamic range. You could try using the "
                         "'range' argument.")

    if smooth is not None:
        if gaussian_filter is None:
            raise ImportError("Please install scipy for smoothing")
        H = gaussian_filter(H, smooth)

    # Compute the density levels.
    Hflat = H.flatten()
    inds = np.argsort(Hflat)[::-1]
    Hflat = Hflat[inds]
    sm = np.cumsum(Hflat)
    sm /= sm[-1]
    V = np.empty(len(levels))
    for i, v0 in enumerate(levels):
        try:
            V[i] = Hflat[sm <= v0][-1]
        except:
            V[i] = Hflat[0]
    V.sort()
    m = np.diff(V) == 0
    if np.any(m):
        logging.warning("Too few points to create valid contours")
    while np.any(m):
        V[np.where(m)[0][0]] *= 1.0 - 1e-4
        m = np.diff(V) == 0
    V.sort()

    # Compute the bin centers.
    X1, Y1 = 0.5 * (X[1:] + X[:-1]), 0.5 * (Y[1:] + Y[:-1])

    # Extend the array for the sake of the contours at the plot edges.
    H2 = H.min() + np.zeros((H.shape[0] + 4, H.shape[1] + 4))
    H2[2:-2, 2:-2] = H
    H2[2:-2, 1] = H[:, 0]
    H2[2:-2, -2] = H[:, -1]
    H2[1, 2:-2] = H[0]
    H2[-2, 2:-2] = H[-1]
    H2[1, 1] = H[0, 0]
    H2[1, -2] = H[0, -1]
    H2[-2, 1] = H[-1, 0]
    H2[-2, -2] = H[-1, -1]
    X2 = np.concatenate([
        X1[0] + np.array([-2, -1]) * np.diff(X1[:2]),
        X1,
        X1[-1] + np.array([1, 2]) * np.diff(X1[-2:]),
    ])
    Y2 = np.concatenate([
        Y1[0] + np.array([-2, -1]) * np.diff(Y1[:2]),
        Y1,
        Y1[-1] + np.array([1, 2]) * np.diff(Y1[-2:]),
    ])

    if plot_datapoints:
        if data_kwargs is None:
            data_kwargs = dict()
        data_kwargs["color"] = data_kwargs.get("color", color)
        data_kwargs["ms"] = data_kwargs.get("ms", 2.0)
        data_kwargs["mec"] = data_kwargs.get("mec", "none")
        data_kwargs["alpha"] = data_kwargs.get("alpha", 0.1)
        ax.plot(x, y, "o", zorder=-1, rasterized=True, **data_kwargs)

    # Plot the base fill to hide the densest data points.
    if (plot_contours or plot_density) and not no_fill_contours:
        ax.contourf(X2, Y2, H2.T, [V.min(), H.max()],
                    cmap=white_cmap, antialiased=False)

    if plot_contours and fill_contours:
        if contourf_kwargs is None:
            contourf_kwargs = dict()
        contourf_kwargs["colors"] = contourf_kwargs.get("colors", contour_cmap)
        contourf_kwargs["antialiased"] = contourf_kwargs.get("antialiased",
                                                             False)
        ax.contourf(X2, Y2, H2.T, np.concatenate([[0], V, [H.max()*(1+1e-4)]]),
                    **contourf_kwargs)

    # Plot the density map. This can't be plotted at the same time as the
    # contour fills.
    elif plot_density:
        ax.pcolor(X, Y, H.max() - H.T, cmap=density_cmap)

    # Plot the contour edge colors.
    if plot_contours:
        if contour_kwargs is None:
            contour_kwargs = dict()
        contour_kwargs["colors"] = contour_kwargs.get("colors", color)
        ax.contour(X2, Y2, H2.T, V, **contour_kwargs)

    ax.set_xlim(range[0])
    ax.set_ylim(range[1])
# Functions to calculate eddy viscosity according to Scott et al. 2023
# includes plotting functions for viewing the scatter plots and yz slices

import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def scatter_yz(io, x, io2=None, ylim=[-1.15, 1.15], zlim=[-1/1.75, 1/1.75], x_key='S13', y_key='uw_wake', fit=True, xscale = 1):
    """
    Scatter plot at a given x location with the linear fit

    INPUTS:
    io (BudgetIO or DeficitIO)    - budget object to plot
    x (int or 1D array)           - x location of the y-z slice
    (optional)
    io2 (BudgetIO or DeficitIO)   - if used io2 should be the precursor budget and io should be the primary budget
    ylim (1D array)               - y limits of the y-z slice
    zlim  (1D array)              - z limits of the y-z slice
    x_key (string)                - budget key for the x variable
    y_key (string)                - budget key for the y variable
    fit (boolean)                 - True if fit should be plotted 

    OUTPUTS
    fig, ax 
    """

    if isinstance(x, list) or isinstance(x, np.ndarray):
        fig, ax = plt.subplots(1, len(x), layout='tight')

        for i in range(len(x)):
            xid, yid, zid = io.get_xids(x=x[i], y=ylim, z=zlim, return_none=True, return_slice=True)

            if io2:
                x_data = io.budget[x_key][xid,yid,zid] - io2.budget[x_key][xid,yid,zid]
                y_data = io.budget[y_key][xid,yid,zid] - io2.budget[y_key][xid,yid,zid]
            else:
                x_data = io.budget[x_key][xid,yid,zid]
                y_data = io.budget[y_key][xid,yid,zid]

            # scatter plot
            ax[i].scatter(x_data, y_data, alpha=0.5)


            if fit:
                slope, intercept, r_value, p_value, std_err = linregress(x_data.ravel(), y_data.ravel())
                # plot linear fit
                ax[i].plot(x_data.ravel(), x_data.ravel()*slope + intercept, linewidth=3, color='red')
                # ax[i].text(0,np.max(y_data), 'slope = {}'.format(slope))

            ax[i].set_title('x/D = {}'.format(int(x[i]*xscale)))
            ax[i].set_xlabel(x_key)
            ax[i].set_ylabel(y_key)

    else:
        fig, ax = plt.subplots(layout='tight')
        xid, yid, zid = io.get_xids(x=x, y=ylim, z=zlim, return_none=True, return_slice=True)

        if io2:
            x_data = io.budget[x_key][xid,yid,zid] - io2.budget[x_key][xid,yid,zid]
            y_data = io.budget[y_key][xid,yid,zid] - io2.budget[y_key][xid,yid,zid]
        else:
            x_data = io.budget[x_key][xid,yid,zid]
            y_data = io.budget[y_key][xid,yid,zid]

        # scatter plot
        ax.scatter(x_data, y_data, alpha=0.5)


        if fit:
            slope, intercept, r_value, p_value, std_err = linregress(x_data.ravel(), y_data.ravel())
            # plot linear fit
            ax.plot(x_data.ravel(), x_data.ravel()*slope + intercept, linewidth=3, color='red')

        ax.set_title('x/D = {}'.format(int(x*xscale)))
        ax.set_xlabel(x_key)
        ax.set_ylabel(y_key)

    return fig, ax


def slice_yz(io, x, terms, io2=None, ylim=[-1.15, 1.15], zlim=[-1/1.75, 1/1.75], **kwargs):
    """
    Plots y-z slices of given budget terms
    x locations are the columns
    terms are the rows

    INPUTS:
    io (BudgetIO or DeficitIO)    - budget object to plot
    x (1D array)           - x location of the y-z slice
    terms (list of strings)       - list of budget keys to plot
    (optional)
    io2 (BudgetIO or DeficitIO)   - if used io2 should be the precursor budget and io should be the primary budget
    ylim (1D array)               - y limits of the y-z slice
    zlim  (1D array)              - z limits of the y-z slice

    OUTPUTS
    fig, ax 
    """

    fig, ax = plt.subplots(len(terms), len(x), layout='tight')

    if len(terms) > 1 and len(x) > 1:
        for i in range(len(terms)):
            for j in range(len(x)):
                xid, yid, zid = io.get_xids(x=x[j], y=ylim, z=zlim, return_none=True, return_slice=True)
                print(xid)
                if io2:
                    slice = io.budget[terms[i]][xid,yid,zid] - io2.budget[terms[i]][xid,yid,zid]
                else:
                    slice = io.budget[terms[i]][xid,yid,zid]

                if j == 0:
                    clim = [-np.max(slice), np.max(slice)]

                im = ax[i][j].imshow(slice.T, extent=[*ylim, *zlim], origin='lower', clim=clim, **kwargs)
            
                ax[0][j].set_title('x/D = {}'.format(x[j]))

            # create an axes on the right side of ax. The width of cax will be 5%
            # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            divider = make_axes_locatable(ax[i][-1])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(im, cax=cax, label=terms[i])

    elif len(terms) > 1 and len(x) == 1:
        for i in range(len(terms)):
            xid, yid, zid = io.get_xids(x=x, y=ylim, z=zlim, return_none=True, return_slice=True)

            if io2:
                slice = io.budget[terms[i]][xid,yid,zid] - io2.budget[terms[i]][xid,yid,zid]
            else:
                slice = io.budget[terms[i]][xid,yid,zid]

            clim = [-np.max(slice), np.max(slice)]

            im = ax[i].imshow(slice.T, extent=[*ylim, *zlim], origin='lower', clim=clim, **kwargs)

            # create an axes on the right side of ax. The width of cax will be 5%
            # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            divider = make_axes_locatable(ax[i])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(im, cax=cax, label=terms[i])

        ax[0].set_title('x/D = {}'.format(x[0]))




    elif len(terms) == 1 and len(x) > 1:
        for i in range(len(x)):
            xid, yid, zid = io.get_xids(x=x[i], y=ylim, z=zlim, return_none=True, return_slice=True)

            if io2:
                slice = io.budget[terms[0]][xid,yid,zid] - io2.budget[terms[0]][xid,yid,zid]
            else:
                slice = io.budget[terms[0]][xid,yid,zid]

            if i == 0:
                clim = [-np.max(slice), np.max(slice)]
            
            im = ax[i].imshow(slice.T, extent=[*ylim, *zlim], origin='lower', clim=clim, **kwargs)
            
            # create an axes on the right side of ax. The width of cax will be 5%
            # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            divider = make_axes_locatable(ax[-1])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(im, cax=cax, label=terms[0])

            ax[i].set_title('x/D = {}'.format(x[i]))

    else:
        xid, yid, zid = io.get_xids(x=x[0], y=ylim, z=zlim, return_none=True, return_slice=True)

        if io2:
            slice = io.budget[terms[0]][xid,yid,zid] - io2.budget[terms[0]][xid,yid,zid]
        else:
            slice = io.budget[terms[0]][xid,yid,zid]

        im = ax.imshow(slice.T, extent=[*ylim, *zlim], origin='lower', clim = [-np.max(slice), np.max(slice)], **kwargs)

        ax.set_title('x/D = {}'.format(x[0]))


    return fig, ax


def fit_nu_T_slope(io, io2=None, ylim=[-1.15, 1.15], zlim=[-1/1.75, 1/1.75], x_key = 'S13', y_key='uw_wake'):
    """
    Calculates the eddy viscosity from the slope of uw_wake vs S13 in the wake

    INPUTS:
    io (BudgetIO or DeficitIO)    - budget object to plot
    (optional)
    io2 (BudgetIO or DeficitIO)   - if used io2 should be the precursor budget and io should be the primary budget
    ylim (1D array)               - y limits of the y-z slice
    zlim  (1D array)              - z limits of the y-z slice
    x_key (string)                - budget key for the x variable
    y_key (string)                - budget key for the y variable

    OUTPUTS
    slopes (1D array)             - numpy array of slope fits
    slopes_err (1D array)         - numpy array of slope standard errors
    """
    xid, yid, zid = io.get_xids(x=3, y=ylim, z=zlim, return_none=True, return_slice=True)

    slopes = np.zeros(len(io.xLine))
    slope_err = np.zeros(len(io.xLine))

    if io2:
        x_data = io.budget[x_key][:,yid,zid] - io2.budget[x_key][:,yid,zid]
        y_data = io.budget[y_key][:,yid,zid] - io2.budget[y_key][:,yid,zid]
    else:
        x_data = io.budget[x_key][:,yid,zid]
        y_data = io.budget[y_key][:,yid,zid]


    for i in range(len(io.xLine)):
        slope, intercept, r_value, p_value, std_err = linregress(-x_data[i].ravel(), y_data[i].ravel())
        
        slopes[i] = slope
        slope_err[i] = std_err

    return slopes, slope_err

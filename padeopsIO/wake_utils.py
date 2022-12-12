# Additional wake functions (some work-in-progress)
# this is crude but it works
# USAGE: from wake_utils import *

import numpy as np
import csv
import os
import padeopsIO


def wake_centroid_2d(u_hub=None, u_wake_hub=None, y=None, thresh=0): 
    """
    u_hub (nx * ny) : 2D slice at hub height
    """
    
    if u_wake_hub is None: 
        del_u = 1 - u_hub  # assume freestream is normalized to 1
    else: 
        del_u = u_wake_hub
        
    del_u[del_u < thresh] = 0
    
    numer = np.trapz(del_u*y, axis=1)
    denom = np.trapz(del_u, axis=1)
    
    return numer/denom


def wake_centroid_3d(u=None, u_wake=None, y=None, z=None, thresh=0): 
    """
    u (nx * ny * nz) : 3D u-velocity field normalized to geostrophic wind
    u_wake (nx * ny * nz) : 3D wake field with background subtracted already. Supercedes
        `u` if included. 
    y (ny) : y-coordinate axis
    z (nz) : z-coordinate axis
    """
    
    if y is None and z is None: 
        raise ValueError("Either y= or z= must be included. ")
    
    if u_wake is None: 
        del_u = 1 - u  # assume freestream is normalized to 1
    else: 
        del_u = u_wake
        
    # apply thresholding
    del_u[del_u < thresh] = 0
    
    ret = ()
    if y is not None: 
        y = y[np.newaxis, :, np.newaxis]
        numer = np.trapz(np.trapz(del_u*y, axis=1), axis=1)
        denom = np.trapz(np.trapz(del_u, axis=1), axis=1)
        ret += (numer/denom, )
    
    if z is not None: 
        z = z[np.newaxis, np.newaxis, :]
        numer = np.trapz(np.trapz(del_u*z, axis=1), axis=1)
        denom = np.trapz(np.trapz(del_u, axis=1), axis=1)
        ret += (numer/denom, )
        
    if len(ret) == 1: 
        return ret[0]
    else: 
        return ret
    
    
def usq_mean(run, diam=1, xlim=None): 
    """
    Returns the mean squared velocity within the wake region 
    
    Parameters
    ----------
    run : BudgetIO object
    diam : threshold for (square) averaging wake region; constant in x
    xlim : 
    
    Returns
    -------
    ubar_sq : [nx] array of mean squared velocity values
    xlim (opt) : [nx] array xLine associated with the given xlim
    """
    
    sl = run.slice(budget_terms=['ubar'], xlim=xlim, ylim=(run.Ly/2-diam, run.Ly/2+diam), zlim=(run.Lz/2-diam, run.Lz/2+diam))
    ubar_sq = np.mean(np.mean(sl['ubar']**2, axis=1), axis=1)
    
    if xlim is not None: 
        return ubar_sq, sl['x']
    else: 
        return ubar_sq
    
    
def rans_budgets(run, xlim=None, ylim=None, zlim=None, sl=None, 
                 compute_x=True, compute_y=True, compute_z=False, 
                 combine_terms=False, compute_residual=True, 
                 useconstantg=True): 
    """
    Compute RANS budgets from output data with mean fields only. 
    
    Sign convention: advection terms are "moved" to the right-hand side such that they include a negative sign. 
    Pressure gradients already include a negative sign, and density is assumed to be unity. 
    
    Parameters
    ----------
    run (BudgetIO object)
    sl (slice from BudgetIO.slice()) : optional, can be used in place of giving a run and x, y, z limits  (run still required)
    xlim, ylim, zlim : array-like, see BudgetIO.slice()
    compute_x, compute_y, compute_z (bool) : dictates which RANS budgets to compute. Default: x, y True, z False
    combine_terms (bool) : combines advection terms and turbulence terms. Default: False
    compute_residual (bool) : computes the residual pointwise. Default: True
    useconstantg (bool) : uses constant geostrophic forcing given in the `physics` namelist if True. 
        If not, calls _get_inflow()
    
    Returns
    -------
    sl (dict) : slice from run.slice()
    rans_x (dict) : dictionary with the following keys - 
        ['ududx', 'vdudy', 'wdudz', 'dpdx', 'xCor', 'duudx', 'duvdy', 'duwdz']
    rans_y (dict) : dictionary with the following keys - 
        ['udvdx', 'vdvdy', 'wdvdz', 'dpdy', 'yCor', 'dvudx', 'dvvdy', 'dvwdz']
    rans_z (dict) : dictionary with the following keys - 
        ['udwdx', 'vdwdy', 'wdwdz', 'dpdz', 'zCor', 'dwudx', 'dwvdy', 'dwwdz']
        
    # TODO : Clean up this function :)

    """
    
    if sl is None: 
        terms = ['ubar', 'vbar', 'wbar', 'pbar', 'uu', 'uv', 'uw', 'vv', 'vw', 'ww', 'dpdx', 'dpdy', 'dpdz']
        sl = run.slice(budget_terms=terms, xlim=xlim, ylim=ylim, zlim=zlim, round_extent=False)
    else: 
        zlim = [sl['z'][0], sl['z'][-1]]
    
    # constants: 
    Ro = run.Ro
    lat = run.input_nml['physics']['latitude']
    C1 = 2/Ro * np.sin(lat*np.pi/180)
#     Gx = 1 
    # need to retrieve and trim the inflow profiles: 
    if useconstantg: 
        Gx, Gy = [1, 0]  # TODO: GENERALIZE
    else: 
        Gx, Gy = run._get_inflow()
        zids = run.get_xids(z=zlim, return_slice=True)
        Gx = Gx[zids]
        Gy = Gy[zids]

    ret = (sl, )
    
    # x-terms: 
    if compute_x: 
        rans_x = {}
        x_adv_terms = ['ududx', 'vdudy', 'wdudz']
        x_turb_terms = ['duudx', 'duvdy', 'duwdz']

        rans_x['ududx'] = -(sl['ubar'] * partialx(sl['ubar'], run.dx))
        rans_x['vdudy'] = -(sl['vbar'] * partialy(sl['ubar'], run.dy))
        rans_x['wdudz'] = -(sl['wbar'] * partialz(sl['ubar'], run.dz))
        rans_x['dpdx'] = sl['dpdx']  #-partialx(sl['pbar'], run.dx)
        rans_x['xCor'] = -C1*(Gy - sl['vbar'])
        rans_x['duudx'] = -partialx(sl['uu'], run.dx)
        rans_x['duvdy'] = -partialy(sl['uv'], run.dy)
        rans_x['duwdz'] = -partialz(sl['uw'], run.dz)
        
        if compute_residual: 
            rans_x['residual'] = sum(rans_x[key] for key in rans_x.keys())
        
        if combine_terms: 
            rans_xcomb = {key: rans_x[key] for key in rans_x.keys() if key not in (x_adv_terms + x_turb_terms)}
            rans_xcomb['advection'] = sum(rans_x[key] for key in x_adv_terms)
            rans_xcomb['turbulence'] = sum(rans_x[key] for key in x_turb_terms)
            
            ret += (rans_xcomb, )
        else: 
            ret += (rans_x, )

    # y-terms: 
    if compute_y: 
        rans_y = {}
        y_adv_terms = ['udvdx', 'vdvdy', 'wdvdz']
        y_turb_terms = ['dvudx', 'dvvdy', 'dvwdz']

        rans_y['udvdx'] = -(sl['ubar'] * partialx(sl['vbar'], run.dx))
        rans_y['vdvdy'] = -(sl['vbar'] * partialy(sl['vbar'], run.dy))
        rans_y['wdvdz'] = -(sl['wbar'] * partialz(sl['vbar'], run.dz))
        rans_y['dpdy'] = sl['dpdy'] #-partialy(sl['pbar'], run.dy)
        rans_y['yCor'] = C1*(Gx-sl['ubar'])
        rans_y['dvudx'] = -partialx(sl['uv'], run.dx)
        rans_y['dvvdy'] = -partialy(sl['vv'], run.dy)
        rans_y['dvwdz'] = -partialz(sl['vw'], run.dz)
        
        if compute_residual: 
            rans_y['residual'] = sum(rans_y[key] for key in rans_y.keys())
        
        if combine_terms: 
            rans_ycomb = {key: rans_y[key] for key in rans_y.keys() if key not in (y_adv_terms + y_turb_terms)}
            rans_ycomb['advection'] = sum(rans_y[key] for key in y_adv_terms)
            rans_ycomb['turbulence'] = sum(rans_y[key] for key in y_turb_terms)
            
            ret += (rans_ycomb, )
        else: 
            ret += (rans_y, )
        
    # z-terms: 
    if compute_z: 
        rans_z = {}
        z_adv_terms = ['udwdx', 'vdwdy', 'wdwdz']
        z_turb_terms = ['dwudx', 'dwvdy', 'dwwdz']

        rans_z['udwdx'] = -(sl['ubar'] * partialx(sl['wbar'], run.dx))
        rans_z['vdwdy'] = -(sl['vbar'] * partialy(sl['wbar'], run.dy))
        rans_z['wdwdz'] = -(sl['wbar'] * partialz(sl['wbar'], run.dz))
        rans_z['dpdz'] = -partialz(sl['pbar'], run.dz)
        rans_z['zCor'] = np.zeros(sl['ubar'].shape)
        rans_z['dwudx'] = -partialx(sl['uw'], run.dx)
        rans_z['dwvdy'] = -partialy(sl['vw'], run.dy)
        rans_z['dwwdz'] = -partialz(sl['ww'], run.dz)
        
        if compute_residual: 
            rans_z['residual'] = sum(rans_z[key] for key in rans_zkeys())
        
        if combine_terms: 
            rans_zcomb = {key: rans_z[key] for key in rans_z.keys() if key not in (z_adv_terms + z_turb_terms)}
            rans_zcomb['advection'] = sum(rans_z[key] for key in z_adv_terms)
            rans_zcomb['turbulence'] = sum(rans_z[key] for key in z_turb_terms)
            
            ret += (rans_zcomb, )
        else: 
            ret += (rans_z, )

    
    if len(ret) == 1: 
        print("Nothing computed! Check compute_x and compute_y keyword arguments. Returning slice only")
        return ret[0]
    else: 
        return ret
    
    
def get_xids(x=None, y=None, z=None, 
             x_ax=None, y_ax=None, z_ax=None, 
             return_none=False, return_slice=False): 
    """
    COPIED from BudgetIO.py
    
    Translates x, y, and z limits in the physical domain to indices based on x_ax, y_ax, z_ax

    Arguments
    ---------
    x, y, z : float or iterable (tuple, list, etc.) of physical locations to return the nearest index for
    return_none : if True, populates output tuple with None if input is None. Default False. 
    return_slice : if True, returns a tuple of slices instead a tuple of lists. 

    Returns
    -------
    xid, yid, zid : list or tuple of lists with indices for the requested x, y, z, args in the order: x, y, z. 
        If, for example, y and z are requested, then the returned tuple will have (yid, zid) lists. 
        If only one value (float or int) is passed in for e.g. x, then an integer will be passed back in xid. 
    """

    ret = ()

    # iterate through x, y, z, index matching for each term
    for s, s_ax in zip([x, y, z], [x_ax, y_ax, z_ax]): 
        if s is not None: 
            if s_ax is None: 
                raise AttributeError('Axis keyword not providede')
                
            if hasattr(s, '__iter__'): 
                xids = [np.argmin(np.abs(s_ax-xval)) for xval in s]
            else: 
                xids = np.argmin(np.abs(s_ax-s))

            if return_slice:  # append slices to the return tuple
                ret = ret + (slice(np.min(xids), np.max(xids)+1), )

            else:  # append index list to the return tuple
                ret = ret + (xids, )

        elif return_none:  # fill with None or slice(None)
            if return_slice: 
                ret = ret + (slice(None), )

            else: 
                ret = ret + (None, )

    if len(ret)==1: 
        return ret[0]  # don't return a length one tuple 
    else: 
        return ret

# --------------- NUMERICS ----------------


# 2D partial derivatives for x, y
def partialx_2d(field, dx): 
    """
    2nd order central differencing in x
    """
    dfdx = np.zeros(field.shape)
    dfdx[1:-1, :] = (field[2:, :] - field[:-2, :]) / (2*dx)  # central differencing
    dfdx[0, :] = (-3*field[0, :] + 4*field[1, :] - field[2, :]) / (2*dx)  # forward differencing
    dfdx[-1, :] = (3*field[-1, :] - 4*field[-2, :] + field[-3, :]) / (2*dx)  # backward differencing
    
    return dfdx
    

def partialy_2d(field, dy): 
    """
    2nd order central differencing in y
    """
    dfdy = np.zeros(field.shape)
    dfdy[:, 1:-1] = (field[:, 2:] - field[:, :-2]) / (2*dy)  # central differencing
    dfdy[:, 0] = (-3*field[:, 0] + 4*field[:, 1] - field[:, 2]) / (2*dy)  # forward differencing
    dfdy[:, -1] = (3*field[:, -1] - 4*field[:, -2] + field[:, -3]) / (2*dy)  # backward differencing
    
    return dfdy


# 3D partial derivatives for x, y, z
def partialx(field, dx): 
    """
    2nd order central differencing in x
    """
    dfdx = np.zeros(field.shape)
    dfdx[1:-1, :, :] = (field[2:, :, :] - field[:-2, :, :]) / (2*dx)  # central differencing
    dfdx[0, :, :] = (-3*field[0, :, :] + 4*field[1, :, :] - field[2, :, :]) / (2*dx)  # forward differencing
    dfdx[-1, :, :] = (3*field[-1, :, :] - 4*field[-2, :, :] + field[-3, :, :]) / (2*dx)  # backward differencing
    
    return dfdx
    

def partialy(field, dy): 
    """
    2nd order central differencing in y
    """
    dfdy = np.zeros(field.shape)
    dfdy[:, 1:-1, :] = (field[:, 2:, :] - field[:, :-2, :]) / (2*dy)  # central differencing
    dfdy[:, 0, :] = (-3*field[:, 0, :] + 4*field[:, 1, :] - field[:, 2, :]) / (2*dy)  # forward differencing
    dfdy[:, -1, :] = (3*field[:, -1, :] - 4*field[:, -2, :] + field[:, -3, :]) / (2*dy)  # backward differencing
    
    return dfdy


def partialz(field, dz): 
    """
    2nd order central differencing in z
    """
    dfdz = np.zeros(field.shape)
    dfdz[:, :, 1:-1] = (field[:, :, 2:] - field[:, :, :-2]) / (2*dz)  # central differencing
    dfdz[:, :, 0] = (-3*field[:, :, 0] + 4*field[:, :, 1] - field[:, :, 2]) / (2*dz)  # forward differencing
    dfdz[:, :, -1] = (3*field[:, :, -1] - 4*field[:, :, -2] + field[:, :, -3]) / (2*dz)  # backward differencing
    
    return dfdz


# second derivatives in x, y, z
def partialx2(field, dx): 
    """
    second derivative in x
    """
    dfdx = np.zeros(field.shape)
    dfdx[1:-1, :, :] = (field[2:, :, :] - 2*field[1:-1, :, :] + field[:2, :, :]) / (dx**2)  # central differencing
    dfdx[0, :, :] = dfdx[1, :, :]  # for now, cheat at the boundaries... 
    dfdx[-1, :, :] = dfdx[-2, :, :]  # for now, cheat at the boundaries... 
    return dfdx
    

def partialy2(field, dy): 
    """
    second derivative in y
    """
    dfdy = np.zeros(field.shape)
    dfdy[:, 1:-1, :] = (field[:, 2:, :] - 2*field[:, 1:-1:, :] + field[:, :-2, :]) / (dy**2)  # central differencing
    dfdy[:, 0, :] = dfdy[:, 1, :]
    dfdy[:, -1, :] = dfdy[:, -2, :]
    return dfdy


def partialz2(field, dz): 
    """
    second derivative in z
    """
    dfdz = np.zeros(field.shape)
    dfdz[:, :, 1:-1] = (field[:, :, 2:] - 2*field[:, :, 1:-1] + field[:, :, :-2]) / (dz**2)  # central differencing
    dfdz[:, :, 0] = dfdz[:, :, 1]
    dfdz[:, :, -1] = dfdz[:, :, -2]
    
    return dfdz


# ------------------- IO ----------------


def read_list(dir_name): 
    """
    Outputs a list of BudgetIO objects from a directory name by reading the associated
    file `./Runs.csv`
    
    Arguments
    ---------
    dir_name (str) : directory path containing a file Runs.csv
    
    Returns
    -------
    run_list : list of BudgetIO objects
    """
    
    # read runs
    # https://stackoverflow.com/questions/24662571/python-import-csv-to-list
    with open(os.path.join(dir_name, 'Runs.csv')) as f:
        reader = csv.reader(f)
        runs = list(reader)

    run_list = []

    for run in runs: 
        run_list.append(padeopsIO.BudgetIO(run[0], padeops=True, verbose=False))
        
    return run_list


def key_search_r(nested_dict, key):
    """
    Copied from budgetkey.py
    
    Performs a recursive search for the dictionary key `key` in any of the dictionaries contained 
    inside `nested_dict`. Returns the value of nested_dict[key] at the first match. 
    
    Parameters
    ----------
    nested_dict (dict-like) : dictionary [possibly of dictionaries]
    key (str) : dictionary key to match
    
    Returns
    -------
    nested_dict[key] if successful, None otherwise. 
    """
    
    try: 
        for k in nested_dict.keys(): 
            if k == key: 
                return nested_dict[k]
            else: 
                res = key_search_r(nested_dict[k], key)
                if res is not None: 
                    return res
        
    except AttributeError as e: 
        return

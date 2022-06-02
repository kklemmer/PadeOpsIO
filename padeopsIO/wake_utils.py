# Additional wake functions (some work-in-progress)
# this is crude but it works
# USAGE: from wake_utils import *

import numpy as np
import csv


def wake_centroid_2d(u_hub, y, thresh=0): 
    """
    u_hub (nx * ny) : 2D slice at hub height
    """
#     thresh = 0 
    del_u = 1 - u_hub  # assume freestream is normalized to 1
    del_u[del_u < thresh] = 0
    
    numer = np.trapz(del_u*y, axis=1)
    denom = np.trapz(del_u, axis=1)
    
    return numer/denom


def wake_centroid_3d(u, u_wake=None, y=None, z=None, thresh=0): 
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
        del_u = uwake
        
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
    
    
def rans_budgets(run, xlim=None, ylim=None, zlim=None, compute_x=True, compute_y=True): 
    """
    Compute RANS budgets from output data with mean fields only. 
    """
    
    terms = ['ubar', 'vbar', 'wbar', 'pbar', 'uu', 'uv', 'uw', 'vv', 'vw', 'ww']
    sl = run.slice(budget_terms=terms, xlim=xlim, ylim=ylim, zlim=zlim)
    
    # constants: 
    Ro = run.Ro
    lat = run.input_nml['physics']['latitude']
    C1 = 2/Ro * np.sin(lat*np.pi/180)
    Gx = 1 

    ret = ()
    
    # x-terms: 
    if compute_x: 
        rans_x = {}

        rans_x['ududx'] = -(sl['ubar'] * partialx(sl['ubar'], run.dx))
        rans_x['vdudy'] = -(sl['vbar'] * partialy(sl['ubar'], run.dy))
        rans_x['wdudz'] = -(sl['wbar'] * partialz(sl['ubar'], run.dz))
        rans_x['dpdx'] = -partialx(sl['pbar'], run.dx)
        rans_x['xCor'] = C1*sl['vbar']
        rans_x['duudx'] = -partialx(sl['uu'], run.dx)
        rans_x['duvdy'] = -partialy(sl['uv'], run.dy)
        rans_x['duwdz'] = -partialz(sl['uw'], run.dz)
        ret += (rans_x, )

    # y-terms: 
    if compute_y: 
        rans_y = {}

        rans_y['udvdx'] = -(sl['ubar'] * partialx(sl['vbar'], run.dx))
        rans_y['vdvdy'] = -(sl['vbar'] * partialy(sl['vbar'], run.dy))
        rans_y['wdvdz'] = -(sl['wbar'] * partialz(sl['vbar'], run.dz))
        rans_y['dpdy'] = -partialy(sl['pbar'], run.dy)
        rans_y['yCor'] = C1*(Gx-sl['ubar'])
        rans_y['dvudx'] = -partialx(sl['uv'], run.dx)
        rans_y['dvvdy'] = -partialy(sl['vv'], run.dy)
        rans_y['dvwdz'] = -partialz(sl['vw'], run.dz)
        ret += (rans_y, )
    
    if len(ret) == 0: 
        raise ValueError("Nothing computed! Check compute_x and compute_y keyword arguments.")
    elif len(ret) == 1: 
        return ret[0]
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
        run_list.append(pio.BudgetIO(run[0], padeops=True, verbose=False))
        
    return run_list



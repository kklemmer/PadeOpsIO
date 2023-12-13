import numpy as np
import os
import re
import glob

from padeopsIO.turbine import Turbine
from padeopsIO.wake_utils import get_xids

def induction_calc(io_prim, io_base, turbine=None):

    if turbine is None:
        turbine = io_prim.turbineArray.turbines[0]

    x = io_prim.xLine + turbine.xloc
    y = io_prim.yLine + turbine.yloc
    z = io_prim.zLine + turbine.zloc

    # calculate fwidth
    h = np.sqrt(io_prim.dx**2 + io_prim.dy**2 + io_prim.dz**2)
    fwidth = 1.5*h

    # get the turbine kernel (currently based off of ADM type 5)
    turbine.get_kernel(x,y,z,fwidth=fwidth, overwrite=True)

    # print(np.where(turbine.kernel > 0))

    # calculate u at the disk from the primary 
    u_disk = np.sum(np.multiply(turbine.kernel,io_prim.budget['ubar']))

    # calculate u_ing at the disk from the precursor 
    u_inf_disk = np.sum(np.multiply(turbine.kernel,io_base.budget['ubar']))

    # calculate the induction
    a = 1 - u_disk/u_inf_disk

    return a
    
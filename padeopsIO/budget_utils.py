# Additional budget functions 

import numpy as np
import os
import padeopsIO

def advection(uj, arr, dx, dy, dz):
    """
    returns advection of arr by uj 
    i.e. u * d/dx(arr) + v * d/dy(arr) + w * d/dz(arr)
    """

    tmp_grad = np.gradient(arr, dx, dy, dz)

    adv = np.sum(np.multiply(uj, tmp_grad), axis=0, keepdims=False)

    return adv

def construct_delta_uiuj(defIO):
    """
    constructs and returns delta ui'uj'
    """

    delta_uiuj = np.zeros([defIO.nx, defIO.ny, defIO.nz, 3, 3])
    delta_uiuj[...,0,0] = defIO.budget['delta_uu']
    delta_uiuj[...,0,1] = defIO.budget['delta_uv']
    delta_uiuj[...,0,2] = defIO.budget['delta_uw']
    delta_uiuj[...,1,0] = defIO.budget['delta_uv']
    delta_uiuj[...,1,1] = defIO.budget['delta_vv']
    delta_uiuj[...,1,2] = defIO.budget['delta_vw']
    delta_uiuj[...,2,0] = defIO.budget['delta_uw']
    delta_uiuj[...,2,1] = defIO.budget['delta_vw']
    delta_uiuj[...,2,2] = defIO.budget['delta_ww']

    return delta_uiuj

def construct_delta_ui_base_uj(defIO):
    """
    constructs and return delta ui' base uj'
    """

    delta_ui_base_uj = np.zeros([defIO.nx, defIO.ny, defIO.nz, 3, 3])
    delta_ui_base_uj[...,0,0] = defIO.budget['delta_u_base_u']
    delta_ui_base_uj[...,0,1] = defIO.budget['delta_u_base_v']
    delta_ui_base_uj[...,0,2] = defIO.budget['delta_u_base_w']
    delta_ui_base_uj[...,1,0] = defIO.budget['base_u_delta_v']
    delta_ui_base_uj[...,1,1] = defIO.budget['delta_v_base_v']
    delta_ui_base_uj[...,1,2] = defIO.budget['delta_v_base_w']
    delta_ui_base_uj[...,2,0] = defIO.budget['base_u_delta_w']
    delta_ui_base_uj[...,2,1] = defIO.budget['base_v_delta_w']
    delta_ui_base_uj[...,2,2] = defIO.budget['delta_w_base_w']

    return delta_ui_base_uj

def construct_uiuj(budgetIO):
    """
    constructs and returns ui'uj' for BudgetIO object
    """    

    uiuj = np.zeros([budgetIO.nx, budgetIO.ny, budgetIO.nz, 3, 3])
    uiuj[...,0,0] = budgetIO.budget['uu']
    uiuj[...,0,1] = budgetIO.budget['uv']
    uiuj[...,0,2] = budgetIO.budget['uw']
    uiuj[...,1,0] = budgetIO.budget['uv']
    uiuj[...,1,1] = budgetIO.budget['vv']
    uiuj[...,1,2] = budgetIO.budget['vw']
    uiuj[...,2,0] = budgetIO.budget['uw']
    uiuj[...,2,1] = budgetIO.budget['vw']
    uiuj[...,2,2] = budgetIO.budget['ww']

    return uiuj

def construct_duidxj(io):
    """
    constructs and returns duidxj for BudgetIO or DeficitIO
    """

    duidxj = np.zeros([io.nx, io.ny, io.nz, 3, 3])
    duidxj[...,0,0] = io.budget['dUdx']
    duidxj[...,0,1] = io.budget['dUdy']
    duidxj[...,0,2] = io.budget['dUdz']
    duidxj[...,1,0] = io.budget['dVdx']
    duidxj[...,1,1] = io.budget['dVdy']
    duidxj[...,1,2] = io.budget['dVdz']
    duidxj[...,2,0] = io.budget['dWdx']
    duidxj[...,2,1] = io.budget['dWdy']
    duidxj[...,2,2] = io.budget['dWdz']

    return duidxj
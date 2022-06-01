import numpy as np
import os
import re
import warnings
import glob

import padeopsIO.budgetIO as pio

class YawIO(pio.BudgetIO): 
    """
    Class that extends BudgetIO which adds helper functions for reading turbine power, yaw, etc. 
    """
    
    def __init__(self, dir_name, **kwargs): 
        """
        Calls the constructor of BudgetIO
        """
        
        super().__init__(dir_name, **kwargs)
        
        if self.associate_nml: 
            self.yaw = self.input_nml['ad_coriolisinput']['yaw']
            self.uInflow = self.input_nml['ad_coriolisinput']['uinflow'] * np.cos(self.yaw*np.pi/180.)
            self.vInflow = self.input_nml['ad_coriolisinput']['uinflow'] * -np.sin(self.yaw*np.pi/180.)
        
        if self.verbose: 
            print("Initialized YawIO object")
        
        
    def read_turb_power(self, tidx=None, turb=1, steady=True): 
        """
        Reads the turbine power from the output file *.pow. 

        tidx (int) : time ID to read turbine power from. Default: calls self.unique_budget_tidx()
        turb (int) : Turbine number. Default 1
        steady (bool) : Averages results if True. If False, returns an array containing the contents of `*.pow`. 
        """
        
        if tidx is None: 
            if self.associate_budget: 
                tidx = self.unique_budget_tidx()
            else: 
                tidx = self.unique_tidx(return_last=True)
        
        fname = self.dir_name + '/Run{:02d}_t{:06d}_turbP{:02}.pow'.format(self.runid, tidx, turb)
        power = np.genfromtxt(fname, dtype=float)

        if steady: 
            return np.mean(power)

        return power  # this is an array


    def read_turb_vel(self, tidx=None, turb=1, steady=True, u=True, v=True): 
        """
        Reads the turbine power from the output file *.pow. 

        tidx (int) : time ID to read turbine power from. Default: calls self.unique_budget_tidx()
        turb (int) : Turbine number. Default 1
        steady (bool) : Averages results if True. If False, returns an array containing the contents of `*.pow`. 
        u, v (bool) : dictates whether to return u, v, or both. Default: u=True, v=True
        """
        
        if tidx is None: 
            if self.associate_budget: 
                tidx = self.unique_budget_tidx()
            else: 
                tidx = self.unique_tidx(return_last=True)
        
        ret = ()
        
        for i, ui in enumerate((u, v)): 
            if ui: 
                if i == 0: 
                    u_string = "U"
                else: 
                    u_string = "V"
                    
                fname = self.dir_name + '/Run{:02d}_t{:06d}_turb{:s}{:02}.vel'.format(self.runid, tidx, u_string, turb)
                uturb = np.genfromtxt(fname, dtype=float)

                if steady: 
                    ret += (np.mean(uturb), )

                else: 
                    ret += (uturb, )  # this is an array
                    
        if len(ret) == 0: 
            raise ValueError("u or v must be True, function cannot return nothing")
        if len(ret) == 1: 
            return ret[0]
        else: 
            return ret
        
        


import numpy as np
import os
import re
import warnings
import glob

try: 
    import f90nml
except ImportError: 
    warnings.warn("Could not find package f90nml in system path. ")
    # TODO - handle this somehow
    f90nml = None

import padeopsIO.budgetkey as budgetkey  # defines key pairing
import padeopsIO.inflow as inflow  # interface to retrieve inflow profiles

class TurbineArray(): 
    """
    # TODO fix class docstring
    """
    
    def __init__(self, turb_dir=None, num_turbines=None, ADM_type=2, init_dict=None, verbose=False): 
        """
        Constructor function for a TurbineArray class
        
        Parameters
        ----------
        turb_dir (path-like) : path to a turbine array directory in PadeOps
        num_turbines (int) : optional, number of turbines. Default: number 
            of turbine files in turbine directory
        ADM_type (int) : ADM type. Default: 2. 
        verbose (bool) : additional print statements
        
        Returns
        -------
        TurbineArray class instance
        """

        # this initializes from a dictionary output by todict()
        if init_dict is not None: 
            self.fromdict(init_dict)

            self.verbose = verbose
            if self.verbose: 
                print("TurbineArray: Initialized from", turb_dir)

            return
        
        self.turb_dir = turb_dir
        self.ADM_type = ADM_type
        self.verbose = verbose
        
        # begin reading in turbines
        filenames = os.listdir(turb_dir)
        
        if num_turbines is not None: 
            self.num_turbines = num_turbines
            
            if num_turbines != len(filenames):  
                warnings.warn("Not all turbines in the specified directory will be used")

                if self.verbose: 
                    print("\tRequested {:d} turbines, but found {:d} files.".format(num_turbines, len(filenames)))
        else: 
            self.num_turbines = len(filenames)
            
        
        # for now, each turbine can simply be a dictionary appended to a list
        self.array = []
        if self.verbose: 
            print("Reading turbines from the following files:\n", filenames)
            
        for i, filename in enumerate(filenames): 
            if i >= self.num_turbines: 
                break  # only read in up to num_turbines turbines! 
                
            turb_nml = f90nml.read(os.path.join(turb_dir, filename))
            self.array.append(turb_nml)
            
            if self.verbose: 
                print("\tTurbineArray: added turbine to array form file", filename)
        
        if self.num_turbines == 1: 
            if self.verbose: 
                print("\tAdding convenience variables...")
            
            # initialize some defaults (these may not be correct!)
            self.thickness = 1.5
            self.usecorrection = False
            self.filterwidth = 0.5
            
            # make the variables more accessible
            for key in self.array[0]['actuator_disk'].keys(): 
                self.__dict__[key] = self.array[0]['actuator_disk'][key]
        
        if self.verbose: 
            print("TurbineArray: Initialized from", turb_dir)
    
    def fromdict(self, init_dict): 
        """
        Converts a dictionary object given by todict() back into a TurbineArray object. 
        """
        for key in init_dict.keys(): 
            self.__dict__[key] = init_dict[key]

    
    def todict(self): 
        """
        Converts self.__dict__ into a dictionary with no namelists. 
        """

        return self.__dict__.copy()
        
        

if __name__ == "__main__": 
    """
    TODO - add unit tests to class
    """
    print("turbineArray.py: No unit tests included yet. ")
    

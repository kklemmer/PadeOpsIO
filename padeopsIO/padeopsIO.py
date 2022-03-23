from tkinter import W
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

import budgetkey  # defines key pairing
import inflow  # interface to retrieve inflow profiles


class BudgetIO(): 
    """
    # TODO fix class docstring
    """

    key = budgetkey.get_key()
    

    def __init__(self, dir_name, **kwargs):  # outputdir_name, runid=1, tidx=0, Lx=256, Ly=256, Lz=128): 
        """
        Class initiator. Creates different instance variables depending on the keyword arguments given. 
        
        Every instance needs a directory name. If this object is reading information from output files dumped
        by PadeOps, then this is the directory where those files are stored. This object may also read information
        from a local subset of saved data.  
        
        The BudgetIO class will read source files from PadeOps if the following keyword arguments are given: 
            runid (int)
            Lx (int)
            Ly (int)
            Lz (int)
            tidx (int) - optional, default is zero. 

        Alternatively, BudgetIO will try to initialize from source files if kwarg `padeops` is given. 
            
        If not all of those keyword arguments are present, then the directory name will (attempt to) read from 
        budgets of saved .npz files. 

        Regardless of the method of reading budget files, __init__ will initialize the following fields: 
        RUN INFORMATION: 
            filename, dir_name, 
        DOMAIN VARIABLES: 
            Lx, Ly, Lz, nx, ny, nz, dx, dy, dz, xLine, yLine, zLine, 
        TURBINE VARIABLES: 
            nTurb, 
        PHYSICS: 
            Re, Ro, 
        BUDGET VARIABLES: 
            last_tidx, last_n, 
        
        """
        
        # print statements? default False
        if 'verbose' in kwargs and kwargs['verbose']: 
            self.verbose = True 
        else: 
            self.verbose = False
        
        self.dir_name = dir_name
        
        # all files associated with this case will begin with <filename>
        if not 'filename' in kwargs: 
            # defaults to the directory name, split on non-word characters
            dir_list = re.split('\W+', dir_name)
            # pick the last non-empty string
            self.filename = next(s for s in reversed(dir_list) if s)
        else: 
            self.filename = kwargs['filename']

        self.filename_budgets = self.filename + '_budgets.npz'  # standardize this
        
        # ========== Associate files ==========

        # if we are given the required keywords, try to initialize from PadeOps source files
        self.associate_padeops = False
        self.associate_npz = False
        self.associate_nml = False
        self.associate_budget = False
        self.associate_grid = False

        if all(x in kwargs for x in ['runid', 'Lx', 'Ly', 'Lz']) or ('padeops' in kwargs): 
            try: 
                self._init_padeops(**kwargs)
                self.associate_padeops = True

                if self.verbose: 
                    print('Initialized BudgetIO at ' + dir_name + ' from PadeOps source files. ')

            except OSError as err: 
                warnings.warn('Attempted to read PadeOps output files, but at least one was missing.')
                print(err)
        
        # if PadeOps fails to associate, try loading files from .npz
        if not self.associate_padeops: 
            self._init_npz(**kwargs)
            self.associate_npz = True

            if self.verbose: 
                print('Initialized BudgetIO at ' + dir_name + ' from .npz files. ')

        self.budget = {}  # empty dictionary

        if 'read_budgets' in kwargs: 
            # if read_budgets passed in as keyword argument, read budgets on initialization
            self.read_budgets(budget_terms=kwargs['read_budgets'])
    

    def _init_padeops(self, **kwargs): 
        """
        Initializes source files to be read from output files in PadeOps. 
        
        Raises OSError if source files cannot be read
        """

        # Begin with keyword arguments - initialize if these were passed in
        if all(x in kwargs for x in ['runid', 'Lx', 'Ly', 'Lz']): 
            self.runid = kwargs['runid']
            self.Lx = kwargs['Lx']
            self.Ly = kwargs['Ly']
            self.Lz = kwargs['Lz']  
                
        # do namelist stuff - this is necessary if Lx, Ly, ... etc were not passed in. 
        try: 
            self._read_inputfile()
        except FileNotFoundError as err:  # TODO fix this
            warnings.warn("BudgetIO: Problem reading input file. ")
            print(err)

        # default time ID for initialization
        if not 'tidx' in kwargs: 
            self.tidx = 0  # tidx defualts to zero
        else: 
            self.tidx = kwargs['tidx']

        # These lines are taken almost verbatim from PadeOpsViz.py

        # read accompanying info file
        # may throw OSError
        info_fname = self.dir_name + '/Run' + '{:02d}'.format(self.runid) + '_info_t' + '{:06d}'.format(self.tidx) + '.out'
        self.info = np.genfromtxt(info_fname, dtype=None)
        self.time = self.info[0]
        
        self.nx = int(self.info[1])
        self.ny = int(self.info[2])
        self.nz = int(self.info[3])

        self._load_grid()

        # object is reading from PadeOps output files directly
        if self.verbose: 
            print('BudgetIO initialized using info files at time:' + '{:.06f}'.format(self.time))

        self.last_tidx = self.last_budget_tidx()  # last tidx in the run with budgets
        self.last_n = self.last_budget_n()  # last tidx with an associated budget

    
    def _read_inputfile(self, **kwargs): 
        """
        Reads the input file (Fortran 90 namelist) associated with the CFD simulation. 

        Dependencies: f90nml, see https://github.com/marshallward/f90nml 
        """
        
        if f90nml is None: 
            warnings.warn('_raed_inputfile(): No namelist reader loaded. ')
            return

        inputfile_ls = []  

        if 'runid' in kwargs.keys(): 
            # if given a runid, search for an inputfile with Run{runid} in the name
            inputfile_ls = glob.glob(self.dir_name + os.sep + 'Run{:02d}*.dat'.format(kwargs['runid']))

            if len(inputfile_ls) == 0 and self.verbose: 
                warnings.warn('_read_inputfile(): runid {:d} requested, \
                    but no matching inputfile found was found.'.format(kwargs['runid']))
                        
        # if no runid given, or if the previous search failed, search all files ending in '*.dat' 
        if len(inputfile_ls) == 0: 
            inputfile_ls = glob.glob(self.dir_name + os.sep + '*.dat')  # for now, just search this. # TODO improve later? 

        # there should only be one input file for each run 
        if len(inputfile_ls) > 1: 
            warnings.warn('_read_inputfile(): Multiple files ending in *.dat found')
            if self.verbose: 
                print("    Found the following files:", inputfile_ls)
        
        if self.verbose: 
            print("Reading namelist file from {}".format(inputfile_ls[0]))
            
        self.input_nml = f90nml.read(inputfile_ls[0])
        
        # aside from reading in the Namelist, which has all of the metadata, also make some
        # parameters more accessible

        # RUN VARIABLES: 
        self.runid = self.input_nml['io']['runid']

        # DOMAIN VARIABLES: 
        self.nx = self.input_nml['input']['nx']
        self.ny = self.input_nml['input']['ny']
        self.nz = self.input_nml['input']['nz']
        self.Lx = self.input_nml['ad_coriolisinput']['Lx']
        self.Ly = self.input_nml['ad_coriolisinput']['Ly']
        self.Lz = self.input_nml['ad_coriolisinput']['Lz']

        if not self.associate_grid: 
            self._load_grid()
        
        # TURBINE VARIABLES: 
        self.nTurb = self.input_nml['windturbines']['num_turbines']

        # PHYSICS: 
        if self.input_nml['physics']['isinviscid']:  # boolean
            self.Re = np.Inf
        else: 
            self.Re = self.input_nml['physics']['Re']

        if self.input_nml['physics']['usecoriolis']: 
            self.Ro = self.input_nml['physics']['Ro']
        else: 
            self.Ro = np.Inf

        self.associate_nml = True  # successfully loaded input file
        
        # BUDGET VARIABLES: 
        self.last_tidx = self.last_budget_tidx()
        self.last_n = self.last_budget_n() 

    
    def _load_grid(self): 
        """
        Creates dx, dy, dz, and xLine, yLine, zLine variables. 
        
        Expects (self.)Lx, Ly, Lz, nx, ny, nz to exist. 
        """

        if self.associate_grid and self.verbose: 
            print("_load_grid(): Grid already exists. ")
            return
        
        self.dx = self.Lx/self.nx
        self.dy = self.Ly/self.ny
        self.dz = self.Lz/self.nz

        # initialize grid
        self.xLine = np.linspace(0,self.Lx-self.dx,self.nx)
        self.yLine = np.linspace(0,self.Ly-self.dy,self.ny)
        self.zLine = np.linspace(self.dz/2,self.Lz-(self.dz/2),self.nz)

        self.associate_grid = True


    def _init_npz(self, **kwargs): 
        """
        Initializes the BudgetIO object by attempting to read .npz files saved from a previous BudgetIO object 
        from write_npz(). 

        Expects target files: 
        One filename including "{filename}_budgets.npz"
        One filename including "_metadata.npz"
        """

        # budget_files = glob.glob(self.dir_name + os.sep + '*_budget.npz')

        # TODO need to initialize metadata fields, establish what budgets exist, etc. 

        if self.verbose: 
            print('  TODO: this did not actually do anything.')
            print('BudgetIO initialized using .npz files.')


    def set_filename(self, filename): 
        """
        Changes the filename associated with this object. 

        Make sure filename_budgets is consistent with __init__()
        """
        self.filename = filename
        self.filename_budgets = filename + "_budgets.npz"
        
        
    def write_npz(self, write_dir=None, budget_terms='default', filename=None, overwrite=False): 
        """
        Saves budgets as .npz files. Each budget receives its own .npz file, with the fourth dimension representing
        the budget term number (minus one, because python indexing starts at zero). 
        
        Budgets are defined in e.g. PadeOps/src/incompressible/budget_time_avg.F90. See budget_key.py
        From a high level: 
            Budget 0: Mean quantities (1st and 2nd order)
            Budget 1: Momentum budget terms
            Budget 2: MKE budget terms
            Budget 3: TKE budget terms
            
        parameters 
        ----------
        write_dir (str) : location to write .npz files. Default: same directory as self.outputdir_name
        budget_terms (dict or str) : dictionary of terms for each budget. Dictionaries should be formatted 
                {<budget #> : [ <list terms>], 
                 <next budget #> : [ <next list terms>]} 
            Alternatively, terms can be a string 'default', which uses the dictionary: 
                {0: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 
                 1: 'all'}
            or budget_terms can also be 'all', which loads all existing budgets and all existing terms for each budget. 
        filename (str) : calls self.set_filename()
        """
        
        if not self.associate_budget: 
            warnings.warn('write_npz(): No budgets associated! ') 
            return 

        # don't unintentionally overwrite files... 
        write_arrs = False  # this variable doesn't actually do anything
        if not os.path.exists(filepath): 
            write_arrs = True

        elif overwrite: 
            warnings.warn("File already exists, overwriting... ")
            write_arrs = True

        else: 
            warnings.warn("Existing files found. Failed to write; try passing overwrite=True to override.")
            return
        
        # declare directory to write to, default to the working directory
        if write_dir is None: 
            write_dir = self.dir_name
            
        # need to parse budget_terms with the key
        key_subset = self._parse_budget_terms(budget_terms)

        # load budgets
        self.read_budgets(key_subset)

        # if `filename` is provided, change this in the object
        # importantly, this needs to be done AFTER reading budgets! 
        if filename is not None: 
            self.set_filename(filename)

        filepath = write_dir + os.sep + self.filename_budgets

        save_arrs = {}
        for key in key_subset: 
            # TODO: We might want to crop the domain of the budgets here to reduce filesize... do that here: 
            save_arrs[key] = self.budget[key]  # only save the requested budgets

        # write npz files! 
        if write_arrs: 
            np.savez(filepath, **save_arrs)
            self._write_metadata()  # TODO
            
            if self.verbose: 
                print("write_npz: Successfully saved the following budgets: ", list(key_subset.keys()))
                print("at " + filepath)
        
        
    def _write_metadata(self): 
        """
        The saved budgets aren't useful on their own unless we also save some information like the mesh
        used in the simulation and (possibly) some other information. That goes here. 
        """
        pass  # TODO


    def _read_metadata(self): 
        """
        Reads the saved .npz written in _write_metadata(). 
        """
        pass  # TODO

    
    def read_budgets(self, budget_terms='default', mmap=None): 
        """
        Accompanying method to write_budgets. Reads budgets saved as .npz files 
        
        Arguments 
        ---------
        budget_terms : list of terms (see ._parse_budget_terms())
        mmap : default None. Sets the memory-map settings in numpy.load(). Expects None, 'r+', 'r', 'w+', 'c'
       """

        # need to parse budget_terms with the key
        key_subset = self._parse_budget_terms(budget_terms)

        # TODO - if a) budgets are already loaded, or 
        #           b) requested budgets do not exist, then we need to remove these from the list.  # this is done :)

        if self.associate_padeops: 
            self._read_budgets_padeops(key_subset)
        
        elif self.associate_npz: 
            self._read_budgets_npz(key_subset, mmap=mmap)
        
        self.associate_budget = True
        if self.verbose and len(key_subset) > 0: 
            print("read_budgets: Successfully loaded budgets. ")
        

    def _read_budgets_padeops(self, key_subset): 
        """
        Uses a method similar to ReadVelocities_Budget() in PadeOpsViz to read and store full-field budget terms. 
        """
        tidx = self.last_tidx 
        n = self.last_n

        # TODO - this needs testing
        for key in key_subset:
            budget, term = BudgetIO.key[key]
            u_fname = self.dir_name + '/Run' + '{:02d}'.format(self.runid) + '_budget' + '{:01d}'.format(budget) + \
                '_term' + '{:02d}'.format(term) + '_t' + '{:06d}'.format(tidx) + '_n' + '{:06d}'.format(n) + '.s3D'
            
            temp = np.fromfile(u_fname, dtype=np.dtype(np.float64), count=-1)
            self.budget[key] = temp.reshape((self.nx,self.ny,self.nz), order='F')  # reshape into a 3D array

        if self.verbose: 
            print('PadeOpsViz loaded the budget fields at time:' + '{:.06f}'.format(tidx))


    def _read_budgets_npz(self, key_subset, mmap=None): 
        """
        Reads budgets written by .write_npz() and loads them into memory
        """

        # load the npz file and keep the requested budget keys
        for key in key_subset: 
            npz = np.load(self.dir_name + os.sep + self.filename_budgets)
            self.budget[key] = npz[key]  

        if self.verbose: 
            print('PadeOpsViz loaded the following budgets from .npz: ', list(key_subset.keys()))


    def _parse_budget_terms(self, budget_terms): 
        """
        Takes a list of budget terms, either keyed in index form (budget #, term #) or in common form (e.g. ['u_bar', 'v_bar'])
        and returns a subset of the `keys` dictionary that matches two together. `keys` dictionary is always keyed in common form. 

        budget_terms can also be a string: 'all', or 'default'. 

        'default' tries to load the following: 
            Budget 0 terms: ubar, vbar, wbar, all Reynolds stresses, and p_bar
            Budget 1 terms: all momentum terms
        'all' checks what budgets exist and tries to load them all. 

        For more information on the bi-directional keys, see budget_key.py
        """

        # add string shortcuts here... # TODO move shortcuts to budgetkey.py? 
        if budget_terms=='default': 
            budget_terms = ['ubar', 'vbar', 'wbar', 
                            'tau11', 'tau12', 'tau13', 'tau22', 'tau23', 'tau33', 
                            'pbar']

        elif budget_terms=='all': 
            budget_terms = self.existing_terms()

        elif type(budget_terms)==str: 
            warnings.warn("keyword argument budget_terms must be either 'default', 'all', or a list.")
            return {}  # empty dictionary

        
        # parse through terms: they are either 1) valid, 2) missing (but valid keys), or 3) invalid (not in BudgetIO.key)

        existing_keys = self.existing_terms()
        existing_tup = [BudgetIO.key[key] for key in existing_keys]  # corresponding associated tuples (#, #)

        valid_keys = [t for t in budget_terms if t in existing_keys]
        missing_keys = [t for t in budget_terms if t not in existing_keys and t in BudgetIO.key]
        invalid_terms = [t for t in budget_terms if t not in BudgetIO.key and t not in BudgetIO.key.inverse]

        valid_tup = [tup for tup in budget_terms if tup in existing_tup]  # existing tuples
        missing_tup = [tup for tup in budget_terms if tup not in existing_tup and tup in BudgetIO.key.inverse]

        # now combine existing valid keys and valid tuples, removing any duplicates

        valid_terms = set(valid_keys + [BudgetIO.key.inverse[tup][0] for tup in valid_tup])  # combine and remove duplicates
        missing_terms = set(missing_keys + [BudgetIO.key.inverse[tup][0] for tup in missing_tup])

        # generate the key
        key_subset = {key: BudgetIO.key[key] for key in valid_terms}

        
        # warn the user if some requested terms did not exist
        if len(key_subset) == 0: 
            warnings.warn('_parse_budget_terms(): No keys being returned; none matched.')

        if len(missing_terms) > 0: 
            warnings.warn('_parse_budget_terms(): Several terms were requested but the following could not be found: \
                {}'.format(missing_terms))

        if len(invalid_terms) > 0: 
            warnings.warn('_parse_budget_terms(): The following budget terms were requested but the following do not exist: \
                {}'.format(invalid_terms))

        return key_subset 


    def _get_inflow(self, offline=False, wInflow=False): 
        """
        Calls the appropriate functions in inflow.py to retrieve the inflow profile for the corresponding flow. 

        Arguments
        ---------
        offline (bool) : if True, uses the target inflow profile prescribed by initialize.F90. Default (False) reads
            the inflow profile from the first x-index of the domain and average over y. 
        wInflow (bool) : if True, returns an array of w-inflow velocities. Default (False) only returns u, v. 
        
        Returns
        -------
        u (array) : [nz x 1] array of u-velocities as a function of height
        v (array) : [nz x 1] array of v-velocities as a function of height
        w (array) : [nz x 1] array of w-velocities as a function of height. Nominally this is all zero. 

        """
        
        # load using InflowParser: 

        if offline: 
            if self.associate_nml: 
                print(dict(self.input_nml['AD_coriolisinput']))
                u, v = inflow.InflowParser.inflow_offline(**dict(self.input_nml['AD_coriolisinput']), zLine=self.zLine)
            
            # reading from the budgets
            else: 
                warnings.warn('_get_inflow: Requested offline inflow, but namelist not associated. Trying online.')
                u, v = inflow.InflowParser.inflow_budgets(self)
                
        else: 
            u, v = inflow.InflowParser.inflow_budgets(self) 

        # return requested information 

        if wInflow: 
            # If this is not nominally zero, then this will need to be fixed. 
            w = np.zeros(self.zLine.shape)
            return u, v, w

        else: 
            return u, v

    
    def calc_wake(self, offline=False, wInflow=False, overwrite=False):
        """
        Computes the wake deficit by subtracting the target inflow from the flow field. 

        Arguments
        ---------
        see _get_inflow()
        overwrite (bool) : re-computes wakes if they are already read in. 

        Returns
        -------
        None (updates self.budget[] with 'uwake' and 'vwake' keys)

        """ 

        target_terms = ['uwake', 'vwake']
        if wInflow: 
            target_terms.append('wwake')

        # check to see if the terms exist already
        if all(t in self.budget.keys() for t in target_terms): 
            if not overwrite: 
                warnings.warn('Wake terms already computed, returning. To compute anyway, use keyword overwrite=True')
                return

        req_terms = ['ubar', 'vbar']  # required budget terms
        if wInflow: 
            req_terms.append('wbar') # might also need w
        
        # Need mean velocity fields to be loaded
        if not all(t in self.budget.keys() for t in req_terms): 
            self.read_budgets(budget_terms=req_terms)

        # retrieve inflow
        if wInflow: 
            u, v, w = self._get_inflow(offline=offline, wInflow=wInflow)
            self.budget['wwake'] = self.budget['wbar'] - w[np.newaxis, np.newaxis, :]
        
        else: 
            u, v = self._get_inflow(offline=offline, wInflow=wInflow)

        self.budget['uwake'] = u[np.newaxis, np.newaxis, :] - self.budget['ubar']
        self.budget['vwake'] = v[np.newaxis, np.newaxis, :] - self.budget['vbar']

        if self.verbose: 
            print("calc_wake(): Computed wake velocities. ")


    def unique_tidx(self): 
        """
        Pulls all the unique tidx values from a directory. 
        """

        if not self.associate_padeops:
            pass  # TODO - is this lost information? is it useful information? 
            
        # retrieves filenames and parses unique integers, returns an array of unique integers
        filenames = os.listdir(self.dir_name)
        runid = self.runid
        
        # searches for the formatting *_t(\d+)* in all filenames
        t_list = [int(re.findall('Run{:02d}.*_t(\d+).*'.format(runid), name)[0]) 
                  for name in filenames if re.findall('Run{:02d}.*_t(\d+).*'.format(runid), name)]
        return np.unique(np.array(t_list))

    
    def last_budget_tidx(self): 
        """
        Pulls all the unique tidx values from a directory. 
        """

        # TODO: fix for .npz

        # retrieves filenames and parses unique integers, returns an array of unique integers
        filenames = os.listdir(self.dir_name)
        runid = self.runid
        
        # searches for the formatting *_t(\d+)* in budget filenames
        t_list = [int(re.findall('Run{:02d}.*budget.*_t(\d+).*'.format(runid), name)[0]) 
                  for name in filenames if re.findall('Run{:02d}.*budget.*_t(\d+).*'.format(runid), name)]
        return np.max(t_list)

    
    def last_budget_n(self): 
        """
        Pulls all unique n from budget terms in a directory and returns the largest value. 
        """

        # TODO: fix for .npz

        filenames = os.listdir(self.dir_name)
        runid = self.runid

        # capturing *_n(\d+)* in filenames
        t_list = [int(re.findall('Run{:02d}.*_n(\d+).*'.format(runid), name)[0]) 
                  for name in filenames if re.findall('Run{:02d}.*_n(\d+).*'.format(runid), name)]
        return np.max(t_list)
    
    
    def existing_budgets(self): 
        """
        Checks file names for which budgets were output.  
        """
        filenames = os.listdir(self.dir_name)

        if self.associate_npz: 
            filename = self.dir_name + os.sep + self.filename_budgets
            with np.load(filename) as npz: 
                t_list = npz.files  # load all the budget filenames
            
            budget_list = [BudgetIO.key[t][0] for t in t_list]

        elif self.associate_padeops: 
            runid = self.runid
            # capturing *_budget(\d+)* in filenames
            budget_list = [int(re.findall('Run{:02d}.*_budget(\d+).*'.format(runid), name)[0]) 
                            for name in filenames if re.findall('Run{:02d}.*_budget(\d+).*'.format(runid), name)]

        else: 
            warnings.warn('existing_budgets(): No associated budget files found. ')
            return []

        return list(np.unique(budget_list))
    
    
    def existing_terms(self, budget=None): 
        """
        Checks file names for a particular budget and returns a list of all the existing terms.  

        Arguments 
        ---------
        budget (integer) : optional, default None. If provided, searches a particular budget for existing terms. 
            Otherwise, will search for all existing terms. `budget` can also be a list of integers. 
            Budget 0: mean statistics
            Budget 1: momentum
            Budget 2: MKE
            Budget 3: TKE
            Budget 5: Wake deficit

        Returns
        -------
        t_list (list) : list of tuples of budgets found

        """

        t_list = []
        
        # if no budget is given, look through all saved budgets
        if budget is None: 
            budget_list = self.existing_budgets()
        
        else: 
            # convert to list if integer is given
            if type(budget) != list: 
                budget_list = [budget]
            else: 
                budget_list = budget

        # find budgets matching .npz convention in write_npz()
        if self.associate_npz: 
            filename = self.dir_name + os.sep + self.filename_budgets
            with np.load(filename) as npz: 
                all_terms = npz.files

            if budget is None:  # i.e. requesting all budgets
                return all_terms  # we can stop here without sorting through each budget
            
            tup_list = [BudgetIO.key[t] for t in all_terms]  # list of associated tuples
            t_list = []  # this is the list to be built and returned

            for b in budget_list: 
                t_list += [tup for tup in tup_list if tup[0] == b]

        # find budgets by name matching with PadeOps output conventions
        elif self.associate_padeops: 

            filenames = os.listdir(self.dir_name)
            runid = self.runid
            
            tup_list = []
            # loop through budgets
            for b in budget_list: 
                # capturing *_term(\d+)* in filenames
                terms = [int(re.findall('Run{:02d}_budget{:01d}_term(\d+).*'.format(runid, b), name)[0]) 
                        for name in filenames if 
                        re.findall('Run{:02d}_budget{:01d}_term(\d+).*'.format(runid, b), name)]
                tup_list += [((b, term)) for term in set(terms)]  # these are all tuples
            
            # convert tuples to keys
            t_list = [BudgetIO.key.inverse[key][0] for key in tup_list]
        
        else: 
            warnings.warn('existing_terms(): No terms found for budget ' + str(budget))
            return []
        
        return t_list


if __name__ == "__main__": 
    """
    TODO - add unit tests to class
    """
    print("padeopsIO: No unit tests included yet. ")

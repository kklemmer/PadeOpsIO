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

import budgetkey  # defines key pairing


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
        
        The PadeOpsIO class will create a PadeOpsViz instance if the following keyword arguments are given: 
            runid (int)
            Lx (int)
            Ly (int)
            Lz (int)
            tidx (int) - optional, default is zero. 
            
        If not all of those keyword arguments are present, then the directory name will (attempt to) read from 
        budgets of saved .npz files. 

        Regardless of the method of reading budget files, __init__ will initialize the following fields: 
        RUN INFORMATION: 
            filename, dir_name, 
        DOMAIN VARIABLES: 
            Lx, Ly, Lz, nx, ny, nz, dx, dy, dz, xline, yline, zline, 
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
        if all(x in kwargs for x in ['runid', 'Lx', 'Ly', 'Lz'] ): 
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

        self.associate_budget = False
        self.budget = {}  # empty dictionary
    

    def _init_padeops(self, **kwargs): 
        """
        Initializes source files to be read from output files in PadeOps. 
        
        Raises OSError if source files cannot be read
        """

        # Begin with keyword arguments
        self.runid = kwargs['runid']
        self.Lx = kwargs['Lx']
        self.Ly = kwargs['Ly']
        self.Lz = kwargs['Lz']

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

        self.dx = self.Lx/self.nx
        self.dy = self.Ly/self.ny
        self.dz = self.Lz/self.nz

        # initialize grid
        self.xLine = np.linspace(0,self.Lx-self.dx,self.nx)
        self.yLine = np.linspace(0,self.Ly-self.dy,self.ny)
        self.zLine = np.linspace(self.dz/2,self.Lz-(self.dz/2),self.nz)

        # object is reading from PadeOps output files directly
        if self.verbose: 
            print('BudgetIO initialized using info files at time:' + '{:.06f}'.format(self.time))
        
        # do namelist stuff
        try: 
            self._read_inputfile()
        except Exception as err:  # TODO fix this
            warnings.warn("BudgetIO: Problem reading input file. ")
            print(err)

        self.last_tidx = self.last_budget_tidx()  # last tidx in the run with budgets
        self.last_n = self.last_budget_n()  # last tidx with an associated budget

    
    def _read_inputfile(self, **kwargs): 
        """
        Reads the input file (Fortran 90 namelist) associated with the CFD simulation. 

        Dependencies: f90nml, see https://github.com/marshallward/f90nml 
        """
        pass  # TODO


    def _init_npz(self, **kwargs): 
        """
        Initializes the BudgetIO object by attempting to read .npz files saved from a previous BudgetIO object 
        from write_npz(). 

        Does not load budgets from .npz unless kwargs['load_budgets'] is True. 

        Expects target files: 
        At least one filename including "_budget%d"
        One filename including "_metadata.npz"
        """

        budget_files = glob.glob(self.dir_name + os.sep + '*_budget*.npz')

        # TODO need to initialize metadata fields, establish what budgets exist, etc. 

        if self.verbose: 
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

        # don't unintentionally overwrite files... 
        write_arrs = False  # there is probably a better way to do this
        if not os.path.exists(filepath): 
            write_arrs = True

        elif overwrite: 
            warnings.warn("File already exists, overwriting... ")
            write_arrs = True

        else: 
            warnings.warn("Existing files found. Failed to write; try passing overwrite=True to override.")

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
        budget_terms : list of terms (see ._parse_terms())
        mmap : default None. Sets the memory-map settings in numpy.load(). Expects None, 'r+', 'r', 'w+', 'c'
       """

        # need to parse budget_terms with the key
        key_subset = self._parse_budget_terms(budget_terms)

        # TODO - if a) budgets are already loaded, or 
        #           b) requested budgets do not exist, then we need to remove these from the list. 

        if self.associate_padeops: 
            self._read_budgets_padeops(key_subset)
        
        elif self.associate_npz: 
            self._read_budgets_npz(key_subset, mmap=mmap)
        
        self.associate_budget = True
        if self.verbose: 
            print("read_budget: Successfully loaded budgets. ")
        

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
            self.budget[key] = npz[key]  # WARNING: this key might not exist! 

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
            return  # don't return anything

        # figure out if budget_terms contains common form or index form
        if all(term in BudgetIO.key for term in budget_terms): 
            key_subset = {key: BudgetIO.key[key] for key in budget_terms}

        elif all(term in BudgetIO.key.inverse for term in budget_terms):
            key_subset = {BudgetIO.key.inverse[key][0]: key for key in budget_terms}

        else: 
            warnings.warn('Requested invalid budget terms.')
            if self.verbose: 
                print("Could not find the following terms: ", 
                [term for term in budget_terms if term not in BudgetIO.key])
                print("Or, if searching key.inverse, the following terms could not be found: ", 
                [term for term in budget_terms if term not in BudgetIO.key.inverse])

        # TODO: Throw an exception if requested budgets do not exist!! 

        return key_subset 


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
            # capturing *_budget(\d+).*
            budget_list = [int(re.findall('.*_budget(\d+).*', name)[0]) for name in filenames
                            if re.findall('.*_budget(\d+).*', name)]

        elif self.associate_padeops: 
            runid = self.runid
            # capturing *_budget(\d+)* in filenames
            budget_list = [int(re.findall('Run{:02d}.*_budget(\d+).*'.format(runid), name)[0]) 
                            for name in filenames if re.findall('Run{:02d}.*_budget(\d+).*'.format(runid), name)]

        else: 
            warnings.warn('existing_budgets(): No budgets found. ')
            return []

        return np.unique(np.array(budget_list))
    
    
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
                # this does NOT check to make sure all the budgets request actually exist # TODO

        # find budgets matching .npz convention in write_npz()
        if self.associate_npz: 
            filename = self.dir_name + os.sep + self.filename_budgets
            npz = np.load(filename)
            t_list = npz.files

        # find budgets by name matching with PadeOps output conventions
        elif self.associate_padeops: 

            filenames = os.listdir(self.dir_name)
            runid = self.runid

            # loop through budgets
            for budget in budget_list: 
                # capturing *_term(\d+)* in filenames
                terms = [int(re.findall('Run{:02d}_budget{:01d}_term(\d+).*'.format(runid, budget), name)[0]) 
                        for name in filenames if 
                        re.findall('Run{:02d}_budget{:01d}_term(\d+).*'.format(runid, budget), name)]
                [t_list.append((budget, term)) for term in terms]  # append tuples
        
        else: 
            warnings.warn('existing_terms(): No terms found for budget ' + str(budget))
            return []
        
        return t_list


if __name__ == "__main__": 
    """
    TODO - add unit tests to class
    """
    print("padeopsIO: No unit tests included yet. ")

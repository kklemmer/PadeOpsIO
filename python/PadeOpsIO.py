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

import budget_key  # defines key pairing


class BudgetIO(): 
    """
    # TODO fix class docstring
    """

    key = budget_key.get_key()
    
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
            self.filename = dir_name.split(os.sep)[-1]
        
        # if we are given the following keywords, try to 
        self.associate_padeops = False
        self.associate_npz = False
        if all(x in kwargs for x in ['runid', 'Lx', 'Ly', 'Lz'] ): 
            try: 
                self._init_padeops(**kwargs)
                self.associate_padeops = True
            except OSError as err: 
                warnings.warn('Attempted to read PadeOps output files, but at least one was missing.')
                print(err)

            # self.viz = PadeOpsViz(dir_name, kwargs['runid'], tidx, kwargs['Lx'], kwargs['Ly'], kwargs['Lz'])
        
        if not self.associate_padeops: 
            self._init_npz(**kwargs)
            self.associate_npz = True
    

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
        except Exception as err: 
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

        # 
        if self.verbose: 
            print('BudgetIO initialized using npz files.')


    def set_filename(self, filename): 
        """
        Changes the filename associated with this object. 
        """
        self.filename = filename
        
        
    def write_npz(self, write_dir=None, budget_terms='default', filename=None): 
        """
        Saves budgets as .npz files. Each budget receives its own .npz file, with the fourth dimension representing
        the budget term number (minus one, because python indexing starts at zero). 
        
        Budgets are defined in e.g. PadeOps/src/incompressible/budget_time_avg.F90. See BudgetIO.key (dictionary)
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
        
        if not self.associate_padeops: 
            warnings.warn('Not reading from PadeOps output files, skipping writing files. ') 
            return
        
        # declare directory to write to, default to the working directory
        if write_dir is None: 
            write_dir = self.dir_name
            
        if filename is not None: 
            self.set_filename(filename)

        # parse which budget terms to load: 
        # TODO: Move the budget parsing into its own function. 
        if budget_terms=='default': 
            budget_terms = {0: list(range(1, 11)), 
                            1: self.existing_terms(1)}

        elif budget_terms=='all': 
            budget_terms = {i: self.existing_terms(i) for i in self.existing_budgets()}

        elif type(budget_terms)==str: 
            warnings.warn("keyword argument budget_terms must be either 'default', 'all', or a dictionary.")
            return  # don't write anything

        # now budget_terms is a dictionary; load those budgets and write them to disk: 

        for budget in list(budget_terms):  # iterate thru dictionary
            # load budgets
            self.ReadVelocities_budget(budget, terms=budget_terms[budget])  # TODO fix this 

            # put budgets into dictionary form 
            save_arrs = {str(i): self.viz.budget[:, :, :, budget_terms[0].index(i)] for i in budget_terms[0]}

            # save the loaded budgets as npz
            filepath = write_dir + os.sep + self.filename + '_budget{:1d}.npz'.format(budget)
            np.savez(filepath, **save_arrs)
            print("Successfully saved {:}".format(filepath))
        
        self._write_metadata()  # TODO
        
        
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

    
    def read_budgets(self, budget_terms='default', mmap=None): 
        """
        Accompanying method to write_budgets. Reads budgets saved as .npz files 
        """
        
        if self.associate_padeops: 
            self._read_budgets_padeops(self, budget_terms)
        
        elif self.associated_npz: 
            self._read_budgets_npz(self, budget_terms, mmap=mmap)
        
        if self.verbose: 
            print("read_budget: Loaded budgets. ")
        

    def _read_budgets_padeops(self, budget_terms): 
        """
        Uses a method similar to ReadVelocities_Budget() in PadeOpsViz to read and store full-field budget terms. 

        Arguments 
        ---------
        budget_terms (see ._parse_terms())
        """
        tidx = self.last_tidx 
        n = self.last_n

        # need to parse budget_terms with the key

        self.budget = np.zeros((self.nx, self.ny, self.nz, np.size(terms))) 
        ind = 0
        for i in terms:
            u_fname = self.dir_name + '/Run' + '{:02d}'.format(self.runid) + '_budget' + '{:01d}'.format(budget_number) + \
            '_term' + '{:02d}'.format(i) + '_t' + '{:06d}'.format(tidx) + '_n' + '{:06d}'.format(n) + '.s3D'
            temp = np.fromfile(u_fname,dtype=np.dtype(np.float64),count=-1)
            self.budget[:,:,:,ind] = self.u_normfactor*temp.reshape((self.nx,self.ny,self.nz),order='F')
            ind+=1
        print('PadeOpsViz loaded the budget fields at time:' + '{:.06f}'.format(tidx))


    
    def unique_tidx(self): 
        # There is almost certainly a better way to do this
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
        runid = self.runid

        # capturing *_budget(\d+)* in filenames
        t_list = [int(re.findall('Run{:02d}.*_budget(\d+).*'.format(runid), name)[0]) 
                  for name in filenames if re.findall('Run{:02d}.*_budget(\d+).*'.format(runid), name)]
        return np.unique(np.array(t_list))
    
    
    def existing_terms(self, budget): 
        """
        Checks file names for a particular budget and returns a list of all the existing terms.  
        """
        
        filenames = os.listdir(self.dir_name)
        runid = self.runid

        # capturing *_term(\d+)* in filenames
        t_list = [int(re.findall('Run{:02d}_budget{:01d}_term(\d+).*'.format(runid, budget), name)[0]) 
                  for name in filenames if 
                  re.findall('Run{:02d}_budget{:01d}_term(\d+).*'.format(runid, budget), name)]
        return np.unique(np.array(t_list))

    
    # def ReadVelocities_budget(self, budget_number, terms=None): 
    #     """
    #     Simplified method for reading saved budgets from raw data using the last available time step.   
    #     """
        
    #     # TODO FIX ME

    #     if terms is None: 
    #         terms = self.existing_terms(budget_number)
        
    #     self.ReadVelocities_budget(self.last_tidx, self.last_n, budget_number, terms)
        

    
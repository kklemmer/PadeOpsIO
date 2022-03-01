import numpy as np
import os
import re
import warnings

from PadeOpsViz import PadeOpsViz


class PadeOpsIO(): 
    """
    This is a derived class extending PadeOpsViz.PadeOpsViz that assists reading and writing python-compatible
    output files to PadeOps CFD output. 
    """
    
    def __init__(self, outputdir_name, **kwargs):  # outputdir_name, runid=1, tidx=0, Lx=256, Ly=256, Lz=128): 
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
        """
                
        self.outputdir_name = outputdir_name
        
        # all files associated with this case will begin with <filename>
        if not 'filename' in kwargs: 
            self.filename = outputdir_name.split(os.sep)[-1]
        
        # create a visualization object if supplied sufficient keyword args
        if all(x in kwargs for x in ['runid', 'Lx', 'Ly', 'Lz'] ): 
            
            # default time ID for initialization
            if not 'tidx' in kwargs: 
                tidx = 0
            else: 
                tidx = kwargs['tidx']
                
            self.viz = PadeOpsViz(outputdir_name, kwargs['runid'], tidx, kwargs['Lx'], kwargs['Ly'], kwargs['Lz'])
            self.last_tidx = self.last_budget_tidx()  # last tidx in the run with budgets
            self.last_budget = self.last_n()  # last tidx with an associated budget
        
        else: 
            self.viz = None
                    
        
    def set_filename(self, filename): 
        """
        Changes the filename associated with this object. 
        """
        self.filename = filename
        
        
    def write_budgets(self, write_dir=None, budget_terms='default', filename=None): 
        """
        Saves budgets as .npz files. Each budget receives its own .npz file, with the fourth dimension representing
        the budget term number (minus one, because python indexing starts at zero). 
        
        Budgets are defined in e.g. PadeOps/src/incompressible/budget_time_avg.F90 
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
            or terms can also be 'all', which loads all existing budgets and all existing terms for each budget. 
        filename (str) : calls self.set_filename()
        """
        
        if self.viz is None: 
            warnings.warn('No PadeOpsViz object linked!') 
            return
        
        # declare directory to write to, default to the working directory
        if write_dir is None: 
            write_dir = self.outputdir_name
            
        if filename is not None: 
            self.set_filename(filename)

        # parse which budget terms to load: 
        if budget_terms=='default': 
            budget_terms = {0: list(range(1, 11)), 
                            1: self.existing_terms(1)}

        elif budget_terms=='all': 
            budget_terms = {i: self.existing_terms(i) for i in self.existing_budgets()}

        elif type(budget_terms)==str: 
            warn("keyword argument budget_terms must be either 'default', 'all', or a dictionary.")
            return  # don't write anything

        # now budget_terms is a dictionary; load those budgets and write them to disk: 

        for budget in list(budget_terms):  # iterate thru dictionary
            # load budgets
            self.ReadVelocities_budget(budget, terms=budget_terms[budget])

            # put budgets into dictionary form 
            save_arrs = {str(i): self.viz.budget[:, :, :, budget_terms[0].index(i)] for i in budget_terms[0]}

            # save the loaded budgets as npz
            filepath = write_dir + os.sep + self.filename + '_budget{:1d}.npz'.format(budget)
            np.savez(filepath, **save_arrs)
            print("Successfully saved {:}".format(filepath))
        
        self._write_otherstuff()  # TODO
        
        
    def _write_otherstuff(self): 
        """
        The saved budgets aren't useful on their own unless we also save some information like the mesh
        used in the simulation and (possibly) some other information. That goes here. 
        """
        pass  # TODO

    
    def read_budgets(self, read_dir, budget_terms='default', mmap=None): 
        """
        Accompanying method to write_budgets. Reads budgets saved as .npz files 
        """
        pass
        
        
    
    def unique_tidx(self): 
        # There is almost certainly a better way to do this
        """
        Pulls all the unique tidx values from a directory. 

        Arguments
        ---------
        dir_name (str) : directory to search for files in
        runid (int) : RunID from the input file, default 1
        """
        # retrieves filbenames and parses unique integers, returns an array of unique integers
        filenames = os.listdir(self.outputdir_name)
        runid = self.viz.runid
        
        # searches for the formatting *_t%i*
        t_list = [int(re.findall('Run{:02d}.*_t(\d+).*'.format(runid), name)[0]) 
                  for name in filenames if re.findall('Run{:02d}.*_t(\d+).*'.format(runid), name)]
        return np.unique(np.array(t_list))

    
    def last_budget_tidx(self): 
        """
        Pulls all the unique tidx values from a directory. 

        Arguments
        ---------
        dir_name (str) : directory to search for files in
        runid (int) : RunID from the input file, default 1
        """
        # retrieves filbenames and parses unique integers, returns an array of unique integers
        filenames = os.listdir(self.outputdir_name)
        runid = self.viz.runid
        
        # searches for the formatting *_t%i*
        t_list = [int(re.findall('Run{:02d}.*budget.*_t(\d+).*'.format(runid), name)[0]) 
                  for name in filenames if re.findall('Run{:02d}.*budget.*_t(\d+).*'.format(runid), name)]
        return np.max(t_list)

    
    def last_n(self): 
        """
        Pulls all unique n from budget terms in a directory and returns the largest value. 
        """
        filenames = os.listdir(self.outputdir_name)

        t_list = [int(re.findall('Run{:02d}.*_n(\d+).*'.format(self.viz.runid), name)[0]) 
                  for name in filenames if re.findall('Run{:02d}.*_n(\d+).*'.format(self.viz.runid), name)]
        return np.unique(np.array(t_list))[-1]
    
    
    def existing_budgets(self): 
        """
        Checks file names for which budgets were output.  
        """
        
        filenames = os.listdir(self.outputdir_name)
        runid = self.viz.runid

        t_list = [int(re.findall('Run{:02d}.*_budget(\d+).*'.format(runid), name)[0]) 
                  for name in filenames if re.findall('Run{:02d}.*_budget(\d+).*'.format(runid), name)]
        return np.unique(np.array(t_list))
    
    
    def existing_terms(self, budget): 
        """
        Checks file names for a particular budget and returns a list of all the existing terms.  
        """
        
        filenames = os.listdir(self.outputdir_name)
        runid = self.viz.runid

        t_list = [int(re.findall('Run{:02d}_budget{:01d}_term(\d+).*'.format(runid, budget), name)[0]) 
                  for name in filenames if 
                  re.findall('Run{:02d}_budget{:01d}_term(\d+).*'.format(runid, budget), name)]
        return np.unique(np.array(t_list))

    
    def ReadVelocities_budget(self, budget_number, terms=None): 
        """
        Simplified method for reading saved budgets from raw data using the last available time step.   
        """
        
        if terms is None: 
            terms = self.existing_terms(budget_number)
        
        self.viz.ReadVelocities_budget(self.last_tidx, self.last_budget, budget_number, terms)
        

    
import numpy as np
import os
import re
import warnings
import glob

import padeopsIO.budgetIO as pio
import padeopsIO.deficitkey as deficitkey  # defines key pairing


class DeficitIO(pio.BudgetIO): 

    key = deficitkey.get_key()
    key_xy = deficitkey.get_key_xy()

    """
    Class that extends BudgetIO to read deficit budgets 
    """
    
    def __init__(self, dir_name, **kwargs): 
        """
        Calls the constructor of BudgetIO
        """
        
        super().__init__(dir_name, **kwargs)
                
        if self.verbose: 
            print("Initialized DeficitIO object")

    
    def existing_budgets(self): 
        """
        Checks file names for which budgets were output.  
        """
        filenames = os.listdir(self.dir_name)

        if self.associate_padeops: 
            runid = self.runid
            # capturing *_budget(\d+)* in filenames
            budget_list = [int(re.findall('Run{:02d}.*_deficit_budget(\d+).*'.format(runid), name)[0]) 
                           for name in filenames 
                           if re.findall('Run{:02d}.*_deficit_budget(\d+).*'.format(runid), name)]
            
        else: 
            if self.associate_npz: 
                filename = self.dir_name + os.sep + self.filename_budgets + '.npz'
                with np.load(filename) as npz: 
                    t_list = npz.files  # load all the budget filenames
            if self.associate_mat: 
                filename = self.dir_name + os.sep + self.filename_budgets + '.mat'
                ret = loadmat(filename)
                t_list = [key for key in ret if key[0] != '_']  # ignore `__header__`, etc. 
            
            budget_list = [self.key[t][0] for t in t_list]

        if len(budget_list) == 0: 
            warnings.warn('existing_budgets(): No associated budget files found. ')
        
        if 0 in budget_list: 
            budget_list.append(5)  # wake budgets can be recovered from mean budgets

        return list(np.unique(budget_list))

    def existing_terms(self, budget=None, include_wakes=False): 
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
            Budget 4: Reynolds stress
            Budget 5: Wake deficit
        include_wakes (bool) : Includes wakes in the returned budget terms if True, default False. 

        Returns
        -------
        t_list (list) : list of tuples of budgets found

        """

        t_list = []

        budget4_comp_dict = {11 : 0,
                             22 : 10,
                             33 : 20,
                             13 : 30, 
                             23 : 40}
        
        
        # if no budget is given, look through all saved budgets
        if budget is None: 
            budget_list = self.existing_budgets()
        
        else: 
            # convert to list if integer is given
            if type(budget) != list: 
                budget_list = [budget]
            else: 
                budget_list = budget

        # find budgets by name matching with PadeOps output conventions
        if self.associate_padeops: 

            filenames = os.listdir(self.dir_name)
            runid = self.runid
            
            tup_list = []
            # loop through budgets
            for b in budget_list: 
                # capturing *_term(\d+)* in filenames
                terms = [int(re.findall('Run{:02d}_deficit_budget{:01d}_term(\d+).*'.format(runid, b), name)[0]) 
                        for name in filenames if 
                        re.findall('Run{:02d}_deficit_budget{:01d}_term(\d+).*'.format(runid, b), name)]
                tup_list += [((b, term)) for term in set(terms)]  # these are all tuples

                # reynolds stress budgets
                if b == 4:
                    for component in budget4_comp_dict:
                        terms = [int(re.findall('Run{:02d}_deficit_budget{:01d}_{:01d}_term(\d+).*'.format(runid, b, component), name)[0]) + budget4_comp_dict[component]
                                 for name in filenames if 
                                 re.findall('Run{:02d}_deficit_budget{:01d}_{:01d}_term(\d+).*'.format(runid, b, component), name)]
                        tup_list += [((b, term)) for term in set(terms)]  # these are all tuples
                    
                
                # wake budgets: 
                wake_budgets = (1, 2, 3)
                if include_wakes and b == 5:  
                    terms = [int(re.findall('Run{:02d}_deficit_budget{:01d}_term(\d+).*'.format(runid, 0), name)[0]) 
                            for name in filenames if 
                            re.findall('Run{:02d}_deficit_budget{:01d}_term(\d+).*'.format(runid, 0), name)]  # read from mean budgets

                    tup_list += [((b, term)) for term in wake_budgets if term in terms]
            
            # convert tuples to keys
            t_list = [self.key.inverse[key][0] for key in tup_list]
        # find budgets matching .npz convention in write_npz()
        else: 
            if self.associate_npz: 
                filename = self.dir_name + os.sep + self.filename_budgets + '.npz'
                with np.load(filename) as npz: 
                    all_terms = npz.files
                
            elif self.associate_mat: 
                filename = self.dir_name + os.sep + self.filename_budgets + '.mat'
                ret = loadmat(filename)
                all_terms = [key for key in ret if key[0] != '_']  # ignore `__header__`, etc. 

            else: 
                raise AttributeError('existing_budgets(): How did you get here? ')

            if budget is None:  # i.e. requesting all budgets
                return all_terms  # we can stop here without sorting through each budget
            
            tup_list = [self.key[t] for t in all_terms]  # list of associated tuples
            t_list = []  # this is the list to be built and returned

            for b in budget_list: 
                t_list += [tup for tup in tup_list if tup[0] == b]

        # else: 
        if len(t_list) == 0: 
            warnings.warn('existing_terms(): No terms found for budget ' + str(budget))

        return t_list
    
    def _read_budgets_padeops(self, key_subset, tidx): 
        """
        Uses a method similar to ReadVelocities_Budget() in PadeOpsViz to read and store full-field budget terms. 
        """

        budget4_components = [11, 22, 33, 13, 23]
        
        if tidx is None: 
            tidx = self.budget_tidx
            
        elif tidx not in self.all_budget_tidx: 
            # find the nearest that actually exists
            tidx_arr = np.array(self.all_budget_tidx)
            closest_tidx = tidx_arr[np.argmin(np.abs(tidx_arr-tidx))]
            
            print("Requested budget tidx={:d} could not be found. Using tidx={:d} instead.".format(tidx, closest_tidx))
            tidx = closest_tidx 
            
        # update self.time and self.tidx: 
#         self.tidx = tidx
        
#         info_fname = self.dir_name + '/Run{:02d}_info_t{:06d}.out'.format(self.runid, self.tidx)
#         self.info = np.genfromtxt(info_fname, dtype=None)
#         self.time = self.info[0]

        # these lines are almost verbatim from PadeOpsViz.py
        for key in key_subset:
            budget, term = DeficitIO.key[key]
            if budget==4:
                component = budget4_components[int(np.floor((term-1)/10))]

                if term > 10:
                    term = term % 10
                    if term == 0:
                        term = 10
                searchstr =  self.dir_name + '/Run{:02d}_deficit_budget{:01d}_{:02d}_term{:02d}_t{:06d}_*.s3D'.format(self.runid, budget, component, term, tidx)
                u_fname = glob.glob(searchstr)[0]  

            else:
                searchstr =  self.dir_name + '/Run{:02d}_deficit_budget{:01d}_term{:02d}_t{:06d}_*.s3D'.format(self.runid, budget, term, tidx)
                u_fname = glob.glob(searchstr)[0]  
            
            self.budget_n = int(re.findall('.*_t\d+_n(\d+)', u_fname)[0])  # extract n from string
            self.budget_tidx = tidx
            
            temp = np.fromfile(u_fname, dtype=np.dtype(np.float64), count=-1)
            self.budget[key] = temp.reshape((self.nx,self.ny,self.nz), order='F')  # reshape into a 3D array

        if self.verbose and len(key_subset) > 0: 
            print('PadeOpsViz loaded the budget fields at time:' + '{:.06f}'.format(tidx))

    def unique_budget_tidx(self, return_last=True): 
        """
        Pulls all the unique tidx values from a directory. 
        
        Parameters
        ----------
        return_last (bool) : If False, returns only the largest TIDX associated with budgets. 
            Else, returns an entire list of unique tidx associated with budgets. Default True
        """

        # TODO: fix for .npz

        # retrieves filenames and parses unique integers, returns an array of unique integers
        filenames = os.listdir(self.dir_name)
        runid = self.runid
        
        # searches for the formatting *_t(\d+)* in budget filenames
        t_list = [int(re.findall('Run{:02d}.*deficit_budget.*_t(\d+).*'.format(runid), name)[0]) 
                  for name in filenames 
                  if re.findall('Run{:02d}.*deficit_budget.*_t(\d+).*'.format(runid), name)]
        
        if len(t_list) == 0: 
            raise ValueError("No budget times found. ")
        
        if return_last: 
            return np.max(t_list)
        else: 
            return np.unique(t_list)

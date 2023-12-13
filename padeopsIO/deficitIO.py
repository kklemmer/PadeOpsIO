import numpy as np
import os
import re
import warnings
import glob

import padeopsIO.budgetIO as pio
import padeopsIO.deficitkey as deficitkey  # defines key pairing

from padeopsIO.budget_utils import *


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

    def grad_stress_calc(self, tidx=None, Lref=1):
        """
        Calculates the velocity and reynolds stress gradients 
        """

        self.budget['ddxk_delta_uiuj'] = np.zeros([self.nx, self.ny, self.nz, 3, 3, 3])
        self.budget['ddxk_delta_ui_base_uj'] = np.zeros([self.nx, self.ny, self.nz, 3, 3, 3])

        tmp_delta_uiuj = np.zeros([self.nx, self.ny, self.nz, 3, 3]) 
        tmp_delta_ui_base_uj = np.zeros([self.nx, self.ny, self.nz, 3, 3]) 

        tmp_delta_uiuj[:,:,:,0,0] = self.budget['delta_uu']
        tmp_delta_uiuj[:,:,:,0,1] = self.budget['delta_uv']
        tmp_delta_uiuj[:,:,:,0,2] = self.budget['delta_uw'] 
        tmp_delta_uiuj[:,:,:,1,1] = self.budget['delta_vv']
        tmp_delta_uiuj[:,:,:,1,2] = self.budget['delta_vw']
        tmp_delta_uiuj[:,:,:,2,2] = self.budget['delta_ww']
         
        tmp_delta_ui_base_uj[:,:,:,0,0] = self.budget['delta_u_base_u']
        tmp_delta_ui_base_uj[:,:,:,0,1] = self.budget['delta_u_base_v']
        tmp_delta_ui_base_uj[:,:,:,0,2] = self.budget['delta_u_base_w']
        tmp_delta_ui_base_uj[:,:,:,1,0] = self.budget['base_u_delta_v'] 
        tmp_delta_ui_base_uj[:,:,:,1,1] = self.budget['delta_v_base_v']
        tmp_delta_ui_base_uj[:,:,:,1,2] = self.budget['delta_v_base_w']
        tmp_delta_ui_base_uj[:,:,:,2,0] = self.budget['base_u_delta_w']
        tmp_delta_ui_base_uj[:,:,:,2,1] = self.budget['base_v_delta_w']
        tmp_delta_ui_base_uj[:,:,:,2,2] = self.budget['delta_w_base_w']
        
        for j in range(3):
            for k in range(3):
                print(np.shape(np.gradient(tmp_delta_uiuj[:,:,:,j,k], self.xLine*Lref, self.yLine*Lref, self.zLine*Lref)))
                self.budget['ddxi_delta_uiuj'][:,:,:,:,j,k] = np.transpose(np.gradient(tmp_delta_uiuj[:,:,:,j,k], self.xLine*Lref, self.yLine*Lref, self.zLine*Lref), [1,2,3,0]) 
                self.budget['ddxi_delta_ui_base_uj'][:,:,:,:,j,k] = np.transpose(np.gradient(tmp_delta_ui_base_uj[:,:,:,j,k], self.xLine*Lref, self.yLine*Lref, self.zLine*Lref), [1,2,3,0]) 


        return

    def tke_budget_calc(self, pre, prim):
        """
        Calculates terms in the TKE budget grouped together
        """

        # read in necessary terms for deficit budget
        def_budget_terms = [term for term in self.key if self.key[term][0] == 3]
        self.read_budgets(budget_terms=def_budget_terms)
        def_budget0_terms = ['delta_u', 'delta_v', 'delta_w', 'delta_uu', 'delta_vv', 'delta_ww', 
                                'delta_uv', 'delta_vw', 'delta_uw', 'delta_u_base_u',
                                'delta_u_base_v', 'base_u_delta_v', 'delta_u_base_w', 'base_u_delta_w',
                                'delta_v_base_v', 'delta_v_base_w', 'base_v_delta_w','delta_w_base_w', 
                                'delta_tau11', 'delta_tau12', 'delta_tau13', 'delta_tau22', 'delta_tau23',
                                'delta_tau33']
        self.read_budgets(budget_terms=def_budget0_terms)

        # read in necessary terms for precursor and primary budgets
        budget_terms = [term for term in pre.key if pre.key[term][0] == 3]
        pre.read_budgets(budget_terms=budget_terms)
        prim.read_budgets(budget_terms=budget_terms) 

        budget0_terms = ['ubar', 'vbar', 'wbar', 'uu', 'vv', 'ww', 'uv', 'uw', 'vw', 'tau11', 'tau12', 
                        'tau22', 'tau23', 'tau33', 'tau13']
        pre.read_budgets(budget_terms=budget0_terms)
        prim.read_budgets(budget_terms=budget0_terms)

        # calculate velocity gradients
        pre.grad_calc()
        prim.grad_calc()
        self.grad_calc()

        # define TKE
        pre.budget['TKE'] = 0.5*(pre.budget['uu'] + pre.budget['vv'] + pre.budget['ww'])
        prim.budget['TKE'] = 0.5*(prim.budget['uu'] + prim.budget['vv'] + prim.budget['ww'])

        self.budget['TKE_wake'] = prim.budget['TKE'] - pre.budget['TKE']
        self.budget['delta_TKE'] = 0.5*(self.budget['delta_uu'] + self.budget['delta_vv'] + self.budget['delta_ww'])
        self.budget['delta_ui_base_ui'] = self.budget['delta_u_base_u'] + self.budget['delta_v_base_v'] + self.budget['delta_w_base_w']

        # make stress tensors
        delta_ui_base_uj = construct_delta_ui_base_uj(self)
        # print(delta_ui_base_uj)

        base_uiuj = construct_uiuj(pre)

        # make velocity gradient tensors
        base_duidxj = construct_duidxj(pre)

        delta_duidxj = construct_duidxj(self)


        # calculate TKE wake budget (prim - pre)
        for term in budget_terms:
                self.budget[term + "_wake"] = prim.budget[term] - pre.budget[term]

        # advection terms
        dx = self.dx
        dy = self.dy
        dz = self.dz
        self.budget['mixed_TKE_base_adv'] = advection([pre.budget['ubar'], pre.budget['vbar'], pre.budget['wbar']], self.budget['delta_ui_base_ui'], dx, dy, dz)
        self.budget['mixed_TKE_delta_adv'] = advection([self.budget['delta_u'], self.budget['delta_v'], self.budget['delta_w']], self.budget['delta_ui_base_ui'], dx, dy, dz)

        self.budget['mixed_TKE_adv_delta_base_k'] = advection([self.budget['delta_u'], self.budget['delta_v'], self.budget['delta_w']], pre.budget['TKE'], dx, dy, dz)
        self.budget['mixed_TKE_adv_delta_base'] = self.budget['mixed_TKE_base_adv'] - self.budget['TKE_adv_delta_base']

        self.budget['mixed_TKE_adv'] = self.budget['TKE_adv_wake'] - self.budget['TKE_adv'] - self.budget['TKE_adv_delta_base']

        # production terms
        self.budget['mixed_TKE_prod_base_base_delta'] = -np.sum(np.multiply(base_uiuj, delta_duidxj), axis=(3,4))
        self.budget['mixed_TKE_prod_base_delta_delta'] = -np.sum(np.multiply(np.transpose(delta_ui_base_uj, [0,1,2,4,3]), delta_duidxj), axis=(3,4))
        self.budget['mixed_TKE_prod_base_delta_base'] = -np.sum(np.multiply(np.transpose(delta_ui_base_uj, [0,1,2,4,3]), base_duidxj), axis=(3,4))
        self.budget['mixed_TKE_prod_delta_base_base'] = -np.sum(np.multiply(delta_ui_base_uj, base_duidxj), axis=(3,4))

        self.budget['mixed_TKE_prod'] = self.budget['TKE_shear_production_wake'] - self.budget['TKE_production'] - self.budget['TKE_prod_delta_base']

        # transport terms
        self.budget['mixed_TKE_p_transport'] = self.budget['TKE_p_transport_wake'] - self.budget['TKE_p_transport']
        self.budget['mixed_TKE_SGS_transport'] = self.budget['TKE_SGS_transport_wake'] - self.budget['TKE_SGS_transport']
        self.budget['mixed_TKE_turb_transport'] = self.budget['TKE_turb_transport_wake'] - self.budget['TKE_turb_transport'] \
                                                    - self.budget['TKE_turb_transport_delta_base']
        
        # buoyancy
        self.budget['mixed_TKE_buoyancy'] = self.budget['TKE_buoyancy_wake'] - self.budget['TKE_buoyancy']

        # dissipation
        self.budget['mixed_TKE_dissipation'] = self.budget['TKE_dissipation_wake'] - self.budget['TKE_dissipation']



    def non_dim_vel(self):
        vel_keys = ['delta_u', 'delta_v', 'delta_w']

        Ug = key_search_r(self.input_nml, 'g_geostrophic')

        for key in vel_keys:
            self.budget[key] /= Ug
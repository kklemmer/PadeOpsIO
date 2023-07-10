"""
Module containing a bi-directional has thable for the mapping and inverse-mapping
for budget terms in index-tuple form (budget #, term #) defined in PadeOps and common
form as defined in get_key()

Kirby Heck
March 2022
"""

class bidict(dict):
    """
    Bi-directional dictionary, courtesy of: 
    https://stackoverflow.com/questions/3318625/how-to-implement-an-efficient-bidirectional-hash-table

    Used to look up a mapping and inverse mapping (not necessarily unique) between key parings. 
    """
    def __init__(self, *args, **kwargs):
        super(bidict, self).__init__(*args, **kwargs)
        self.inverse = {}
        for key, value in self.items():
            self.inverse.setdefault(value,[]).append(key) 

    def __setitem__(self, key, value):
        if key in self:
            self.inverse[self[key]].remove(key) 
        super(bidict, self).__setitem__(key, value)
        self.inverse.setdefault(value,[]).append(key)        

    def __delitem__(self, key):
        self.inverse.setdefault(self[key],[]).remove(key)
        if self[key] in self.inverse and not self.inverse[self[key]]: 
            del self.inverse[self[key]]
        super(bidict, self).__delitem__(key)


def key_search_r(nested_dict, key):
    """
    Performs a recursive search for the dictionary key `key` in any of the dictionaries contained 
    inside `nested_dict`. Returns the value of nested_dict[key] at the first match. 
    
    Parameters
    ----------
    nested_dict (dict-like) : dictionary [possibly of dictionaries]
    key (str) : dictionary key to match
    
    Returns
    -------
    nested_dict[key] if successful, None otherwise. 
    """
    
    try: 
        for k in nested_dict.keys(): 
            if k == key: 
                return nested_dict[k]
            else: 
                res = key_search_r(nested_dict[k], key)
                if res is not None: 
                    return res
        
    except AttributeError as e: 
        return
    
    
def get_key(): 
    """
    Returns a bidirectional hash table between colloquial string values for budget terms and 
    a tuple-look up ordered (budget #, term #) as defined in e.g. budget_time_avg.F90. 
    """
    key = {  # BUDGET 0 TERMS: (1st and second order averages, scalars excluded)
        'delta_u': (0, 1), 
        'delta_v': (0, 2), 
        'delta_w': (0, 3), 
        'delta_p': (0, 4), 
        'delta_uu': (0, 5), 
        'delta_uv': (0, 6), 
        'delta_uw': (0, 7), 
        'delta_vv': (0, 8), 
        'delta_vw': (0, 9), 
        'delta_ww': (0, 10), 
        'delta_u_base_u': (0, 11), 
        'delta_u_base_v': (0, 12), 
        'base_u_delta_v': (0, 13), 
        'delta_u_base_w': (0, 14), 
        'base_u_delta_w': (0, 15), 
        'delta_v_base_v': (0, 16), 
        'delta_v_base_w': (0, 17), 
        'base_v_delta_w': (0, 18), 
        'delta_w_base_w': (0, 19),
        'delta_tau11'   : (0, 20),
        'delta_tau12'   : (0, 21),
        'delta_tau13'   : (0, 22), 
        'delta_tau22'   : (0, 23), 
        'delta_tau23'   : (0, 24),
        'delta_tau33'   : (0, 25),
        'delta_T'       : (0, 26),
        'delta_uT'      : (0, 27),
        'delta_vT'      : (0, 28),
        'delta_wT'      : (0, 29),
        'delta_TT'      : (0, 30),
        # BUDGET 1 TERMS: () 
        'xAdv_total': (1,1),
        'xAdv_base_delta_mean': (1, 2),  # x-advection
        'xAdv_delta_delta_mean': (1, 3),  # x-advection
        'xAdv_delta_base_mean': (1, 4),  # x-advection
        'xAdv_base_delta_fluc': (1, 5),  # x-advection
        'xAdv_delta_delta_fluc': (1, 6),  # x-advection
        'xAdv_delta_base_fluc': (1, 7),  # x-advection
        'dpdx': (1, 8),  # x-pressure gradient
        'xSGS': (1, 9),  # x-sub grid stresses
        'xAD': (1, 10),   # x-Actuator disk
        'yAdv_total': (1,11),
        'yAdv_base_delta_mean': (1, 12),  # y-advection
        'yAdv_delta_delta_mean': (1, 13),  # y-advection
        'yAdv_delta_base_mean': (1, 14),  # y-advection
        'yAdv_base_delta_fluc': (1, 15),  # y-advection
        'yAdv_delta_delta_fluc': (1, 16),  # y-advection
        'yAdv_delta_base_fluc': (1, 17),  # y-advection
        'dpdy': (1, 18), 
        'ySGS': (1, 19), 
        'yAD' : (1, 20),
        'zAdv_total': (1,21),
        'zAdv_base_delta_mean': (1, 22),  # z-advection
        'zAdv_delta_delta_mean': (1, 23),  # z-advection
        'zAdv_delta_base_mean': (1, 24),  # z-advection
        'zAdv_base_delta_fluc': (1, 25),  # z-advection
        'zAdv_delta_delta_fluc': (1, 26),  # z-advection
        'zAdv_delta_base_fluc': (1, 27),  # z-advection
        'dpdz': (1, 28), 
        'zSGS': (1, 29),
        'zB': (1, 30),
        'xCor': (1, 31), # x-coriolis
        'xGeo': (1, 32), # x-geostrophic pressure grad. 
        'yCor': (1, 33), 
        'yGeo': (1, 34),
        # BUGDET TERMS 2
        'MKE_TKE_loss': (2, 1), 
        'MKE_adv': (2, 2), 
        'MKE_turb_transport': (2, 3), 
        'MKE_p_transport': (2, 4), 
        'MKE_SGS_transport': (2, 5), 
        'MKE_dissipation': (2, 6), 
        'MKE_AD': (2, 7), 
        'MKE_geostrophic': (2, 8), 
        'MKE_coriolis': (2, 9),
        'MKE_buoyancy': (2,10),
        'MKE_TKE_loss_delta_delta': (2,11),
        'MKE_TKE_loss_base_delta': (2,12),
        'MKE_TKE_loss_delta_base': (2,13),
        'MKE_turb_transport_delta_delta': (2,14),
        'MKE_turb_transport_base_delta': (2,15),
        'MKE_turb_transport_delta_base': (2,16),
        'MKE_adv_base_delta': (2,17),
        'MKE_adv_delta_delta': (2,18),
        'MKE_adv_delta_base': (2,19)
        }
    
    return bidict(key)


def get_key_xy(): 
    """
    Returns a bidirectional hash table between colloquial string values for budget terms and 
    a tuple-look up ordered (budget #, term #) as defined in e.g. budget_time_avg.F90. 
    """
    key = {  # BUDGET 0 TERMS: (1st and second order averages, scalars excluded)
        'ubar': (0, 1), 
        'vbar': (0, 2), 
        'Tbar': (0, 3), 
        'uu': (0, 4), 
        'uv': (0, 5), 
        'uw': (0, 6), 
        'vv': (0, 7), 
        'vw': (0, 8), 
        'ww': (0, 9), 
        'uT': (0, 10), 
        'vT': (0, 11), 
        'wT': (0, 12), 
        'TT': (0, 13), 
        'tau13': (0, 14), 
        'tau23': (0, 15), 
        'q3': (0, 16), 
        'pbar': (0, 17), 
        'tau11': (0, 18), 
        'tau12': (0, 19), 
        'tau22': (0, 20), 
        'tau33': (0, 21), 
        # BUDGET 1 TERMS: (momentum)
        'DuDt': (1, 1),  # x-advection
        'xSGS': (1, 2),  # x-pressure gradient
        'xVisc': (1, 3),  # x-sub grid stresses
        'xCor': (1, 4),   # x-Actuator disk
        'xGeo': (1, 5),  
        'xAD': (1, 6), 
        'DvDt': (1, 7), 
        'ySGS': (1, 8), 
        'yVisc': (1, 9), 
        'yCor': (1, 10),
        'yGeo': (1, 11), # x-coriolis
        'DwDt': (1, 12), # x-geostrophic pressure grad. 
        'zSGS': (1, 13), 
        'dpdz': (1, 14), 
        # BUDGET 2 TERMS: (MKE)  TODO - improve the naming keys
        'MKE_TKE_loss': (2, 1), 
        'MKE_adv': (2, 2), 
        'MKE_dissipation': (2, 3), 
        'MKE_SGS_tau_transport': (2, 4), 
        'MKE_AD': (2, 5), 
        'MKE_geostrophic': (2, 6), 
        'MKE_coriolis': (2, 7),
        # BUDGET 3 TERMS: (TKE)
        'TKE_shear_production': (3, 1),
        'TKE_turb_transport': (3, 2),  # this term is inclusive of Dk/Dt
        'TKE_p_transport': (3, 3), 
        'TKE_SGS_transport': (3, 4), 
        'TKE_dissipation': (3, 5), 
        'TKE_AD': (3, 6), 
        'TKE_coriolis': (3, 7), 
        'TKE_buoyancy': (3, 8), 
        # BUDGET 4 TERMS:
        'uu_shear_production': (4, 1), # uu
        'uu_turb_transport': (4, 2),
        'uu_p_strain': (4, 3),
        'uu_p_transport': (4, 4),
        'uu_dissipation': (4, 5),
        'uu_SGS_transport': (4, 6),
        'uu_buoyancy': (4, 7),
        'uu_coriolis': (4, 8),
        'uu_AD': (4, 9),
        'uw_shear_production': (4, 10), # uw
        'uw_turb_transport': (4, 11),
        'uw_p_strain': (4, 12),
        'uw_p_transport': (4, 13),
        'uw_dissipation': (4,14),
        'uw_SGS_transport': (4, 15),
        'uw_buoyancy': (4, 16),
        'uw_coriolis': (4, 17),
        'uw_AD': (4, 18),
        'vw_shear_production': (4, 19), # vw
        'vw_turb_transport': (4, 20),
        'vw_p_strain': (4, 21),
        'vw_p_transport': (4, 22),
        'vw_dissipation': (4, 23),
        'vw_SGS_transport': (4, 24),
        'vw_buoyancy': (4, 25),
        'vw_coriolis': (4, 26),
        'vw_AD': (4, 27),
        'ww_shear_production': (4, 28), #ww
        'ww_turb_transport': (4, 29),
        'ww_p_strain': (4, 30),
        'ww_p_transport': (4, 31),
        'ww_dissipation': (4, 32),
        'ww_SGS_transport': (4, 33),
        'ww_buoyancy': (4, 34),
        'ww_coriolis': (4, 35),
        'ww_AD': (4, 36),
        # BUDGET 5 TERMS: Wake deficit
        'uwake': (5, 1), 
        'vwake': (5, 2), 
        'wwake': (5, 3)
        }
    
    return bidict(key)


def key_labels(): 
    """
    Returns a dictionary that assigns a label to each budget key. 
    """
    key = {  # BUDGET 0 TERMS: (1st and second order averages, scalars excluded)
        'ubar': '$\\bar{u}/U$', 
        'vbar': '$\\bar{v}/U$',
        'wbar': '$\\bar{w}/U$',
        'uu': "$\\overline{u'u'}/U^2$",
        'uv': "$\\overline{u'v'}/U^2$",
        'uw': "$\\overline{u'w'}/U^2$",
        'vv': "$\\overline{v'v'}/U^2$",
        'vw': "$\\overline{v'w'}/U^2$",
        'ww': "$\\overline{w'w'}/U^2$",
        'pbar': "$\\bar{p}$", 
        'tau11': "$\\tau_{11}$", 
        'tau12': "$\\tau_{12}$", 
        'tau13': "$\\tau_{13}$", 
        'tau22': "$\\tau_{22}$", 
        'tau23': "$\\tau_{23}$", 
        'tau33': "$\\tau_{33}$", 
        'pu': "$\\overline{p'u'}$", 
        'pv': "$\\overline{p'v'}$", 
        'pw': "$\\overline{p'w'}$", 
        'uk': "$\\overline{k'u'}$", 
        'vk': "$\\overline{k'v'}$", 
        'wk': "$\\overline{k'w'}$", 
        'ujtau1j': (0, 23), 
        'ujtau2j': (0, 24), 
        'ujtau3j': (0, 25), 
        'Tbar': "$\\bar{T}$", 
        'uT': "$\\overline{T'u'}$", 
        'vT': "$\\overline{T'v'}$", 
        'wT': "$\\overline{T'w'}$", 
        'TT': "$\\overline{T'T'}$", 
        # BUDGET 1 TERMS: (momentum)
        'DuDt': "$Du/Dt$",  # x-advection
        'dpdx': "$\\partial p/\\partial x$",  # x-pressure gradient
        'xSGS': "$x$-SGS",  # x-sub grid stresses
        'xAD': "AD$_x$",   # x-Actuator disk
        'DvDt': "$Dv/Dt$",  
        'dpdy': "$\\partial p/\\partial y$", 
        'ySGS': "$y$-SGS", 
        'DwDt': "$Dw/Dt$", 
        'dpdz': "$\\partial p/\\partial z$", 
        'zSGS': '$z$-SGS',
        'xCor': '$x$-Coriolis', # x-coriolis
        'xGeo': '$x$-Geostrophic', # x-geostrophic pressure grad. 
        'yCor': '$y$-Coriolis', 
        'yGeo': '$y$-Geostrophic', 
        'yAD': "AD$_y$", 
        # BUDGET 2 TERMS: (MKE)  TODO - improve the naming keys
        'MKE_TKE_loss': (2, 1), 
        'MKE_adv': (2, 2), 
        'MKE_tau_transport': (2, 3), 
        'MKE_p_transport': (2, 4), 
        'MKE_SGS_transport': (2, 5), 
        'MKE_dissipation': (2, 6), 
        'MKE_AD': (2, 7), 
        'MKE_geostrophic': (2, 8), 
        'MKE_coriolis': (2, 9),
        # BUDGET 3 TERMS: (TKE)
        'TKE_shear_production': "$\\mathcal{P}_k$",
        'DkDt': "$Dk/Dt$", 
        'TKE_turb_transport': "$\\mathcal{T}_t$", 
        #'TKE_p_strain': (3, 3), 
        'TKE_p_transport': "$\\mathcal{T}_p$", 
        'TKE_SGS_transport': "$\\mathcal{T}_{SGS}$", 
        'TKE_dissipation': "$\\varepsilon$", 
        'TKE_buoyancy': (3, 7), 
        'TKE_coriolis': (3, 8), 
        'TKE_AD': (3, 9), 
        # BUDGET 4 TERMS: TODO
        # BUDGET 5 TERMS: Wake deficit
        'uwake': "$\\Delta u_{wake}$",
        'vwake': "$\\Delta v_{wake}$",
        'wwake': "$\\Delta w_{wake}$", 
        # SPATIAL COORDINATES
        'x': "$x/D$", 
        'y': "$y/D$", 
        'z': "$z/D$", 
        't': "$t$"
        # RANS BUDGETS - TODO
#         'ududx': '$\\bar{u}
        }
    
    return key

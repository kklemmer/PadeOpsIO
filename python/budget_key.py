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


def get_key(): 
    """
    Returns a bidirectional has table between colloquial string values for budget terms and 
    a tuple-look up ordered (budget #, term #) as defined in e.g. budget_time_avg.F90. 
    """
    key = {'u_bar': (0, 1), 
           'v_bar': (0, 2), 
           'w_bar': (0, 3), 
           'uu': (0, 4), 
           'uv': (0, 5), 
           'uw': (0, 6), 
           'vv': (0, 7), 
           'vw': (0, 8), 
           'ww': (0, 9), 
           'p_bar': (0, 10) }
           # TODO - finish term matching
    
    return bidict(key)

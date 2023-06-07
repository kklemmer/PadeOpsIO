import re

def parser(filename): 
    """
    Parses a namelist into nested dictionaries. Each key to the returned dictionary 
    is the name of a namelist in `filename`, and the value of that key is itself a 
    dictionary of namelist variables and associated values. 

    Parameters 
    ----------
    filename (path-like) : path to a namelist text file

    Returns
    -------
    Namelists (nested dictionary)
    """

    Namelists = {}
    Namelist = {}  # active namelist
    Active = None
    with open(filename, 'r') as f: 
        lines = f.readlines()

        for line in lines: 
            res = re.search('&(\S+)', line)
            if res is not None: 
                # print('Namelist name: ', res.group(1))
                Active = res.group(1)  # active list
                continue

            res = re.search('(\S+)\s+=\s+(\S+)', line)
            if res is not None: 
                tmp = re.sub('[dD]', 'e', res.group(2))
                key = res.group(1)
                try: 
                    value = float(tmp)
                except ValueError as e: 
                    value = res.group(2)
                # TODO: Handle boolean variables too

                # print(key, value)
                Namelist[key] = value
                continue

            res = re.search(r'/', line)
            if res is not None: 
                # print('Found end of Namelist')
                if Active is not None: 
                    Namelists[Active] = Namelist.copy()
                    Namelist = {}
                    Active = None

    return Namelists

def writer(filepath, Namelists): 
    """
    Writes a namelist formatted in the same way that is output by nml_parser()

    """
    with open(filepath, 'w') as f: 
        for key in Namelists.keys(): 
            f.write('&' + key + '\n')
            for var in Namelists[key].keys(): 
                value = Namelists[key][var]
                if isfloat(value): 
                    value = '{:.08e}'.format(value)
                f.write('{:''<16} = {:s}\n'.format(var, value))
            f.write('/\n')


def isfloat(value): 
    """
    Checks if the input can be cast into a float. 
    https://www.programiz.com/python-programming/examples/check-string-number
    """
    try: 
        float(value)
        return True
    except ValueError: 
        return False

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
                # tmp = re.sub('[dD]', 'e', res.group(2))
                key = res.group(1)
                value = cast_str_to_X(res.group(2))
                # try: 
                #     value = float(tmp)
                # except ValueError as e: 
                #     value = res.group(2)
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
                tmp = Namelists[key][var]
                value = cast_to_str(tmp)
                f.write('{:''<24} = {:s}\n'.format(var, value))
            f.write('/\n')


def isfloat(value): 
    """
    Checks if the input can be cast into a float. 
    https://www.programiz.com/python-programming/examples/check-string-number

    Deprecated 06/07/2023
    """
    try: 
        float(value)
        return True
    except ValueError: 
        return False


def cast_str_to_X(value): 
    """
    Tries to cast the input value to an integer, then a float. If both, then returns the original string. 
    """

    # try integer casting: 
    try: 
        tmp = int(value)
        return tmp
    except ValueError: 
        pass

    # integer casting failed, try float
    try: 
        tmp1 = re.sub('[dD]', 'e', value)  # Fortran uses d, D as well as e, E for exponential notation
        tmp2 = float(tmp1)
        return tmp2
    except ValueError: 
        return value  # returns original value
    

def cast_to_str(value): 
    """
    Takes a value of unknown type (value) and casts to a fortran namelist-readable string. 
    """

    # note: this is NOT a try/except statement because we want floats (that are also integers)
    # to be parsed as floats, not integers. 
    if type(value) == int: 
        return str(value)

    try: 
        tmp = float(value)
        return '{:.08e}'.format(tmp)
    except ValueError: 
        return value  # return the original string 


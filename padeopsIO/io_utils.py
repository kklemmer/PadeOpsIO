# Additional IO functions (some work-in-progress)
# USAGE: from io_utils import *

import csv
import re
import os
import padeopsIO
from numpy import array as array


def read_list(dir_name): 
    """
    Outputs a list of BudgetIO objects from a directory name by reading the associated
    file `./Runs.csv`
    
    Arguments
    ---------
    dir_name (str) : directory path containing a file Runs.csv
    
    Returns
    -------
    run_list : list of BudgetIO objects
    """
    
    # read runs
    # https://stackoverflow.com/questions/24662571/python-import-csv-to-list
    with open(os.path.join(dir_name, 'Runs.csv')) as f:
        reader = csv.reader(f)
        runs = list(reader)

    run_list = []

    for run in runs: 
        run_list.append(padeopsIO.BudgetIO(run[0], padeops=True, verbose=False))
        
    return run_list


def key_search_r(nested_dict, key):
    """
    Copied from budgetkey.py
    
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
    

def query_logfile(filename, search_terms=['TIDX'], 
                  fsearch=r'({:s}).*\s+([-+]?(\d+(\.\d*)?|\.\d+)([dDeE][-+]?\d+)?)'): 
    """
    Queries the PadeOps output log file for text lines printed out by temporalhook.F90. 
    
    By default, the search looks for TERM followed by any arbitrary characters, then at least 1 
    character of white space followed by a number of format %e (exponential). Casts the resulting
    string into a float. 
    
    Parameters
    ----------
    filename (path) : string or path to an output log file. 
    search_terms (list) : list of search terms. Default: length 1 list ['TIDX']. 
        Search terms are case sensitive. 
    fsearch (string) : string format for regex target. 
    
    Returns
    -------
    ret (dict) : dictionary of {search_term: [list of matched values]}
    """
    
    ret = {key: [] for key in search_terms}
    
    search_str = ''
    # build the regex match
    for term in search_terms: 
        search_str += term + '|'

    search = fsearch.format(search_str[:-1])  # do not include last pipe in search_str

    with open(filename, 'r') as f:
        lines = f.readlines()
        
        for k, line in enumerate(lines): 
            match = re.search(search, line)
            if match is not None: 
                key = match.groups()[0]  # this was the matched keyword
                # append the associated matched value, cast into a float
                ret[key].append(float(match.groups()[1]))
                
    # convert lists to array: 
    return {key: array(ret[key]) for key in ret.keys()}

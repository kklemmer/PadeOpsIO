import os 
import sys
import numpy as np

# sys.path.append(r'./padeops_io/')  # why is this our cwd? 
from padeopsIO import BudgetIO

# testIO = BudgetIO('some_dir', verbose=False)
dir_name = r'/scratch/08445/tg877441/AD_coriolis_shear/new_geostrophic/r1s0v1' 
testIO = BudgetIO(dir_name, Lx=50, Ly=20, Lz=10, runid=1, verbose=True)

budgets = testIO.existing_budgets()
testIO.existing_terms(budgets[0])
testIO.read_budgets(budget_terms='all')
testIO.write_npz(budget_terms='all')

print(testIO.__dict__)

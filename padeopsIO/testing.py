import os 
import sys
import numpy as np

# sys.path.append(r'./padeops_io/')  # why is this our cwd? 
from padeopsIO import BudgetIO

# testIO = BudgetIO('some_dir', verbose=False)
testIO = BudgetIO(r'./test_data', Lx=1, Ly=1, Lz=1, runid=1, verbose=True)

budgets = testIO.existing_budgets()
testIO.existing_terms(budgets[0])
testIO.read_budgets(budget_terms='all')
testIO.write_npz(budget_terms='all')

print(testIO.__dict__)

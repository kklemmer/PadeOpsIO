import os 
import sys
import numpy as np

sys.path.append(r'./python/')  # why is this our cwd? 
from PadeOpsIO import BudgetIO

# testIO = BudgetIO('some_dir', verbose=False)
testIO = BudgetIO('random_dir', Lx=1, Ly=1, Lz=1, runid=1, verbose=True)

print(testIO.__dict__)

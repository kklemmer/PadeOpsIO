import numpy as np
import os
import re
import glob

from padeopsIO.io_utils import key_search_r
from padeopsIO.wake_utils import get_xids

class Turbine(): 
    """
    Constructs a Turbine object with properties from the turbine namelist from the input files. 
    If a grid is provided to the turbine object, then the forcing kernel can also be generated. 
    """
    
    req_vars = ['xloc', 'yloc', 'zloc', 'ct', 'diam']  # required variables 
    opt_vars = ['yaw', 'tilt', 'filterwidth', 'usecorrection']  # optional add'l variables

    
    def __init__(self, nml, n=-1, verbose=False, sort='xloc'): 
        """
        Initialization for the Turbine object, which takes a nested dictionary `nml` (or f90nml object)
        to generate a Turbine object. 
        
        Parameters
        ----------
        nml (nested dictionary) : should include, at minimum, the following: 
            xloc, yloc, zloc, ct (= C_T', as in Calaf, et al. (2010)), diam
        """
        
        self.verbose = verbose
        self.input_nml = nml  # save the whole namelist
                
        for req in Turbine.req_vars: 
            ret = key_search_r(nml, req)
            if ret is not None: 
                self.__dict__[req] = ret  # access these more easily
            else: 
                raise AttributeError('Turbine() input namelist missing argument', req)
                
        for opt in Turbine.opt_vars: 
            ret = key_search_r(nml, opt)
            if ret is not None: 
                self.__dict__[opt] = ret  
        
        if self.verbose: 
            print('Turbine(): Successfully initialized turbine at x={:.3f}, y={:.3f}, z={:.3f}'.format(self.xloc, self.yloc, self.zloc))
            
        # default sort: by xloc, then by yloc
        self.sort_by = sort
        
        # turbine number; not included in the input file 
        self.n = n
        
    
    def set_sort(self, sort): 
        if sort in self.__dict__: 
            self.sort_by = sort
        else: 
            raise ValueError('Turbine.set_sort(): given variable does not exist for the turbine object')
            
            
    def get_correction(self, return_correction=False): 
        """
        Computes the correction factor M as defined by Taylor expansion of Shapiro, et al. (2019)
        
        Parameters
        ----------
        return_correction (bool) : Optional, returns correction factor if True. Default: False
        """
        fwidth = key_search_r(self.input_nml, 'filterwidth')
        use_corr = key_search_r(self.input_nml, 'usecorrection')
        
        if fwidth is None or not use_corr: 
            M = 1. 
        else: 
            M = 1/(1. + self.ct/2 * fwidth/np.sqrt(3*np.pi)/self.diam)
            
        if return_correction: 
            return M
        else: 
            self.M = M
        
    
    def get_kernel(self, x, y, z, ADM_type=5, fwidth=None, buffer_fact=3, return_kernel=False): 
        """
        Compute a 3D forcing kernel based off of the coordinate axes xLine, yLine, zLine given that
        these axes define a coordinate system consistent with the turbine xloc, yloc, zloc. 
        
        Parameters
        ----------
        xLine, yLine, zLine (array-like) : 1D arrays of coordinate axes
        ADM_type (int) : integer for the ADM type, consistent with igrid in PadeOps. 
            Default is 5 (Shapiro, et al. (2019))
        fwidth (float) : filter width (smoothing kernel factor)
        buff_fact (float) : grid partition factor 
        return_kernel (bool) : returns the 3D forcing kernel if True. Default: False
            
        Returns
        -------
        kernel (nx * ny * nz) : forcing kernel, only if return_kernel=True
        """
        
        # begin ADM 5
        if ADM_type == 5: 
            
            if fwidth is None: 
                fwidth = key_search_r(self.input_nml, 'filterwidth')
                if fwidth is None:  # still None, raise error
                    raise ValueError("Turbine.get_kernel(): no filterwidth found in the turbine input file")
            
            # get control points: 
            xcs, ycs, zcs = self._get_ctrl_pts(x, y, z)  
            C1 = (6/np.pi/fwidth**2)**(1.5)  # normalizing constant

            kernel = np.zeros((len(x), len(y), len(z)))
            
            # partition grid to a region local to the ADM
            buff = buffer_fact*fwidth
            xmin = min(xcs) - buff
            xmax = max(xcs) + buff
            ymin = min(ycs) - buff
            ymax = max(ycs) + buff
            zmin = min(zcs) - buff
            zmax = max(zcs) + buff

            xids, yids, zids = get_xids(x_ax=x, y_ax=y, z_ax=z, 
                                        x=[xmin, xmax], y=[ymin, ymax], z=[zmin, zmax], 
                                        return_slice=True)

            # build the kernel only on the partitioned grid subset around the ADM (with some buffer)
            X, Y, Z = np.meshgrid(x[xids], y[yids], z[zids], indexing='ij')

            # built the kernel with the Greens function (can be slow!): 
            for xc, yc, zc in zip(xcs, ycs, zcs): 
                kernel[xids, yids, zids] += C1*np.exp(-6./fwidth**2 * ((X-xc)**2 + (Y-yc)**2 + (Z-zc)**2))

            kernel[kernel < 1e-10] = 0  # set these identically to zero, mirroring PadeOps implementation
            
            dx = x[1]-x[0]; dy = y[1]-y[0]; dz = z[1]-z[0]
            kernel /= np.sum(kernel)*dx*dy*dz  # normalize such that the kernel integrates to 1
            # note that np.sum(kernel)*dx*dy*dz ≈ number of control points 
            
            if return_kernel:  # return if requested
                return kernel
            else: 
                self.kernel = kernel
                if self.verbose: 
                    print('Turbine.get_kernel(): computed kernel function')
        
        else: 
            raise ValueError("Turbine.get_kernel(): No other ADM types are currently set up")
        
    
    def _get_ctrl_pts(self, x, y, z): 
        """
        Helper function to the 3D integration of get_kernel for an unyawed, untilted ADM Type 5 (Assumes
        the disk normal vector is in the x-direction). 
        Control points are centered at the ADM location and spaced equally to the grid spacing. 
        """
        dx = x[1]-x[0]
        dy = y[1]-y[0]
        dz = z[1]-z[0]
        
        R = self.diam/2
        y_per_R = np.ceil(R/dy)
        z_per_R = np.ceil(R/dz)  # points in y, z for the temporary grid
        
        Y, Z = np.meshgrid(
            np.arange(-y_per_R, y_per_R+1)*dy, 
            np.arange(-z_per_R, z_per_R+1)*dz, indexing='ij'
        )
        mask = np.ravel((Y**2 + Z**2) < R**2)  # unravel the mask matrix
        
        yravel = np.ravel(Y)
        zravel = np.ravel(Z)  # unravel the 2D grid matrices Y, Z
        
        xctrl = np.zeros(sum(mask)) + self.xloc
        yctrl = yravel[mask] + self.yloc
        zctrl = zravel[mask] + self.zloc
        
        return xctrl, yctrl, zctrl
        
        
    def __lt__(self, other): 
        """
        Allows the turbine array to be sorted, uses self.sort_by to evaluate the sort order. 
        """
        if self.sort_by == 'xloc': 
            if self.xloc == other.xloc: 
                return self.yloc < other.yloc  # tie-breaker
            else: 
                return self.xloc < other.xloc
            
        else: 
            return self.__dict__[self.sort_by] < other.__dict__[self.sort_by]
        
    
    def __str__(self): 
        return "Turbine object at x={:.3f}, y={:.3f}, z={:.3f}".format(self.xloc, self.yloc, self.zloc)
            
        
        
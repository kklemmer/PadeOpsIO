import matplotlib.pyplot as plt
import numpy as np
import os
import sys

import padeopsIO.budgetkey as budgetkey
import padeopsIO.deficitkey as deficitkey


class PlotIO(): 
    """
    Interface class for plotting and saving figures using a BudgetIO object. 
    """

    keylab = budgetkey.key_labels()

    x_lab = '$x$'
    y_lab = '$y$'
    z_lab = '$z$'
    fs = 12
    plt.rcParams.update({'font.size': fs})

    cbar = "coolwarm"
    dpi = 100
    # add other constants here? Note: these are mutable by the user

    momentum_budget_colors = {"Adv_total": "tab:blue", 
                              "dpd": "tab:orange",
                              "SGS": "tab:green",
                              "AD": "tab:red",
                              "Cor": "tab:purple",
                              "Geo": "tab:brown",
                              "B": "tab:gray"}

    budget_colors = {"shear_production": "tab:blue",
                     "buoyancy" : "tab:gray",
                     "AD" : "tab:cyan",
                     "adv" : "tab:orange", #
                     "p_transport" : "tab:purple",  #
                     "SGS_transport" : "tab:brown",   #
                     "turb_transport" : "tab:green", #
                     "p_strain" : "tab:red", #
                     "dissipation" : "tab:pink", #
                     "coriolis" : 'tab:olive'} #

    def set_fontsize(fontsize): 
        """
        Sets the font size
        """
        PlotIO.fs = fontsize
        plt.rcParams.update({'font.size': fontsize})


    def plot_xy(self, io, z, xlim=None, ylim=None, budget_terms=None, field_terms=None, 
                ax=None, xlabel=None, ylabel=None, tidx=None, xscale=None,**kwargs): 
        """
        Plots a slice in the xy-plane. 
        
        If `budget_terms` is a list, then plot_xy will recursively call itself for each term individually
        and produce figures of each separately. 
        
        Arguments
        ---------
        io (BudgetIO obj) : BudgetIO object linked to padeops data
        budget_terms (list or str) : terms to plot
        z (float) : z-location to pull slices from
        xlim, ylim : slice bounds, see BudgetIO.slice()
        
        Returns 
        -------
        
        """

        if not tidx:
            try:
                tidx = io.last_tidx
            except:
                print(" ")
            try:
                tidx = io.last_budget_tidx
            except:
                print(" ")
        if budget_terms:
            slices = io.slice(budget_terms=budget_terms, xlim=xlim, ylim=ylim, zlim=z, tidx=tidx)
            terms_list = budget_terms
        elif field_terms:
            slices = io.slice(field_terms=field_terms, xlim=xlim, ylim=ylim, zlim=z, tidx=tidx)
            terms_list = field_terms

        if not xscale:
            xscale=1

        for term in terms_list: 
            fig, ax = plt.subplots()
            im = ax.imshow(slices[term].T, extent=slices['extent']*xscale, origin='lower',**kwargs)

            common_cbar(fig, im, ax, label=term)

            ax.set_xlabel(PlotIO.x_lab)
            ax.set_ylabel(PlotIO.y_lab)

            plt.show()

        return ax

    
    def plot_xz(self, io, y, xlim=None, zlim=None, budget_terms=None, field_terms=None, 
                ax=None, xlabel=None, ylabel=None, tidx=None, xscale=None, **kwargs): 
        """
        Plots a slice in the xz-plane. 
        
        If `budget_terms` is a list, then plot_xy will recursively call itself for each term individually
        and produce figures of each separately. 
        
        Arguments
        ---------
        io (BudgetIO obj) : BudgetIO object linked to padeops data
        budget_terms (list or str) : terms to plot
        y (float) : y-location to pull slices from
        xlim, zlim : slice bounds, see BudgetIO.slice()
        
        Returns 
        -------
        
        """
        
        if not tidx:
            try:
                tidx = io.last_tidx
            except:
                print(" ")
            try:
                tidx = io.last_budget_tidx
            except:
                print(" ")

        if budget_terms:
            slices = io.slice(budget_terms=budget_terms, xlim=xlim, ylim=y, zlim=zlim, tidx=tidx)
            terms_list = budget_terms
        elif field_terms:
            slices = io.slice(field_terms=field_terms, xlim=xlim, ylim=y, zlim=zlim, tidx=tidx)
            terms_list = field_terms

        if not xscale:
            xscale=1

        for term in terms_list: 
            fig, ax = plt.subplots()
            im = ax.imshow(slices[term].T, extent=slices['extent']*xscale, origin='lower', **kwargs)

            common_cbar(fig, im, ax, label=term)

            ax.set_xlabel(PlotIO.x_lab)
            ax.set_ylabel(PlotIO.z_lab)

            plt.show()

        return ax
    
    def plot_yz(self, io, x, ylim=None, zlim=None, budget_terms=None, field_terms=None, 
                ax=None, xlabel=None, ylabel=None, tidx=None, xscale=None,**kwargs): 
        """
        Plots a slice in the yz-plane. 
        
        If `budget_terms` is a list, then plot_xy will recursively call itself for each term individually
        and produce figures of each separately. 
        
        Arguments
        ---------
        io (BudgetIO obj) : BudgetIO object linked to padeops data
        budget_terms (list or str) : terms to plot
        x (float) : x-location to pull slices from
        ylim, zlim : slice bounds, see BudgetIO.slice()
        
        Returns 
        -------
        
        """


        if not tidx:
            try:
                tidx = io.last_tidx
            except:
                print(" ")
            try:
                tidx = io.last_budget_tidx
            except:
                print(" ")
        if budget_terms:
            slices = io.slice(budget_terms=budget_terms, xlim=x, ylim=ylim, zlim=zlim, tidx=tidx)
            terms_list = budget_terms
        elif field_terms:
            slices = io.slice(field_terms=field_terms, xlim=x, ylim=ylim, zlim=zlim, tidx=tidx)
            terms_list = field_terms

        if not xscale:
            xscale=1

        for term in terms_list: 
            fig, ax = plt.subplots()
            im = ax.imshow(slices[term].T, extent=slices['extent']*xscale, origin='lower', **kwargs)

            common_cbar(fig, im, ax, label=term)

            ax.set_xlabel(PlotIO.y_lab)
            ax.set_ylabel(PlotIO.z_lab)
            plt.show()

        return ax

    def bar_plot(self, io, terms, coords, color_list=None, alpha_list=None, labels=None, compute_residual=True):
        """
        Plots volume-averaged terms in a bar chart. 
        
        Arguments
        ---------
        io (BudgetIO of DeficitIO obj) : BudgetIO or DeficitIo object linked to padeops data
        terms (list or str) : terms to plot
        coords (tuple) : tuple of x, y, and z limits
        color_list (list) : list of plot colors in order of terms
        alpha_list (list) : list of alpha values in order of terms

        Returns 
        -------
        fig, ax
        """

        if not isinstance(terms, list):
            if terms == "x-mom":
                comp_str = ['DuDt', 'x']
            
                keys = [key for key in io.key if io.key[key][0] == 1]
                keys_tmp = [key for key in keys if comp_str[0] in key or comp_str[1] in key]
                terms = [key for key in keys_tmp if "fluc" not in key and "mean" not in key]

                print(terms)


            elif terms == "y-mom":
                comp_str = ['DvDt', 'y']
            
                keys = [key for key in io.key if io.key[key][0] == 1]
                keys_tmp = [key for key in keys if comp_str[0] in key or comp_str[1] in key]
                terms = [key for key in keys_tmp if "fluc" not in key and "mean" not in key]

            elif terms == "z-mom":
                comp_str = ['DwDt', 'z']
            
                keys = [key for key in io.key if io.key[key][0] == 1]
                keys_tmp = [key for key in keys if comp_str[0] in key or comp_str[1] in key]
                terms = [key for key in keys_tmp if "fluc" not in key and "mean" not in key]  

            elif terms == "MKE":
                terms = [key for key in io.key if io.key[key][0] == 2 and io.key[key][1] <= 10]

            elif terms == "TKE":
                terms = [key for key in io.key if io.key[key][0] == 3 and io.key[key][1] <= 8]
            else:
                print("Please enter a valid budget type. Options include: x-mom, y-mom, z-mom, and MKE.")
                return
            
        if not color_list:
            color_list = plt.rcParams['axes.prop_cycle'].by_key()['color']     
            print(len(color_list))

        if not alpha_list:
            alpha_list = np.ones(len(terms))      

        xid, yid, zid = io.get_xids(x=coords[0], y=coords[1], z=coords[2], return_none=True, return_slice=True)
        residual = 0

        fig, ax = plt.subplots()
        
        i = 0
        for term in terms:
            ax.bar(i, np.trapz(np.trapz(np.trapz(io.budget[term][xid,yid,zid], 
                                        x=io.zLine[zid], axis=2),
                                        x=io.yLine[yid], axis=1),
                                        x=io.xLine[xid], axis=0), linewidth=0.5, color=color_list[i], alpha=alpha_list[i])
            
            residual = residual + io.budget[term]
            
            i+=1
             
        if compute_residual:   
            ax.bar(len(terms), np.trapz(np.trapz(np.trapz(residual[xid,yid,zid], 
                                            x=io.zLine[zid], axis=2),
                                            x=io.yLine[yid], axis=1),
                                            x=io.xLine[xid], axis=0), linewidth=0.5, color='tab:gray')
            if not labels:
                labels = terms.copy()
                labels.append('Res')
                
        
            ax.set_xticks(range(len(terms)+1))
        else:
            ax.set_xticks(range(len(terms)))      


        print(terms)
        print(range(len(terms)+1))
        
        ax.set_xticklabels(labels)

        return fig, ax

    def plot_budget(self, io, budget, coords=None, fig=None, ax=None, alpha=1):
        '''
        Plots the MKE budget
        
        Arguments
        ---------
        io (BudgetIO of DeficitIO obj) : BudgetIO or DeficitIo object linked to padeops data
        coords (tuple) : tuple of x, y, and z limits
        Returns
        -------
        fig : figure object 
        ax : axes object
        '''

        if budget == "x-mom":
            comp_str = ['DuDt', 'x']
        
            keys = [key for key in io.key if io.key[key][0] == 1]
            keys_tmp = [key for key in keys if comp_str[0] in key or comp_str[1] in key]
            keys = [key for key in keys_tmp if "fluc" not in key and "mean" not in key]

            colors = self.momentum_budget_colors


        elif budget == "y-mom":
            comp_str = ['DvDt', 'y']
        
            keys = [key for key in io.key if io.key[key][0] == 1]
            keys_tmp = [key for key in keys if comp_str[0] in key or comp_str[1] in key]
            keys = [key for key in keys_tmp if "fluc" not in key and "mean" not in key]

            colors = self.momentum_budget_colors

        elif budget == "z-mom":
            comp_str = ['DwDt', 'z']
        
            keys = [key for key in io.key if io.key[key][0] == 1]
            keys_tmp = [key for key in keys if comp_str[0] in key or comp_str[1] in key]
            keys = [key for key in keys_tmp if "fluc" not in key and "mean" not in key]  

            colors = self.momentum_budget_colors

        elif budget == "MKE":
            keys = [key for key in io.key if io.key[key][0] == 2 and io.key[key][1] <= 10]
            try: 
                io.key['MKE_adv_delta_base']
                keys.append('MKE_adv_delta_base')
            except:
                print("")

        elif budget == "TKE":
            keys = [key for key in io.key if io.key[key][0] == 3 and io.key[key][1] <= 8]
            colors = self.budget_colors
        else:
            print("Please enter a valid budget type. Options include: x-mom, y-mom, z-mom, and MKE.")
            return
        
        if not coords:
            xid = (slice(0, len(io.xLine)), )
            yid = (slice(0, len(io.yLine)), )
            zid = (slice(0, len(io.zLine)), ) 
        else:
            xid, yid, zid = io.get_xids(x=coords[0], y=coords[1], z=coords[2], return_none=True, return_slice=True)

        if not fig and not ax:
            fig, ax = plt.subplots()

        residual = 0
        
        for key in keys:
            # color = [color_value for color_key, color_value in colors if color_key in key]

            ax.plot(np.mean(np.mean(io.budget[key][xid,yid,zid], axis=1), axis=0), io.zLine, label=key, alpha=alpha)
            residual += io.budget[key]

        ax.plot(np.mean(np.mean(residual[xid,yid,zid], axis=1), axis=0), 
            io.zLine, label='Residual', linestyle='--', color='black',alpha=alpha)

        ax.set_ylabel('$z/D$')
            
        return fig, ax

    def plot_var_in_time(self, io, variables, coords=None, increment=1, timeDim=1000, zScale=1):
        
        print(variables)

        if not variables:
            print("Need to set variables.")
            return

        fig, ax = plt.subplots(1, len(variables))

        tidxs = io.unique_tidx()
        times = io.unique_times()

        colors = plt.cm.rainbow(np.linspace(0, 1, len(tidxs)))

        for i in range(0,len(tidxs),increment):
            io.read_fields(tidx=tidxs[i])
            
            if len(variables) == 1:
                ax.plot(np.mean(io.field[variables[0]], axis=(0,1)), io.zLine*zScale, label = np.round(times[i]*timeDim/3600,1), 
                           color=colors[i])
            else:
                ax[0].plot(np.mean(io.field[variables[0]], axis=(0,1)), io.zLine*zScale, label = np.round(times[i]*timeDim/3600,1), 
                           color=colors[i])
                for k in range(1,len(variables),increment):
                    ax[k].plot(np.mean(io.field[variables[k]], axis=(0,1)), io.zLine*zScale, color=colors[i])
            

        sm = plt.cm.ScalarMappable(cmap='rainbow', norm=plt.Normalize(vmin=times[0]*timeDim/3600, vmax=times[-1]*timeDim/3600))
        plt.colorbar(sm, label = 'Hour')

        ax[0].set_ylabel('z/D')


        for i in range(len(variables)):
            ax[i].set_xlabel(variables[i])

        return fig, ax


        # ax[0].set_xlim([300, 302])
            
        plt.show()
 
# ----------- additional helper functions ------------

def common_cbar(fig, image, ax=None, location='right', label=None, height=1., width=0.02): 
    """
    Helper function to add a colorbar
    
    colorbar: https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph 
 
    parameters
    ----------
    fig : figure object to add a colorbar to
    image : AxesImage object or similar that has a color scale
    ax : (optional) Axes object to reference for size of colorbar
    location : string, # TODO
    height : height of colorbar as a fraction of the given axis's height (or selected, if no axis given). 
        Default 1
    width : width of colorbar as a fraction of the give axis's height
        Default 0.02
    """
    if ax is None:  # need to pick an axis from the figure...
        # ax = fig.axes[0]  
        # TODO: make this better... maybe choose the largest one? Maybe make a new one each time? ugh
        ax = common_axis(fig)
    
    h = ax.get_position().height  # axis height
    cax = fig.add_axes([ax.get_position().x1+h*width, 
                        ax.get_position().y0+h/2*(1-height), 
                        h*width, 
                        height*h])  # [start x, start y, width, height]
    cbar = fig.colorbar(image, cax=cax)
    
    if label is not None: 
        cbar.set_label(label)
    
    return cbar  # this probably will never be used


def common_axis(fig, xlabel=None, ylabel=None, title=None): 
    """
    Helper function to format common axes

    format common axes: 
    https://stackoverflow.com/questions/16150819/common-xlabel-ylabel-for-matplotlib-subplots
    """
    # create a new axis spanning the whole figure
    ax_new = fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    return ax_new  # return the newly created axis


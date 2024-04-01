import numpy as np

import matplotlib.pyplot as plt
import matplotlib.legend_handler
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec


from matplotlib.colors import Normalize

from matplotlib.cm import ScalarMappable

import scipy.interpolate

import padeopsIO.budgetkey as budgetkey
import padeopsIO.deficitkey as deficitkey

import seaborn as sns

import padeopsIO as pio


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
                              "Dt": "tab:blue",
                              "adv_mean": "tab:blue",
                              "dpd": "tab:orange",
                              "SGS": "tab:green",
                              "AD": "tab:red",
                              "Cor": "tab:purple",
                              "Geo": "tab:brown",
                              "B": "tab:gray",
                              "turb": "tab:pink"}

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
                terms_not_in_key=None, ax=None, xlabel=None, ylabel=None, tidx=None, xscale=None, fig=None, cbar=True, nonDim=1, **kwargs): 
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
        elif terms_not_in_key:
            slices = io.slice(terms_not_in_key=terms_not_in_key, xlim=xlim, ylim=ylim, zlim=z, tidx=tidx)
            terms_list = terms_not_in_key

        if not xscale:
            xscale=1

        for term in terms_list:
            if ax is None:
                fig, ax = plt.subplots()
        
            im = ax.imshow(nonDim*slices[term].T, extent=slices['extent']*xscale, origin='lower',**kwargs)

            if cbar:
                common_cbar(fig, im, ax, label=term)

            ax.set_xlabel(PlotIO.x_lab)
            ax.set_ylabel(PlotIO.y_lab)

          #  plt.show()

        return ax

    
    def plot_xz(self, io, y, xlim=None, zlim=None, budget_terms=None, field_terms=None, 
                terms_not_in_key=None, ax=None, xlabel=None, ylabel=None, tidx=None, xscale=None, fig=None, cbar=True, nonDim=1, **kwargs): 
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
        elif terms_not_in_key:
            slices = io.slice(terms_not_in_key=terms_not_in_key, xlim=xlim, ylim=y, zlim=zlim, tidx=tidx)
            terms_list = terms_not_in_key

        if not xscale:
            xscale=1

        for term in terms_list: 
            if ax is None:
                fig, ax = plt.subplots()

            im = ax.imshow(nonDim*slices[term].T, extent=slices['extent']*xscale, origin='lower', **kwargs)

            if cbar:
                common_cbar(fig, im, ax, label=term)

            ax.set_xlabel(PlotIO.x_lab)
            ax.set_ylabel(PlotIO.z_lab)

#            plt.show()

        return ax
    


    def plot_yz(self, io, x, ylim=None, zlim=None, budget_terms=None, field_terms=None, 
                terms_not_in_key=None, ax=None, xlabel=None, ylabel=None, tidx=None, xscale=None, fig=None, cbar=True, nonDim=1,  **kwargs): 
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
        elif terms_not_in_key:
            slices = io.slice(terms_not_in_key=terms_not_in_key, xlim=x, ylim=ylim, zlim=zlim, tidx=tidx)
            terms_list = terms_not_in_key

        if not xscale:
            xscale=1

        for term in terms_list: 
            if ax is None:
                fig, ax = plt.subplots()

            im = ax.imshow(nonDim*slices[term].T, extent=slices['extent']*xscale, origin='lower', **kwargs)

            if cbar:
                common_cbar(fig, im, ax, label=term)

            ax.set_xlabel(PlotIO.y_lab)
            ax.set_ylabel(PlotIO.z_lab)
            #plt.show()

        return ax


    def bar_plot(self, io, terms, coords, color_list=None, alpha_list=None, 
                 labels=None, compute_residual=True, mask = None, fig=None, ax=None, nonDim=1):
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
            terms, colors = self.get_terms(io, terms)
            
        if not color_list:
            color_list = sns.color_palette("tab10", len(terms))
            
        if not alpha_list:
            alpha_list = np.ones(len(terms))      

        print(len(terms))
        print(len(color_list))

        xid, yid, zid = io.get_xids(x=coords[0], y=coords[1], z=coords[2], return_none=True, return_slice=True)
        residual = 0

        if mask is None:
            mask = np.ones([io.nx, io.ny, io.nz])


        if not fig and not ax:
            fig, ax = plt.subplots()

        # if mask is not None:

        i = 0
        for term in terms:
            # ax.bar(i, np.trapz(mask[xid,yid,zid] * io.budget[term][xid,yid,zid] * nonDim, dx=[io.dx,io.dy,io.dz]), 
                # linewidth=0.25, color=color_list[i], alpha=alpha_list[i])
            ax.bar(i, np.trapz(np.trapz(np.trapz(mask[xid,yid,zid] * io.budget[term][xid,yid,zid], 
                                        x=io.zLine[zid], axis=2),
                                        x=io.yLine[yid], axis=1),
                                        x=io.xLine[xid], axis=0)*nonDim, linewidth=0.25, color=color_list[i], alpha=alpha_list[i])
            
            residual = residual + io.budget[term]
            i+=1
             
        if compute_residual:   
            # ax.bar(len(terms), np.trapz(mask[xid,yid,zid] * residual[xid,yid,zid] * nonDim) \
            #         * io.dx * io.dy * io.dz, linewidth=0.25, color='black')
            ax.bar(len(terms), np.trapz(np.trapz(np.trapz(mask[xid,yid,zid] * residual[xid,yid,zid], 
                                            x=io.zLine[zid], axis=2),
                                            x=io.yLine[yid], axis=1),
                                            x=io.xLine[xid], axis=0)*nonDim, linewidth=0.25, color='tab:gray')
            if not labels:
                labels = terms.copy()
                labels.append('Res')
                
        
            ax.set_xticks(range(len(terms)+1))
        else:
            ax.set_xticks(range(len(terms)))      
        
        ax.set_xticklabels(labels)

        return fig, ax

    def bar_plot_streamtube(self, io, terms, xcoords, color_list=None, alpha_list=None, 
                 labels=None, compute_residual=True, mask = None, fig=None, ax=None, nonDim=1):
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
            keys, colors = self.get_terms(io, terms)
            
        if not color_list:
            color_list = plt.rcParams['axes.prop_cycle'].by_key()['color']     
            
        if not alpha_list:
            alpha_list = np.ones(len(terms))      

        x1 = np.argmin(np.abs(io.xLine - xcoords[0]))
        x2 = np.argmin(np.abs(io.xLine - xcoords[1]))

        print(x1)
        print(x2)
        residual = 0


        if not fig and not ax:
            fig, ax = plt.subplots()

        i = 0
        for term in terms:
            tmp = np.sum(mask * io.budget[term], (1,2)) * io.dy * io.dz
            ax.bar(i, np.sum(tmp[x1:x2] * nonDim) * io.dx, 
                    linewidth=0.25, color=color_list[i], alpha=alpha_list[i])
            # ax.bar(i, np.trapz(np.trapz(np.trapz(mask[xid,yid,zid] * io.budget[term][xid,yid,zid], 
            #                             x=io.zLine[zid], axis=2),
            #                             x=io.yLine[yid], axis=1),
            #                             x=io.xLine[xid], axis=0)*nonDim, linewidth=0.25, color=color_list[i], alpha=alpha_list[i])
            
            residual = residual + io.budget[term]
            i+=1
             
        if compute_residual:   
            tmp_res = np.sum(mask * residual, (1,2)) * io.dy * io.dz
            ax.bar(len(terms), np.sum(residual[x1:x2] * nonDim) * io.dx,
                    linewidth=0.25, color='tab:gray')
            # ax.bar(len(terms), np.trapz(np.trapz(np.trapz(mask[xid,yid,zid] * residual[xid,yid,zid], 
            #                                 x=io.zLine[zid], axis=2),
            #                                 x=io.yLine[yid], axis=1),
            #                                 x=io.xLine[xid], axis=0)*nonDim, linewidth=0.25, color='tab:gray')
            if not labels:
                labels = terms.copy()
                labels.append('Res')
                
        
            ax.set_xticks(range(len(terms)+1))
        else:
            ax.set_xticks(range(len(terms)))      
        
        ax.set_xticklabels(labels)

        return fig, ax
    
    def grouped_bar_plot(self, io_list, terms, coords, color_list=None, alpha_list=None, 
                         labels=None, compute_residual=True, fig=None, ax=None, legend_labels=None, nonDim=[1,1,1,1,1],
                         wake=False):
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
        # define width of bars based on number of ios
        width = 0.8/len(io_list)

        # assumes all the terms are the same for now
        # run through all ios in order to make sure RANS budgets are calculated
        if not isinstance(terms, list):
            for io in io_list:
                terms_tmp, colors = self.get_terms(io, terms)

            terms = terms_tmp

        # option to append "_wake" to terms
        if not isinstance(wake, list):
            wake = [False] * len(io_list)

        if not color_list:
            color_list = sns.color_palette("tab10", len(terms))    

        if not alpha_list:
            alpha_list = np.ones(len(terms))      

        xid_list = []
        yid_list = []
        zid_list = []
        res_list = []
        for io in io_list:
            xid, yid, zid = io.get_xids(x=coords[0], y=coords[1], z=coords[2], return_none=True, return_slice=True)
            residuals = np.zeros([xid.stop-xid.start, yid.stop-yid.start, zid.stop-zid.start])
            
            xid_list.append(xid)
            yid_list.append(yid)
            zid_list.append(zid)

            res_list.append(residuals)
        if not fig and not ax:
            fig, ax = plt.subplots()
        
        i = 0
        for term in terms:
            for j in range(len(io_list)):
                if wake[j]:
                    term = term + "_wake"
                xid = xid_list[j]
                yid = yid_list[j]
                zid = zid_list[j]
                ax.bar(i + j*width, nonDim[j]*np.trapz(np.trapz(np.trapz(io_list[j].budget[term][xid,yid,zid], 
                                            x=io_list[j].zLine[zid], axis=2),
                                            x=io_list[j].yLine[yid], axis=1),
                                            x=io_list[j].xLine[xid], axis=0), width=width, color=color_list[j], alpha=alpha_list[i])
                res_list[j] = res_list[j] + io_list[j].budget[term][xid,yid,zid]
            i+=1
     
        legend = []

        if compute_residual:   
            for j in range(len(io_list)):
                xid = xid_list[j]
                yid = yid_list[j]
                zid = zid_list[j]
                legend.append(ax.bar(len(terms) + j*width, nonDim[j]*np.trapz(np.trapz(np.trapz(res_list[j], 
                                                x=io_list[j].zLine[zid], axis=2),
                                                x=io_list[j].yLine[yid], axis=1),
                                                x=io_list[j].xLine[xid], axis=0), width=width, color=color_list[j]))
            if not labels:
                labels = terms.copy()
                labels.append('Res')
                
        
            ax.set_xticks(range(len(terms)+1))
        else:
            ax.set_xticks(range(len(terms)))      
        
        ax.set_xticklabels(labels)

        if not legend_labels:
            legend_labels = ['unstable','neutral', 'stable']
        ax.legend(handles=legend, labels=legend_labels)

        return fig, ax

    def grouped_bar_plot2(self, io_list, terms, coords, color_list=None, alpha_list=None, 
                         labels=None, compute_residual=True,  masks=None, fig=None, ax=None, 
                         legend_labels=None, nonDim=[1,1,1,1,1], wake = False):
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
        # define width of bars based on number of ios
        width = 0.8/len(io_list)

        # assumes all the terms are the same for now
        # run through all ios in order to make sure RANS budgets are calculated
        if not isinstance(terms, list):
            for io in io_list:
                terms_tmp, colors = self.get_terms(io, terms)

            terms = terms_tmp

        # option to append "_wake" to terms
        if not isinstance(wake, list):
            wake = [False] * len(io_list)

        if masks is None:
            masks = []
            for io in io_list:
                masks.append(np.ones([io.nx, io.ny, io.nz]))

        if not color_list:
            color_list = sns.color_palette("tab10", len(terms)) 

        if not alpha_list:
            alpha_list = np.ones(len(terms))      

        xid_list = []
        yid_list = []
        zid_list = []
        res_list = []
        for io in io_list:
            xid, yid, zid = io.get_xids(x=coords[0], y=coords[1], z=coords[2], return_none=True, return_slice=True)
            residuals = np.zeros([xid.stop-xid.start, yid.stop-yid.start, zid.stop-zid.start])

            xid_list.append(xid)
            yid_list.append(yid)
            zid_list.append(zid)

            res_list.append(residuals)
        if not fig and not ax:
            fig, ax = plt.subplots()
        
        i = 0
        legend = []            
        for term in terms:
            print(term)
            for j in range(len(io_list)):
                # if wake[j]:
                #     term = term + "_wake"
                # else:
                #     term = term.replace('_wake', "")
                xid = xid_list[j]
                yid = yid_list[j]
                zid = zid_list[j]
                if i==0:
                    legend.append(ax.bar(i + j*width, nonDim[j]*np.trapz(np.trapz(np.trapz(io_list[j].budget[term][xid,yid,zid] \
                                            * masks[j][xid,yid,zid], 
                                            x=io_list[j].zLine[zid], axis=2),
                                            x=io_list[j].yLine[yid], axis=1),
                                                       x=io_list[j].xLine[xid], axis=0), width=width, color=color_list[i], alpha=alpha_list[j]))
                else:
                    ax.bar(i + j*width, nonDim[j]*np.trapz(np.trapz(np.trapz(io_list[j].budget[term][xid,yid,zid] \
                                            * masks[j][xid,yid,zid], 
                                            x=io_list[j].zLine[zid], axis=2),
                                            x=io_list[j].yLine[yid], axis=1),
                                                       x=io_list[j].xLine[xid], axis=0), width=width, color=color_list[i], alpha=alpha_list[j])
                res_list[j] = res_list[j] + io_list[j].budget[term][xid,yid,zid]
            i+=1
     

        if compute_residual:   
            for j in range(len(io_list)):
                xid = xid_list[j]
                yid = yid_list[j]
                zid = zid_list[j]
                ax.bar(len(terms) + j*width, nonDim[j]*np.trapz(np.trapz(np.trapz(res_list[j] * masks[j][xid,yid,zid], 
                                                x=io_list[j].zLine[zid], axis=2),
                                                x=io_list[j].yLine[yid], axis=1),
                                                x=io_list[j].xLine[xid], axis=0), width=width, color='black', alpha=alpha_list[j])
            if not labels:
                labels = terms.copy()
                labels.append('Res')
                
        
            ax.set_xticks(range(len(terms)+1))
        else:
            ax.set_xticks(range(len(terms)))      
        
        ax.set_xticklabels(labels)

        if not legend_labels:
            legend_labels = ['unstable','neutral', 'stable']
        ax.legend(handles=legend, labels=legend_labels)

        return fig, ax

    def grouped_bar_plot2_avg(self, io_list, terms, coords, color_list=None, alpha_list=None, 
                         labels=None, compute_residual=True,  masks=None, fig=None, ax=None, 
                         legend_labels=None, nonDim=[1,1,1,1,1], wake = False):
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
        # define width of bars based on number of ios
        width = 0.8/len(io_list)

        # assumes all the terms are the same for now
        # run through all ios in order to make sure RANS budgets are calculated
        if not isinstance(terms, list):
            for io in io_list:
                terms_tmp, colors = self.get_terms(io, terms)

            terms = terms_tmp

        # option to append "_wake" to terms
        if not isinstance(wake, list):
            wake = [False] * len(io_list)

        if masks is None:
            masks = []
            for io in io_list:
                masks.append(np.ones([io.nx, io.ny, io.nz]))

        if not color_list:
            color_list = sns.color_palette("tab10", len(terms)) 

        if not alpha_list:
            alpha_list = np.ones(len(terms))      

        xid_list = []
        yid_list = []
        zid_list = []
        res_list = []
        for io in io_list:
            xid, yid, zid = io.get_xids(x=coords[0], y=coords[1], z=coords[2], return_none=True, return_slice=True)
            residuals = np.zeros([xid.stop-xid.start, yid.stop-yid.start, zid.stop-zid.start])

            xid_list.append(xid)
            yid_list.append(yid)
            zid_list.append(zid)

            res_list.append(residuals)
        if not fig and not ax:
            fig, ax = plt.subplots()
        
        i = 0
        legend = []            
        for term in terms:
            print(term)
            for j in range(len(io_list)):
                # if wake[j]:
                #     term = term + "_wake"
                # else:
                #     term = term.replace('_wake', "")
                xid = xid_list[j]
                yid = yid_list[j]
                zid = zid_list[j]
                if i==0:
                    legend.append(ax.bar(i + j*width, nonDim[j]*np.sum(io_list[j].budget[term][xid,yid,zid] \
                                            * masks[j][xid,yid,zid])/np.sum(masks[j][xid,yid,zid]), 
                                         width=width, color=color_list[i], alpha=alpha_list[j]))
                else:
                    ax.bar(i + j*width, nonDim[j]*np.sum(io_list[j].budget[term][xid,yid,zid] \
                                            * masks[j][xid,yid,zid])/np.sum(masks[j][xid,yid,zid]), width=width, color=color_list[i], alpha=alpha_list[j])
                res_list[j] = res_list[j] + io_list[j].budget[term][xid,yid,zid]
            i+=1
     

        if compute_residual:   
            for j in range(len(io_list)):
                xid = xid_list[j]
                yid = yid_list[j]
                zid = zid_list[j]
                ax.bar(len(terms) + j*width, nonDim[j]*np.sum(res_list[j] * masks[j][xid,yid,zid])/np.sum(masks[j][xid,yid,zid]), 
                                                width=width, color='black', alpha=alpha_list[j])
            if not labels:
                labels = terms.copy()
                labels.append('Res')
                
        
            ax.set_xticks(range(len(terms)+1))
        else:
            ax.set_xticks(range(len(terms)))      
        
        ax.set_xticklabels(labels)

        if not legend_labels:
            legend_labels = ['unstable','neutral', 'stable']
        ax.legend(handles=legend, labels=legend_labels)

        return fig, ax


    def grouped_bar_plot3(self, io_dict, coords, labels=None, color_list=None, mask=False, fig=None, ax=None):
        """
        Plots volume-averaged terms in a bar chart. 
        
        Arguments
        ---------
        io_dict : dictionary containing BudgetIO or Deficit IO objects

        Returns 
        -------
        fig, ax
        """
        # define width of bars based on number of ios
        width = 0.8/len(io_dict)

        if not color_list:
            color_list = sns.color_palette("tab10", len(terms))    
        for io in io_dict:
            xid, yid, zid = io_dict[io]['io'].get_xids(x=coords[0], y=coords[1], z=coords[2], return_none=True, return_slice=True)
            residuals = np.zeros([xid.stop-xid.start, yid.stop-yid.start, zid.stop-zid.start])

            io_dict[io]['xid'] = xid
            io_dict[io]['yid'] = yid
            io_dict[io]['zid'] = zid

            io_dict[io]['residuals'] = residuals

            if mask:
                io_dict[io]['mask'] = io_dict[io]['io'].stream_mask
            else:
                io_dict[io]['mask'] = np.ones([io_dict[io]['io'].nx, io_dict[io]['io'].ny, io_dict[io]['io'].nz])
        if not fig and not ax:
            fig, ax = plt.subplots()
        
        legend = [] 

        j = 0
        for io in io_dict:
            i = 0
            for term in io_dict[io]['terms']:     
                xid = io_dict[io]['xid']
                yid = io_dict[io]['yid']
                zid = io_dict[io]['zid']
                nonDim = io_dict[io]['nonDim']

                if i==0:
                    legend.append(ax.bar(i + j*width, nonDim*np.trapz(np.trapz(np.trapz(io_dict[io]['io'].budget[term][xid,yid,zid] * io_dict[io]['mask'][xid,yid,zid],
                                            x=io_dict[io]['io'].zLine[zid], axis=2),
                                            x=io_dict[io]['io'].yLine[yid], axis=1),
                                            x=io_dict[io]['io'].xLine[xid], axis=0), width=width, color=color_list[i], alpha=io_dict[io]['alpha']))
                else:
                    ax.bar(i + j*width, nonDim*np.trapz(np.trapz(np.trapz(io_dict[io]['io'].budget[term][xid,yid,zid] * io_dict[io]['mask'][xid,yid,zid], 
                                            x=io_dict[io]['io'].zLine[zid], axis=2),
                                            x=io_dict[io]['io'].yLine[yid], axis=1),
                                                       x=io_dict[io]['io'].xLine[xid], axis=0), width=width, color=color_list[i], alpha=io_dict[io]['alpha'])
                io_dict[io]['residuals'] = io_dict[io]['residuals'] + io_dict[io]['io'].budget[term][xid,yid,zid]
                i+=1
            j+=1
     
        k = 0
        for io in io_dict:
            xid = io_dict[io]['xid']
            yid = io_dict[io]['yid']
            zid = io_dict[io]['zid']
            nonDim = io_dict[io]['nonDim']


            ax.bar(len(io_dict[io]['terms']) + k*width, nonDim*np.trapz(np.trapz(np.trapz(io_dict[io]['residuals'] * io_dict[io]['mask'][xid,yid,zid], 
                                            x=io_dict[io]['io'].zLine[zid], axis=2),
                                            x=io_dict[io]['io'].yLine[yid], axis=1),
                                            x=io_dict[io]['io'].xLine[xid], axis=0), width=width, color='black', alpha=io_dict[io]['alpha'])      
            k+=1      
        print(i)       
        
        ax.set_xticks(range(i+1))
        ax.set_xticklabels(labels)

        legend_labels = [io_dict[io]['legend_label'] for io in io_dict]
        ax.legend(handles=legend, labels=legend_labels)

        return fig, ax
    
    def stacked_bar_plot(self, io_list, terms, coords, labels=None, compute_residual=True, fig=None, ax=None,
                         alpha_list=None, color_list=None):
        """
        Plots volume-averaged terms from multiple io objects in a stacked bar chart. 
        
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

        if not fig and not ax:
            fig, ax = plt.subplots()

        # first budget goes at the bottom so bottom_bool is set to false
        # io_b is set to the first io and then terms are added later in the loop
        bottom_bool = False
        j = 0
        io_b = io_list[0]
        residual_b = 0
        for io in io_list:
            if not isinstance(terms, list):
                keys, colors = self.get_terms(io, terms)

            xid, yid, zid = io.get_xids(x=coords[0], y=coords[1], z=coords[2], return_none=True, return_slice=True)
            residual = 0
            
            i = 0
            for term in terms:
                if not bottom_bool:
                    ax.bar(i, np.trapz(np.trapz(np.trapz(io.budget[term][xid,yid,zid], 
                                                x=io.zLine[zid], axis=2),
                                                x=io.yLine[yid], axis=1),
                                                x=io.xLine[xid], axis=0), linewidth=0.5, color=color_list[j])
                else:
                    ax.bar(i, np.trapz(np.trapz(np.trapz(io.budget[term][xid,yid,zid], 
                                x=io.zLine[zid], axis=2),
                                x=io.yLine[yid], axis=1),
                                x=io.xLine[xid], axis=0),
                                bottom = np.trapz(np.trapz(np.trapz(io_b.budget[term][xid,yid,zid], 
                                x=io_b.zLine[zid], axis=2),
                                x=io_b.yLine[yid], axis=1),
                                x=io_b.xLine[xid], axis=0), linewidth=0.5, color=color_list[j])
                    
                    io_b.budget[term] += io.budget[term]
                residual = residual + io.budget[term]
                residual_b = residual_b + io.budget[term]
                
                i+=1
                
            if compute_residual:
                if not bottom_bool:   
                    ax.bar(len(terms), np.trapz(np.trapz(np.trapz(residual[xid,yid,zid], 
                                                    x=io.zLine[zid], axis=2),
                                                    x=io.yLine[yid], axis=1),
                                                    x=io.xLine[xid], axis=0), linewidth=0.5, color=color_list[j])
                    
                else:
                    ax.bar(len(terms), np.trapz(np.trapz(np.trapz(residual[xid,yid,zid], 
                                                    x=io.zLine[zid], axis=2),
                                                    x=io.yLine[yid], axis=1),
                                                    x=io.xLine[xid], axis=0),
                                        bottom=np.trapz(np.trapz(np.trapz(residual_b[xid,yid,zid], 
                                                    x=io_b.zLine[zid], axis=2),
                                                    x=io_b.yLine[yid], axis=1),
                                                    x=io_b.xLine[xid], axis=0),linewidth=0.5, color=color_list[j])
                residual_b = residual
                if not labels:
                    labels = terms.copy()
                    labels.append('Res')
                    
            
                ax.set_xticks(range(len(terms)+1))
            else:
                ax.set_xticks(range(len(terms)))      
            
            ax.set_xticklabels(labels)

            j+=1
            bottom_bool = True

        return fig, ax

    def plot_integrated_budget(self, io, budget, dims=None, coords=None, fig=None, ax=None):
        '''
        Plot integrated budget
        '''

        if not isinstance(budget, list):
            keys,_ = sefl.get_terms(io, budget)

        if not coords:
            xid = slice(0, len(io.xLine))
            yid = slice(0, len(io.yLine))
            zid = slice(0, len(io.zLine)) 
        else:
            xid, yid, zid = io.get_xids(x=coords[0], y=coords[1], z=coords[2], return_none=True, return_slice=True)

        if not dims:
            dims = ['x','y','z']


    def plot_budget(self, io, budget, coords=None, fig=None, ax=None, alpha=1, zScale=1, nonDim=1, color_list=plt.rcParams['axes.prop_cycle'].by_key()['color'], **kwargs):
        '''
        Plot x and y averaged budgets
        
        Arguments
        ---------
        io (BudgetIO of DeficitIO obj) : BudgetIO or DeficitIo object linked to padeops data
        coords (tuple) : tuple of x, y, and z limits
        Returns
        -------
        fig : figure object 
        ax : axes object
        '''

        if not isinstance(budget,list):
            keys, _ = self.get_terms(io, budget)
        
        if not coords:
            xid = slice(0, len(io.xLine))
            yid = slice(0, len(io.yLine))
            zid = slice(0, len(io.zLine)) 
        else:
            xid, yid, zid = io.get_xids(x=coords[0], y=coords[1], z=coords[2], return_none=True, return_slice=True)


        if not fig and not ax:
            fig, ax = plt.subplots()

        residual = 0
        
        i=0
        for key in keys:
            # color = [color_value for color_key, color_value in colors if color_key in key]
            color_list = plt.rcParams['axes.prop_cycle'].by_key()['color']     

            ax.plot(np.mean(io.budget[key][xid,yid,zid], axis=(0,1))*nonDim, io.zLine[zid]*zScale, label=key, alpha=alpha,  color=color_list[i], **kwargs)

            residual += io.budget[key]
            i+=1

        ax.plot(np.mean(np.mean(residual[xid,yid,zid], axis=1), axis=0)*nonDim, 
            io.zLine[zid]*zScale, label='Residual', linestyle='--', color='black',alpha=alpha)

        ax.set_ylabel('$z/D$')
            
        return fig, ax

    def plot_budget_in_x(self, io, budget, coords=None, fig=None, ax=None, alpha=1, color_list=plt.rcParams['axes.prop_cycle'].by_key()['color'], 
                         xScale=1, nonDim=1, linestyle='-', mask=None, smooth=False, window=12, labels=None, res_bool=True, **kwargs):
        '''
        Plots a given budget
        
        Arguments
        ---------
        io (BudgetIO of DeficitIO obj) : BudgetIO or DeficitIo object linked to padeops data
        coords (tuple) : tuple of x, y, and z limits
        Returns
        -------
        fig : figure object 
        ax : axes object
        '''

        if not isinstance(budget, list):
            keys, _ = self.get_terms(io, budget)
        else:
            keys = budget

        if mask is None:
            mask = np.ones([io.nx, io.ny, io.nz])
            
        if labels is None:
            labels = keys

        if not coords:
            xid = slice(0, len(io.xLine))
            yid = slice(0, len(io.yLine))
            zid = slice(0, len(io.zLine))
        else:
            xid, yid, zid = io.get_xids(x=coords[0], y=coords[1], z=coords[2], return_none=True, return_slice=True)


        if not fig and not ax:
            fig, ax = plt.subplots()

        residual = 0

        print(keys)
        
        i=0
        for key in keys:
            budget_term = (1/np.count_nonzero(mask[:,yid,zid])) * np.sum(np.multiply(io.budget[key][:,yid,zid], mask[:,yid,zid]), axis=(1,2))*nonDim
            if smooth:
                budget_term = smooth_data(budget_term, window)
            ax.plot(io.xLine[xid]*xScale, budget_term[xid],
                    label=labels[i], alpha=alpha,  color=color_list[i], linestyle=linestyle, **kwargs)

            residual += io.budget[key]
            i+=1


        if res_bool:
            res = (1/np.count_nonzero(mask[:,yid,zid])) * np.sum(np.multiply(residual[:,yid,zid], mask[:,yid,zid]), axis=(1,2))*nonDim
            if smooth:
                res = smooth_data(res, window)
            ax.plot(io.xLine[xid]*xScale, res[xid],
                    label=r'$\rm{Res}$', linestyle='--', color='black',alpha=alpha, **kwargs)

        ax.set_xlabel('$x/D$')
            
        return fig, ax


    def plot_var_in_time(self, io, variables, coords=None, increment=1, timeDim=1000, zScale=1):

        if not variables:
            print("Need to set variables.")
            return

        fig, ax = plt.subplots(1, len(variables))

        if not coords:
            xid = slice(0, len(io.xLine))
            yid = slice(0, len(io.yLine))
            zid = slice(0, len(io.zLine))
        else:
            xid, yid, zid = io.get_xids(x=coords[0], y=coords[1], z=coords[2], return_none=True, return_slice=True)


        tidxs = io.unique_tidx()
        times = io.unique_times()

        colors = plt.cm.rainbow(np.linspace(0, 1, len(tidxs)))

        for i in range(0,len(tidxs),increment):
            io.read_fields(tidx=tidxs[i])
            
            if len(variables) == 1:
                ax.plot(np.mean(io.field[variables[0]][xid,yid,zid], axis=(0,1)), io.zLine[zid]*zScale, label = np.round(times[i]*timeDim/3600,1), 
                           color=colors[i])
            else:
                ax[0].plot(np.mean(io.field[variables[0]][xid,yid,zid], axis=(0,1)), io.zLine[zid]*zScale, label = np.round(times[i]*timeDim/3600,1), 
                           color=colors[i])
                for k in range(1,len(variables),increment):
                    ax[k].plot(np.mean(io.field[variables[k]][xid,yid,zid], axis=(0,1)), io.zLine[zid]*zScale, color=colors[i])
            

        sm = plt.cm.ScalarMappable(cmap='rainbow', norm=plt.Normalize(vmin=times[0]*timeDim/3600, vmax=times[-1]*timeDim/3600))
        plt.colorbar(sm, label = 'Hour')

        ax[0].set_ylabel('z/D')


        for i in range(len(variables)):
            ax[i].set_xlabel(variables[i])

        return fig, ax


    # Plot slices of the data at the given coordinates
    def plot_slices(self, io, xloc, key, ylim=None, zlim=None, ax=None, cbar=False, nonDim=1, vmin=-1, vmax=1):

        if ylim is None:
            ylim = [io.yLine[0], io.yLine[-1]]
        
        if zlim is None:
            zlim = [io.zLine[0], io.zLine[-1]]

        xid, yid, zid = io.get_xids(x=xloc, y=ylim, z=zlim, return_none=True, return_slice=True)

        x = io.xLine
        y = io.yLine[yid]
        z = io.zLine[zid]
        data = io.budget[key][...,yid,zid]

        if ax is None:
            ax = plt.figure().add_subplot(111, projection='3d')
        # Take slices interpolating to allow for arbitrary values
        data_x = scipy.interpolate.interp1d(x, data, axis=0, kind='quadratic')(xloc)

        # Plot X slice
        
        num_levels = 22
        midpoint = 0
        levels = np.linspace(vmin, vmax, num_levels)
        midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
        vals = np.interp(midp, [vmin, midpoint, vmax], [0, 0.5, 1])
        colorchoice = plt.cm.bwr(vals)
        cmap, norm = colors.from_levels_and_colors(levels, colorchoice)  
        xs, ys, zs = data.shape
        xplot = ax.plot_surface(xloc, y[:, np.newaxis], z[np.newaxis, :], rstride=1, cstride=1, facecolors=cmap(norm(data_x*nonDim)), shade=False)

        if cbar:
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            plt.colorbar(sm, ax=ax, shrink=0.1, label='$\overline{\Delta u}\; (m/s)$')
            
        ax.invert_yaxis()
        
        return xplot
    

    def get_terms(self, io, budget):
        """
        Returns the terms for the selected budget
        """

        if 'x' in budget:
            comp_str = ['DuDt', 'x']

            colors = self.momentum_budget_colors


        elif 'y' in budget:
            comp_str = ['DvDt', 'y']

            colors = self.momentum_budget_colors

        elif 'z' in budget:
            comp_str = ['DwDt', 'z']

            colors = self.momentum_budget_colors

        elif budget == "MKE":
            keys = [key for key in io.key if io.key[key][0] == 2 and io.key[key][1] <= 10]
            if isinstance(io, pio.DeficitIO):
                keys.append('MKE_adv_delta_base')

            colors = self.budget_colors

        elif budget == "TKE":
            keys = [key for key in io.key if io.key[key][0] == 3 and io.key[key][1] <= 11]
            if 'TKE_adv_delta_base' in io.key.keys():
                keys.insert(2,'TKE_adv_delta_base')
            if 'TKE_turb_transport_delta_base' in io.key.keys():
                keys.insert(2,'TKE_turb_transport_delta_base')
            if 'TKE_prod_delta_delta' in io.key.keys():
                # keys.remove('TKE_production')
                keys.remove('TKE_prod_delta_delta')
                keys.remove('TKE_prod_base_delta')
            colors = self.budget_colors
        else:
            print("Please enter a valid budget type. Options include: x-mom, y-mom, z-mom, and MKE.")
            return

        if 'mom' in budget:
            keys = [key for key in io.key if io.key[key][0] == 1]
            keys_tmp = [key for key in keys if comp_str[0] in key or comp_str[1] in key]
            keys = [key for key in keys_tmp if "fluc" not in key and "mean" not in key]  

        if 'model' in budget:
            keys = [key for key in io.key if io.key[key][0] == 1]
            keys_tmp = [key for key in keys if comp_str[0] in key or comp_str[1] in key]
            keys = [key for key in keys_tmp if "fluc" not in key and "mean" not in key]  

            keys.insert(0, comp_str[1] + 'adv_mean_x')
            keys.insert(1, comp_str[1] + 'adv_mean_x')
            keys.insert(2, comp_str[1] + 'adv_mean_yz_base')
            keys.append(comp_str[1] + 'adv_mean_yz_delta_total')
            
        if 'RANS' in budget:
            # gather all the keys in budget 1 for x
            keys = [key for key in io.key if io.key[key][0] == 1]
            keys = [key for key in keys if comp_str[1] in key]


            # check for RANS terms for BudgetIO object
            # if not there, run rans_calc()
            if comp_str[1] + 'Adv_total' in io.key.keys():

                keys.remove(comp_str[1] + 'Adv_total')

                if comp_str[1] + 'turb' not in io.key.keys():
                    io.budget[comp_str[1] + 'turb'] = io.budget[comp_str[1] + 'Adv_base_delta_fluc'] \
                                                + io.budget[comp_str[1] + 'Adv_delta_delta_fluc'] \
                                                + io.budget[comp_str[1] + 'Adv_delta_base_fluc']

                if comp_str[1] + 'adv_mean' not in io.key.keys():
                    io.budget[comp_str[1]+ 'adv_mean'] = io.budget[comp_str[1] + 'Adv_base_delta_mean'] \
                                + io.budget[comp_str[1] + 'Adv_delta_delta_mean']

                keys.remove(comp_str[1] + 'Adv_base_delta_fluc')
                keys.remove(comp_str[1] + 'Adv_delta_delta_fluc')
                keys.remove(comp_str[1] + 'Adv_delta_base_fluc')
                keys.remove(comp_str[1] + 'Adv_base_delta_mean')
                keys.remove(comp_str[1] + 'Adv_delta_delta_mean')

            elif comp_str[1] + 'turb' not in io.key.keys() or comp_str[1] + 'adv_mean' not in io.key.keys():
                io.rans_calc()

            keys.insert(0,comp_str[1] + 'adv_mean')
            keys.append(comp_str[1] + 'turb')


        return keys, colors

    def yz_slices(self, io_list, terms, labels, xloc, nonDim_list,norm, cbar_label=r"$\mathcal{P}_{ijk}\;D/U^3$", residual=False,cmap='PuOr_r',  **kwargs):
    
        fig = plt.figure(figsize=(16,5))

        if residual:
            gs1 = gridspec.GridSpec(len(io_list), len(terms)+1)
            for io in io_list:
                io.budget['res'] = np.zeros([io.nx, io.ny, io.nz])
        else:
            gs1 = gridspec.GridSpec(len(io_list), len(terms))
        gs1.update(wspace=0.1, hspace=0.0) # set the spacing between axes. 
    

        axes = []
        i=0
        for term in terms:
            for j in range(len(io_list)):
                ax1 = fig.add_subplot(gs1[i+j*len(terms)])
                self.plot_yz(io_list[j], xloc, ylim=[-2.5,2.5], zlim=[-90/126, 3], terms_not_in_key=[term], fig=fig, ax=ax1, nonDim=nonDim_list[j], 
                                cbar = False, interpolation='bilinear', cmap=cmap, norm=norm, **kwargs)
                
                if residual:
                    io_list[j].budget['res'] += io_list[j].budget[term]

                if i!=0:
                    ax1.set_ylabel("")
                    ax1.set_yticklabels([])
                else:
                    ax1.set_ylabel(r'$z/D$')

                if j!=0:
                    ax1.set_xlabel(r'$y/D$')
                else:
                    ax1.set_xlabel("")
                    ax1.set_xticklabels([])
                    ax1.set_title(labels[i])

                axes.append(ax1)
        
            i+=1

        if residual:
            for i in range(len(io_list)):
                ax1 = fig.add_subplot(gs1[(len(terms)+1)*(i+1)])

                pioPlot.plot_yz(io_list[i], xloc, ylim=[-2.5,2.5], zlim=[-90/126, 3], terms_not_in_key=['res'], fig=fig, ax=ax1, nonDim=nonDim_list[i], 
                                cbar = False, interpolation='bilinear', cmap=cmap, norm=norm,  **kwargs)


        sm = plt.cm.ScalarMappable(cmap=cmap)
    
        sm.set_norm(norm)

        cb = fig.colorbar(sm, ax = axes, label=cbar_label, shrink=0.9, 
                          pad=0.025, spacing='proportional')
        cb.ax.set_yscale('linear')

        return fig, axes

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

def make_handle_tuples(ax):
    """
    Helper function to make double column figure legend
    Input
    ax - axis object
    Output
    handles_tuples - tuples of handles to be used by the matplotlib.legend_handler
    labels - legend labels
    Usage:
    from matplotlib.legend_handler HandleTuple

    <...plot figure...>

    handle_tuples, labels = pio.padeplots.make_handle_tuples(ax)

    fig.legend(handle_tuples, labels, handler_map={tuple: HandlerTuple(ndivide=None)})
    """
    handles, labels = ax.get_legend_handles_labels()
    
    handles_tuples = []
    length = int(len(handles)/2)
    
    for i in range(length):
        tmp_tuple = (handles[i], handles[i + length])
        handles_tuples.append(tmp_tuple)
        
    return handles_tuples, labels[:length]
    

def smooth_data(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode="same")
    return y_smooth

def draw_line(ax, origin, width, height):
    # Create a Rectangle patch
    rect = patches.Rectangle(origin, width=width, height=height, linewidth=0.5, edgecolor='gray', facecolor='white')

    # Add the patch to the Axes
    ax.add_patch(rect)

def make_colorbar(clim, cmap):
    
    norm = Normalize(vmin=clim[0], vmax=clim[1])

    cm = ScalarMappable(norm, cmap=cmap)

    return cm

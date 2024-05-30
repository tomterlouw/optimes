import pandas as pd 
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D
from scipy.interpolate import pchip_interpolate
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.colors import SymLogNorm
from matplotlib.colors import LinearSegmentedColormap

from config import (NAME_REF_DB, ASSESSMENT_YEAR, ASSESS_COUNTRY, COST_DATA)
from mapping import (dict_color_n, replacement_dict, name_dict_fig, color_dict_cost, 
                    column_dict_categories, hex_dict_cost, color_dict_env, hex_dict_env, dict_units, cap_mapping)

def make_stack_plot(file, list_prods, list_cons, nr_subplot, x_sub, ax_spec = 0, plot_season = 'winter', new_fig=True, storage_techs=False):
    """
    Creates a stacked plot.
    
    Args:
        file (str): File path.
        list_prods (list): List of energy production units.
        list_cons (list): List of energy consumption units.
        nr_subplot (int): Number of subplots.
        x_sub (int): X-subplot value.
        ax_spec (int): Axis specification.
        plot_season (str): Plot season.
        new_fig (bool): Whether to create a new figure.
        storage_techs: Storage technologies in the stack plot.
        fig: Existing figure.
        axs: Existing axes.

    Returns:
        ax: Axes object.
    """    

    rc = {
        'axes.facecolor':'white',
        'axes.grid' : True,
        'grid.color': '.9',
        'font.family':'Arial',
        'xtick.bottom': False,
        'xtick.top': False,
        'ytick.left': False,
        'ytick.right': False,
            'patch.facecolor': 'b',
            'patch.linewidth': 1,
        'patch.force_edgecolor': False
        }
    plt.rcParams.update(rc) 
    plt.rcParams['axes.axisbelow'] = True
    plt.rcParams.update({'mathtext.default':  'regular' })
    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 1
    
    global fig, axs
    
    general_font = 18
    
    generators_start = pd.read_excel(file)
    generators_start.set_index('time', inplace=True)

    if plot_season == "summer":
        start_date_plot = "{}-07-22 00:00:00".format(ASSESSMENT_YEAR)
        end_date_plot = "{}-07-29 00:00:00".format(ASSESSMENT_YEAR)
        generators_start = generators_start.loc[start_date_plot:end_date_plot]
        
    elif plot_season == "winter":
        start_date_plot = "{}-01-17 00:00:00".format(ASSESSMENT_YEAR)
        end_date_plot = "{}-01-24 00:00:00".format(ASSESSMENT_YEAR)
        generators_start = generators_start.loc[start_date_plot:end_date_plot] 
        
    generators_init_prod = generators_start.filter(items=list_prods)
    
    if "p_h2_ves" in list_prods:
        generators_init_prod.p_h2_ves[generators_init_prod.p_h2_ves>0] = 0
        generators_init_prod[generators_init_prod<0] = generators_init_prod[generators_init_prod<0].apply(lambda x: x*-1)
        generators_init_prod.rename(columns={'p_h2_ves': 'p_h2_ves_dis'}, inplace=True)
    
    generators_init_prod = generators_init_prod.loc[:, (generators_init_prod != 0).any(axis=0)]
    
    if storage_techs:
        storage_start = generators_start.filter(items=storage_techs)  
    
    generators_init_cons = generators_start.filter(items=list_cons).apply(lambda x: x*-1)
    generators_init_cons = generators_init_cons.loc[:, (generators_init_cons != 0).any(axis=0)]
    generators_init_cons[generators_init_cons>0] = 0
    
    list_cols_prod = []
    list_cols_cons = []

    # Get data from df for the stackplot
    for i, col in enumerate(list(generators_init_prod.columns)):
        list_cols_prod.append(generators_init_prod[generators_init_prod.columns[i]])

    for i, col in enumerate(list(generators_init_cons.columns)):
        list_cols_cons.append(generators_init_cons[generators_init_cons.columns[i]])
    
    generators_init_sum_analyze = pd.concat([generators_init_prod,generators_init_cons],axis=1).sum(axis=1)
    
    if len(generators_init_sum_analyze[generators_init_sum_analyze>1])>0:
        print("CHECK PLEASE")
        print(generators_init_sum_analyze[generators_init_sum_analyze>1])
        
    # Make a figure
    if new_fig:
        fig, axs = plt.subplots(nr_subplot, x_sub, figsize=(20, 26*nr_subplot))

    axe = axs.ravel()
    
    ax = axe[ax_spec]
    lw=2

    colums_s = list(generators_init_prod.columns) + list(generators_init_cons.columns)

    if len(list_cols_prod)>0 and len(list_cols_cons)>0:    
        ax.stackplot(list(generators_init_prod.index), list_cols_prod, colors=[dict_color_n[r] for r in list(generators_init_prod.columns)])
        ax.stackplot(list(generators_init_cons.index), list_cols_cons,  colors=[dict_color_n[r] for r in list(generators_init_cons.columns)])

        # Apply replacements using list comprehension
        list_names = [replacement_dict.get(w, w) for w in list(colums_s)]

        if len(list_names)>3:
            ax.legend(list_names, ncol=5, bbox_to_anchor=(1.02,-0.11), fontsize=general_font-3, framealpha=0.2, frameon=False, borderpad=0.2)
        else:
            ax.legend(list_names, ncol=5, bbox_to_anchor=(0.8,-0.11), fontsize=general_font-3, framealpha=0.2, frameon=False, borderpad=0.2)

        if (plot_season=='summer') or (plot_season=='winter'):
            # Only add lines to plots which don't have to many timeslots, otherwise it is hardly readable
            ax.plot(generators_init_prod.sum(axis=1), color='black', linewidth=1, linestyle='--', label="")
            ax.plot(generators_init_cons.sum(axis=1), color='black', linewidth=1, linestyle='--', label="")
    else:
        print("AXIS NOT THERE AS THERE ARE NO PROD OR CONS COLS")
        print(list_cols_prod,list_cols_cons)
        ax.axis('off')
        
    if (storage_techs) and (storage_start.sum().item() >0) and (plot_season != ""):
        ax2 = ax.twinx()
        ax2.plot(storage_start, color='black', linewidth=lw)

        if "E_batt" == storage_techs[0]:
            ax2.set_ylabel("Battery electricity\nstorage [MWh]", fontsize=general_font-2, color='black')            
        elif "E_heat_stor" == storage_techs[0]:
            ax2.set_ylabel("Residential\nheat storage [MWh]", fontsize=general_font-2, color='black')             
        else:
            ax2.set_ylabel("Hydrogen\nstorage [MWh]", fontsize=general_font-2, color='black')
        ax2.grid(False)    
        ax2.tick_params(axis='both', which='major', labelsize=general_font, rotation = 360)
        ax2.tick_params(axis='both', which='minor', labelsize=general_font, rotation = 47)
    
    ax.set_ylim(min(generators_init_cons.sum(axis=1))*1.1, max(generators_init_prod.sum(axis=1))*1.1)
    ax.grid(which="major", color="grey", linestyle="--", linewidth=0.5)

    ax.set_xlim(generators_init_cons.index[0], generators_init_cons.index[-1])
    ax.set_ylabel("Power [MW]", fontsize=general_font)
        
    ax.tick_params(axis='both', which='major', labelsize=general_font, rotation = 360)
    ax.tick_params(axis='both', which='minor', labelsize=general_font, rotation = 47)
    
    if plot_season == "":
        ax.set_xticklabels(['Jan', "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    
    fig.subplots_adjust(top=0.27, bottom=0.15, wspace = 0.1, hspace=0.72)

    return ax

def plot_weather_profiles(plot_info, weather_data, mean_prov = True):
    """
    Plots data and generation profiles

    Args:
        plot_info (dict): dictionary with data, first one is the data from dataframe, second the name above the figure, and the y_name.
        weather_data (dataframe): dataframe with annual weather data.
    """
    font = 10
    fig, axs = plt.subplots(1, len(plot_info), figsize=(16, 3), sharey=False, sharex=False)
    axe = axs.ravel()

    for i, (data_key, info) in enumerate(zip(plot_info.keys(), plot_info.values())):
        ax = axe[i]

        data = weather_data[data_key]
        data.plot(ax=ax, linewidth=1, color='darkblue')

        ax.set_ylabel(info['y_name'], fontsize=font)
        ax.set_xlabel("")
        ax.set_ylim(data.min() * 0.95, data.max() * 1.05)
        ax.set_title(info['name'], fontsize=font)

        ax.tick_params(axis='both', which='major', labelsize=font, rotation=360)
        ax.tick_params(axis='both', which='minor', labelsize=font, rotation=47)

        if mean_prov:
            ax.text("{}-09-29 23:00:00".format(ASSESSMENT_YEAR), data.max() * 0.95, s="Mean: {}".format(round(data.mean(), 1)), fontsize=font - 1)
            
        ax.set_xticklabels(['Jan', "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])

    plt.savefig("figs/weather_data.png", dpi=100, bbox_inches='tight')

def plot_g_d_profiles(data_mapping, data_frame, mean_prov = True):
    """
    Plots data and generation profiles

    Args:
        data_mapping (dict): dictionary with data, first one is the data from dataframe, second the name above the figure, and the y_name.
        data_frame (dataframe): dataframe with annual data.

    Returns:
        Plotted figure stored on local disk.
    """
    rc = {'figure.figsize':(10,4),
        'axes.facecolor':'white',
        'axes.grid' : True,
        'grid.color': '.9',
        'font.family':'Arial',
        'font.size' : 20,
        'xtick.bottom': False,
        'xtick.top': False,
        'ytick.left': False,
        'ytick.right': False,
            'patch.facecolor': 'b',
            'patch.linewidth': 3,
        'patch.force_edgecolor': False
        }

    plt.rcParams.update(rc) 
    plt.rcParams['axes.axisbelow'] = True
    plt.rcParams.update({'mathtext.default':  'regular' })
    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 1

    font = 12
    fig, axs = plt.subplots(4, 2, figsize=(16, 16), sharey=False, sharex=False)

    axe = axs.ravel()

    for i, (data_key, data_info) in enumerate(zip(data_mapping.keys(), data_mapping.values())):
        ax = axe[i]
        data = data_frame[data_key] 
        data.plot(ax=ax, linewidth=1, color='darkblue')

        ax.set_ylabel(data_info['y_name'], fontsize=font)
        ax.set_xlabel("")

        if data_key == 'electricity_demand' or data_key == 'bld_heat_demand':
            ax.set_ylim(0, data.max() * 1.02)
        else:
            ax.set_ylim(data.min() * 1.02, data.max() * 1.02)

        ax.set_title(data_info['name'], fontsize=font)

        ax.tick_params(axis='both', which='major', labelsize=font, rotation=360)
        ax.tick_params(axis='both', which='minor', labelsize=font, rotation=47)

        if mean_prov:
            ax.text("{}-09-29 23:00:00".format(ASSESSMENT_YEAR), data.max() * 0.95, s="Mean: {}".format(round(data.mean(), 3)), fontsize=font - 1)

        ax.set_xticklabels(['Jan', "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])

    plt.tight_layout()
        
    #plt.savefig("figs/profiles.eps", dpi = 100, bbox_inches='tight')    
    plt.savefig("figs/profiles.png", dpi = 150, bbox_inches='tight')  

# credits to https://matplotlib.org/3.1.3/gallery/specialty_plots/radar_chart.html
# for spide rgraphs
def radar_factory(num_vars, frame='circle'):
    """Create a radar chart with `num_vars` axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle' | 'polygon'}
        Shape of frame surrounding axes.

    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)

    class RadarAxes(PolarAxes):

        name = 'radar'

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # rotate plot such that the first axis is at the top
            self.set_theta_zero_location('N')

        def fill(self, *args, closed=True, **kwargs):
            """Override fill so that line is closed by default"""
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.concatenate((x, [x[0]]))
                y = np.concatenate((y, [y[0]]))
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
            # in axes coordinates.
            if frame == 'circle':
                return Circle((0.5, 0.5), 0.5, edgecolor="black", lw = 2)
            elif frame == 'polygon':
                return RegularPolygon((0.5, 0.5), num_vars,
                                      radius=.6, edgecolor="black", lw = 2)
            else:
                raise ValueError("unknown value for 'frame': %s" % frame)

        def draw(self, renderer):
            """ Draw. If frame is polygon, make gridlines polygon-shaped """
            if frame == 'polygon':
                gridlines = self.yaxis.get_gridlines()
                for gl in gridlines:
                    gl.get_path()._interpolation_steps = num_vars
            super().draw(renderer)


        def _gen_axes_spines(self):
            if frame == 'circle':
                return super()._gen_axes_spines()
            elif frame == 'polygon':
                # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                spine = Spine(axes=self,
                              spine_type='circle',
                              path=Path.unit_regular_polygon(num_vars))
                # unit_regular_polygon gives a polygon of radius 1 centered at
                # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                # 0.5) in axes coordinates.
                spine.set_transform(Affine2D().scale(.5).translate(.5, .5)
                                    + self.transAxes)


                return {'polar': spine}
            else:
                raise ValueError("unknown value for 'frame': %s" % frame)

    register_projection(RadarAxes)
    return theta

def generate_cost_ghg_figure(df_select_cost, df_init_gwp, df_init_cost, nr_scenarios=9, nr_left_subplot=3):
    """
    Generate a dual-subplot bar plot for visualizing the costs and life cycle GHG-emissions in different scenarios.

    Parameters:
    - df_select_cost (DataFrame): DataFrame containing selected cost data.
    - df_init_gwp (DataFrame): DataFrame containing life cycle GHG-emission data.
    - df_init_cost (DataFrame): DataFrame containing initial cost data.
    - nr_scenarios (int, optional): Number of scenarios to display in the plot. Default is 8.
    - nr_left_subplot (int, optional): Number of scenarios to display on the left subplot. Default is 3.

    Output:
    - The generated plot is saved as a PNG file in the "figs" directory with the filename "contri_costs_ghgs.png".
    - The plot is adjusted for layout, and the resulting image is saved.

    Example Usage:
    generate_cost_ghg_figure(df_select_cost, df_init_gwp, df_init_cost, nr_scenarios=8, nr_left_subplot=3)
    """
    # Total figure
    rc = {
        'axes.facecolor':'white',
        'axes.grid' : False,
        'grid.color': '.8',
        'font.family':'Arial',
        'font.size' : 15,
        'xtick.bottom': True,
        'xtick.top': True,
        'ytick.left': True,
        'ytick.right': True,
            'patch.facecolor': 'b',
            'patch.linewidth': 1,
        'patch.force_edgecolor': True,
        'axes.linewidth':2
        }

    plt.rcParams.update(rc) 
    plt.rcParams['axes.axisbelow'] = True
    plt.rcParams.update({'mathtext.default':  'regular' })
    plt.rcParams["axes.edgecolor"] = "black"   

    fig, axes = plt.subplots(1,2, figsize = (60,22), gridspec_kw={'width_ratios': [nr_left_subplot/nr_scenarios, (nr_scenarios-nr_left_subplot)/nr_scenarios]})

    # Rename the color names in the figure
    list_names = [name_dict_fig[w] for w in list(df_select_cost.T.columns)]
    fontsi = 50

    substr = "(BI)"
    for k, ax in enumerate(axes):
        df_select_cost_2, df_init_gwp_2, df_init_cost_2 = [
            df.filter(regex=substr) if k == 0 else df.drop([col for col in df.columns if substr in col], axis=1)
            for df in [df_select_cost, df_init_gwp, df_init_cost]
        ]
        df_select_cost_2.T.plot(kind='bar', ax = ax, position =1.1, width =0.27, stacked = True, zorder=30, legend=True if k==1 else False,
                            color=[color_dict_cost[r] for r in list(df_select_cost_2.T.columns)])
        
        labels = df_select_cost_2.columns
        x = np.arange(len(labels))

        bars = ax.patches
        patterns = [hex_dict_cost[z] for z in list(df_select_cost_2.T.columns)]  # set hatch patterns in the correct order
        hatches = []  # list for hatches in the order of the bars

        for h in patterns:  # loop over patterns to create bar-ordered hatches
            for i in range(int(len(bars) / len(patterns))):
                hatches.append(h)

        for bar, hatch in zip(bars, hatches):  # loop over bars and hatches to set hatches in correct order
            bar.set_hatch(hatch)

        ax2=ax.twinx()
        df_init_gwp_2.T.plot(kind='bar', ax = ax2, width =0.27, stacked = True, zorder=500, legend=False, position=0,
                            color=[color_dict_env[r] for r in list(df_init_gwp_2.T.columns)], alpha=1)
        
        bars = ax2.patches
        patterns = [hex_dict_env[z] for z in list(df_init_gwp_2.T.columns)]  # set hatch patterns in the correct order
        hatches = []  # list for hatches in the order of the bars

        for h in patterns:  # loop over patterns to create bar-ordered hatches
            for i in range(int(len(bars) / len(patterns))):
                hatches.append(h)

        for bar, hatch in zip(bars, hatches):  # loop over bars and hatches to set hatches in correct order
            bar.set_hatch(hatch)

        if k==1:
            ax.legend(labels = list_names, fontsize=fontsi-2, ncol=4, bbox_to_anchor=(1.075,-0.21), frameon=False)

        ax.set_xlabel("")
        ax.set_ylabel("Costs [M€/year]", size= fontsi+5)
        ax2.set_ylabel("Life cycle GHG-emissions [kt CO$_{2}$-eq./year]", size= fontsi+5)
        ax.xaxis.grid(False, which = 'both')
        ax.yaxis.grid(True, which = 'both')
        ax.set_xlim(-0.41,len(list(df_select_cost_2.columns))-0.52)
        
        # Calculate the maximum and minimum values for y-axis limits
        ymax = max(df_select_cost_2.applymap(lambda x: max(x, 0)).abs().sum().max(),
                df_init_gwp_2.applymap(lambda x: max(x, 0)).abs().sum().max())
        ymin = min(df_select_cost_2.applymap(lambda x: min(x, 0)).sum().min(),
                df_init_gwp_2.applymap(lambda x: min(x, 0)).sum().min())

        # Set y-axis limits with a margin
        margin = 0.1  # Adjust the margin as needed
        ax.set_ylim(ymin - margin * abs(ymin), ymax + margin * abs(ymax))
        ax2.set_ylim(ymin - margin * abs(ymin), ymax + margin * abs(ymax))

        ax.hlines(y=0, xmin=-1, xmax=30, color='black', lw=2, zorder=20)

        ax.tick_params(axis='both', which='major', labelsize=fontsi+4, rotation = 360)
        ax.tick_params(axis='both', which='minor', labelsize=fontsi, rotation = 47)

        ax2.tick_params(axis='both', which='major', labelsize=fontsi+4, rotation = 360)
        ax2.tick_params(axis='both', which='minor', labelsize=fontsi, rotation = 47)
        
        ax.set_title("a. Industry (bakery)" if k==0 else "b. Entire multi-energy system", fontsize=fontsi+8, weight='bold')

        list_total_ghgs = list(df_init_gwp_2.sum(axis=0))
        list_total_costs = list(df_init_cost_2.T.total_costs/1e3)

        for x in range(0,int(len(labels))):
            ax.text(x-0.26, 1.20*ymin if k==0 else 1.23*ymin, "Costs", rotation=90, fontsize=fontsi-4, ha='center', va='center')
            ax.text(x+0.1, 1.20*ymin if k==0 else 1.23*ymin, "GHGs", rotation=90, fontsize=fontsi-4, ha='center', va='center')
            ax.vlines(x=x+0.5, ymin=-20, ymax=30, color='black', lw=1, zorder=20)
        
        if k==0:
            auto_x= [2]
            for x in auto_x:
                #ax.axvspan(x-1.5,x-0.5,facecolor='#d2e5c9', alpha=0.45,zorder=-200)
                #ax.axvspan(x-0.5, x+1.5, facecolor='#eeeeee', alpha=0.45,zorder=-200)
                ax.text(x-0.43,ymax*0.95,"Constrained", fontsize=fontsi, weight = 'bold')
                
                loose_dashed="-."
                # Report relative increase of GHGs
                #rel_change=int(round( ((list_total_ghgs[x-1]-list_total_ghgs[x])/list_total_ghgs[x]*100),0) )
                # Add lines for cost cahnge with BAU for -90%
                #ax.vlines(x=x+0.133, ymin=list_total_ghgs[x], ymax=list_total_ghgs[x]*1.5, colors='black', linestyle=loose_dashed, lw=3,zorder=100)  
                #ax.hlines(y=list_total_ghgs[x]*1.5, xmin=x-0.83, xmax=x+0.133, colors='black', linestyle=loose_dashed, lw=3,zorder=100)
                #ax.vlines(x=x-0.85, ymin=list_total_ghgs[x-1], ymax=list_total_ghgs[x]*1.5, colors='black', linestyle=loose_dashed, lw=3,zorder=100)  
                #ax.annotate("{}%".format(rel_change), (x-0.7, list_total_ghgs[x]*1.52), fontsize=fontsi-2, weight='bold',zorder=100)  
                    
                rel_change=int(round( ((list_total_costs[x-1]-list_total_costs[x])/list_total_costs[x]*100),0) )
                # Add lines for cost cahnge with BAU for -90%
                ax.vlines(x=x-0.16, ymin=list_total_costs[x], ymax=list_total_costs[x]*1.5, colors='black', linestyle=loose_dashed, lw=3,zorder=100)  
                ax.hlines(y=list_total_costs[x]*1.5, xmin=x-1.14, xmax=x-0.16, colors='black', linestyle=loose_dashed, lw=3,zorder=100)
                ax.vlines(x=x-1.16, ymin=list_total_costs[x-1], ymax=list_total_costs[x]*1.5, colors='black', linestyle=loose_dashed, lw=3,zorder=100)  
                ax.annotate("{}%".format(rel_change), (x-0.7, list_total_costs[x]*1.52), fontsize=fontsi-2, weight='bold',zorder=100) 
        else:
            list_xs=[2,4]
            for j, x in enumerate(list_xs):
                #ax.axvspan(x-1.5,x-0.5,facecolor='#d2e5c9', alpha=0.45,zorder=-200)
                #ax.axvspan(x-0.5, x+.5, facecolor='#eeeeee' if x ==2 else "white", alpha=0.45,zorder=-200)
                ax.text(x-0.43,ymax*0.95,"Constrained" if x ==2 else "", fontsize=fontsi, weight = 'bold')
                loose_dashed="-."
                if x==2:
                    # Add lines for cost cahnge with BAU for -90%
                    #ax.vlines(x=x+0.132, ymin=list_total_ghgs[x], ymax=list_total_ghgs[x]*1.5, colors='black', linestyle=loose_dashed, lw=3,zorder=100)  
                    #ax.hlines(y=list_total_ghgs[x]*1.5, xmin=x-0.83, xmax=x+0.1225, colors='black', linestyle=loose_dashed, lw=3,zorder=100)
                    #ax.vlines(x=x-0.86, ymin=list_total_ghgs[x-1], ymax=list_total_ghgs[x]*1.5, colors='black', linestyle=loose_dashed, lw=3,zorder=100)  
                    # Report relative increase of GHGs
                    #ax.annotate("{}%".format(int(round( ((list_total_ghgs[x-1]-list_total_ghgs[x])/list_total_ghgs[x]*100),0) )), (x-0.7, list_total_ghgs[x]*1.52), fontsize=fontsi-2, weight='bold')  
                    
                    rel_change=int(round( ((list_total_costs[x-1]-list_total_costs[x])/list_total_costs[x]*100),0) )
                    # Add lines for cost cahnge with BAU for -90%
                    ax.vlines(x=x-0.17, ymin=list_total_costs[x], ymax=list_total_costs[x]*1.5, colors='black', linestyle=loose_dashed, lw=3,zorder=100)  
                    ax.hlines(y=list_total_costs[x]*1.5, xmin=x-1.14, xmax=x-0.16, colors='black', linestyle=loose_dashed, lw=3,zorder=100)
                    ax.vlines(x=x-1.15, ymin=list_total_costs[x-1], ymax=list_total_costs[x]*1.5, colors='black', linestyle=loose_dashed, lw=3,zorder=100)  
                    ax.annotate("{}%".format(rel_change), (x-0.7, list_total_costs[x]*1.5), fontsize=fontsi-2, weight='bold',zorder=100)        

            auto_x= 4
            #ax.axvspan(auto_x-0.5,auto_x+1.5, facecolor='lightgrey', alpha=0.45,zorder=-200)
            ax.text(auto_x-0.3,ymax*0.95,"Off-Grid", fontsize=fontsi, weight = 'bold')

        fig.canvas.draw()
        labels_2 = [item.get_text() for item in ax.get_xticklabels()]
        ax.set_xticklabels(labels_2, weight = 'bold', fontsize=fontsi-4)
        ax.tick_params(axis='x', pad=140)

        # Plot the total values for net cost and GHG emissions
        total_values_c = df_select_cost_2.sum(axis=0)
        ax.plot([i - 0.155 for i in range(len(total_values_c))], total_values_c, marker='D', 
                markerfacecolor='white', markeredgewidth=3.5, markeredgecolor='black',
                linestyle='None', markersize=23, zorder=200)
        total_values_g = df_init_gwp_2.sum(axis=0)
        ax.plot([i + 0.142 for i in range(len(total_values_g))], total_values_g, marker='D', 
                markerfacecolor='white', markeredgewidth=3.5, markeredgecolor='black',
                linestyle='None', markersize=23, zorder=200)
        
        ax.set_zorder(ax2.get_zorder()+1)
        ax.patch.set_visible(False)
        
    fig.subplots_adjust(wspace=0.22)
    plt.savefig(r"figs/contri_costs_ghgs.png", dpi = 300, bbox_inches = 'tight')

#generate_cost_ghg_figure(df_select_cost, df_init_gwp, df_init_cost)

def plot_lca_impact_category(lca_method, df_init_lca, nr_scenarios=9, nr_left_subplot=3):
    """
    Generate a dual-subplot bar plot for visualizing environmental impacts in a life cycle assessment (LCA).

    Parameters:
    - lca_method (tuple): A tuple with different elements.
    - df_init_lca (DataFrame): The input DataFrame with LCA results.
    - nr_scenarios (int, optional): Number of scenarios to display in the plot. Default is 8.
    - nr_left_subplot (int, optional): Number of scenarios to display on the left subplot. Default is 3.
    
    Output:
    - The generated plot is saved as a PNG file in the "figs" directory, with the filename derived from the LCA method and impact category.
    - The plot is adjusted for layout, and the resulting image is saved.

    Example Usage:
    plot_lca_impact_category(lca_method, df_init_lca, nr_scenarios=9, nr_left_subplot=3)
    """
    rc = {
        'axes.facecolor':'white',
        'axes.grid' : False,
        'grid.color': '.8',
        'font.family':'Arial',
        'font.size' : 8,
        'xtick.bottom': True,
        'xtick.top': True,
        'ytick.left': True,
        'ytick.right': True,
            'patch.facecolor': 'b',
            'patch.linewidth': 1,
        'patch.force_edgecolor': True,
        'axes.linewidth':2
        }
    plt.rcParams.update(rc) 
    plt.rcParams['axes.axisbelow'] = True
    plt.rcParams.update({'mathtext.default':  'regular' })
    plt.rcParams["axes.edgecolor"] = "black"
    fontsi=45

    fig, axes = plt.subplots(1, 2, figsize=(60, 20), gridspec_kw={'width_ratios': [nr_left_subplot/nr_scenarios, (nr_scenarios-nr_left_subplot)/nr_scenarios]})
    df_init_sel = df_init_lca.xs(lca_method[1], level='category', axis=0, drop_level=False).droplevel([0, 2, 3], axis=0)
    substr = "(BI)"

    for k, ax in enumerate(axes):
        df_init_sel_2 = df_init_sel.filter(regex=substr) if k == 0 else df_init_sel.drop(df_init_sel.filter(regex=substr).columns, axis=1)
        fontsi = 35

        df_init_sel_2.T.plot(kind='bar', ax=ax, width=0.40, stacked=True, zorder=9, legend=True if k == 1 else False,
                            color=[color_dict_env[r] for r in df_init_sel_2.T.columns])

        bars = ax.patches
        patterns = [hex_dict_env[z] for z in df_init_sel_2.T.columns]
        hatches = [h for h in patterns for _ in range(int(len(bars) / len(patterns)))]

        for bar, hatch in zip(bars, hatches):
            bar.set_hatch(hatch)  

        if k == 1:
            ax.legend(labels=list(df_init_sel.T.columns), fontsize=fontsi, ncol=5, bbox_to_anchor=(0.95, -0.05), frameon=False)

        ylabel = "ReCiPe Endpoint [{}/year]".format(dict_units[lca_method]) if lca_method[1] == 'total' else "Impacts on {} [{}/year]".format(lca_method[1].replace(": ", " "), dict_units[lca_method])
        ax.set_ylabel(ylabel, size=fontsi)
        ax.xaxis.grid(False, which='both')
        ax.yaxis.grid(True, which='both')
        ax.set_xlim(-0.41, len(list(df_init_sel_2.columns))-0.62)
        ax.set_ylim(df_init_sel_2.min(axis=1).min() * 2, df_init_sel_2.sum(axis=0).max()*1.1)

        ymax= df_init_sel_2.applymap(lambda x: max(x, 0)).abs().sum().sum()
        ymin= df_init_sel_2.applymap(lambda x: min(x, 0)).sum().sum()
        ax.set_ylim(ymin*1.02,ymax*1.02)

        ax.set_title("a. Industry (bakery)" if k == 0 else "b. Entire multi-energy system", fontsize=fontsi+8, weight='bold')

        ax.hlines(y=0, xmin=-1, xmax=30, color='black', lw=2, zorder=20)

        ax.tick_params(axis='both', which='major', labelsize=fontsi, rotation=360)
        ax.tick_params(axis='both', which='minor', labelsize=fontsi-3, rotation=47)

        # Increase the font size of the exponent part for the subplot
        ax.yaxis.offsetText.set(size=32) 

    to_disk = "figs/{}.png".format(lca_method[1].replace("/", "_").replace(": ", "_").replace(" ", "_").replace(",", ""))
    fig.subplots_adjust(wspace=0.12)
    plt.savefig(to_disk, dpi=150, bbox_inches='tight')

def generate_radar_plot_lca(df_init_lca, df_init_cost, nr_rows = 3, nr_cols = 2):
    """
    Generate radar plots based on input dataframes.

    Parameters:
        df_init_lca (pd.DataFrame): DataFrame for LCA impacts.
        df_init_cost (pd.DataFrame): DataFrame for cost results.

    Returns:
        None
    """
    # Preprocess data for spider graphs
    total_df_spider = df_init_lca.groupby(['category']).sum().T.rename(columns=column_dict_categories).T
    total_df_spider['maximum'] = np.array(total_df_spider.max(axis=1))
    df_normalized = total_df_spider.div(total_df_spider['maximum'], axis=0).drop(columns=['maximum'])
    df_normalized[df_normalized<0]=0 # This is done to avoid strange results in the spider graph, but can be minus due to avoidance
    df_values = df_normalized.T.to_numpy()

    # Set plotting configurations
    rc = {
        'axes.facecolor': 'white',
        'axes.grid': True,
        'grid.color': '.8',
        'font.family': 'Arial',
        'font.size': 90,
        'xtick.bottom': True,
        'xtick.top': True,
        'ytick.left': True,
        'ytick.right': True,
        'patch.facecolor': 'b',
        'patch.linewidth': 11,
        'patch.force_edgecolor': True,
        'axes.linewidth': 20
    }

    plt.rcParams.update(rc)
    plt.rcParams['axes.axisbelow'] = True
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 7
    plt.rcParams['grid.color'] = '0.3'
    plt.rcParams['patch.linewidth'] = 25

    # Collect data for radar plots
    labels_2 = list(df_normalized.T.columns)
    names = list(df_normalized.columns)
    list_configs = labels_2

    data = [labels_2, ('', df_values)]
    colors = ["darkblue"] * 10

    N = len(data[0])
    theta = radar_factory(N, frame='polygon')
    list_total_costs = list(df_init_cost.T.total_costs / 1e3)

    tickss = ['0', '0.2', '0.4', '0.6', "0.8", '1']
    fig, axs = plt.subplots(nr_cols, nr_rows, figsize=(150,104), subplot_kw=dict(projection='radar'))
    general_font = 170

    # Results per config
    for i, cf in enumerate(names):
        j = i // nr_rows  # Divide by 4 to determine the row index
        k = i % nr_rows  # Remainder when divided by 4 to determine the column index

        data = [labels_2, (names, [df_values[i]])]

        N = len(data[0])
        theta = radar_factory(N, frame='polygon')

        spoke_labels = data.pop(0)
        title, case_data = data[0]

        axs[j, k].text(0.8, 1.60, "{}".format(names[i]), fontsize=general_font + 5, weight='bold')
        axs[j, k].text(0.87, 1.57, "{} M€".format(round(list_total_costs[i], 2)), fontsize=general_font - 5, style='italic')
        axs[j, k].set_ylim(bottom=0)

        for d in case_data:
            line = axs[j, k].plot(theta, d, linewidth=8, color=colors[i])
            axs[j, k].fill(theta, d, alpha=0.2, color=colors[i])

        axs[j, k].set_varlabels("")
        axs[j, k].set_rgrids([0, 0.2, 0.4, 0.6, 0.8, 1], alpha=0.8, angle=15, labels=tickss, fontsize=general_font - 5)
        axs[j, k].grid(linewidth=2)

        #if i > 2:
            #axs[j, k].set_facecolor('#eeeeee')
            #axs[j, k].set_alpha(0.45)

        axs[j, k].set_frame_on(True)
        pi = np.pi
        for b in range(N):
            angle_rad = b / float(N) * 2 * pi

            if angle_rad == 0:
                ha, distance_ax = "center", 1.05
            elif 0 < angle_rad < pi:
                ha, distance_ax = "left", 1.2
            elif angle_rad == pi:
                ha, distance_ax = "center", 1.06
            else:
                ha, distance_ax = "right", 1.2
            axs[j, k].text(angle_rad, distance_ax, labels_2[b], size=general_font - 19,
                           horizontalalignment=ha, verticalalignment="center")

    # Draw horizontal lines at specific coordinates
    ys = [0.497]
    for y in ys:
        line = plt.Line2D([0, 0.99], [y, y], transform=fig.transFigure, color='black', lw=5)
        fig.add_artist(line)

    plt.rc('xtick', labelsize=general_font)
    plt.rc('ytick', labelsize=general_font)
    fig.subplots_adjust(top=0.01, bottom=0, wspace=8.5, hspace=13.5)

    plt.tight_layout()
    plt.savefig("figs/Results_trade_offs_multi.png", dpi=75, bbox_inches='tight')

def generate_ind_radar_plot_lca(col_name, df_init_lca):
    """
    Generate radar plot based on input dataframes.

    Parameters:
        col_name (str): column name of which spider graph to be made.
        df_init_lca (pd.DataFrame): DataFrame for LCA impacts.

    Returns:
        None
    """
    rc = {
        'axes.facecolor':'white',
        'axes.grid' : True,
        'grid.color': '.8',
        'font.family':'Arial',
        'font.size' : 90,
        'xtick.bottom': True,
        'xtick.top': True,
        'ytick.left': True,
        'ytick.right': True,
            'patch.facecolor': 'b',
            'patch.linewidth': 11,
        'patch.force_edgecolor': True,
        'axes.linewidth':20
        }

    plt.rcParams.update(rc) 
    plt.rcParams['axes.axisbelow'] = True
    plt.rcParams.update({'mathtext.default':  'regular' })
    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 7
    plt.rcParams['grid.color'] = '0.3'
    plt.rcParams['patch.linewidth'] = 15

    # Preprocess data for spider graphs
    total_df_spider = df_init_lca.groupby(['category']).sum().T.rename(columns=column_dict_categories).T
    total_df_spider['maximum'] = np.array(total_df_spider.max(axis=1))
    df_normalized = total_df_spider.div(total_df_spider['maximum'], axis=0).drop(columns=['maximum'])[col_name]
    df_normalized[df_normalized<0]=0 # This is done to avoid strange results in the spider graph, but can be minus due to avoidance
    df_values = df_normalized.T.to_numpy()

    # Iterate over data
    fig, ax = plt.subplots(1, 1, figsize=(15, 10), subplot_kw=dict(projection='radar'))
    general_font = 22

    labels = list(df_normalized.index)
    names = list([col_name])

    data = [labels,
            (names, [df_values])
            ]

    N = len(data[0])
    theta = radar_factory(N, frame='polygon')

    spoke_labels = data.pop(0)
    title, case_data = data[0]
    ax.set_ylim(bottom=0)
    ax.set_rgrids([], alpha=0.8, angle=15, labels="", fontsize=general_font - 5)

    for d in case_data:
        line = ax.plot(theta, d, linewidth=8, color='darkblue')
        ax.fill(theta, d, alpha=0.2, color='darkblue')

    ax.set_varlabels("")
    ax.grid(linewidth=1)
    ax.set_frame_on(True)
    pi = np.pi
    for b in range(N):
        angle_rad = b / float(N) * 2 * pi
        if angle_rad == 0:
            ha, distance_ax = "center", 1.1
        elif 0 < angle_rad < pi:
            ha, distance_ax = "left", 1.2
        elif angle_rad == pi:
            ha, distance_ax = "center", 1.1
        else:
            ha, distance_ax = "right", 1.2

    plt.rc('xtick', labelsize=general_font)
    plt.rc('ytick', labelsize=general_font)
    fig.subplots_adjust(top=0.27, bottom=0.15, wspace=6.2, hspace=0.8)
    plt.tight_layout()
    plt.savefig("figs/Results_trade_offs_{}.png".format(col_name.replace("\n","")), dpi = 25, bbox_inches = 'tight')

# Plot additional images
def plot_image(x_coord, y_coord, image, ax=None):
    ax = ax or plt.gca()
    im = OffsetImage(image, zoom=45 / ax.figure.dpi)
    im.image.axes = ax
    ab = AnnotationBbox(im, (x_coord, y_coord), frameon=False, pad=0.0, zorder=-10)
    ax.add_artist(ab)

def plot_pareto_front(df_init_out_gr, add_sep_points=False, name_sep_points=False):
    """
    Plots trade-offs between annual life cycle GHG emissions and costs.

    Parameters:
    - df_init_out_gr (pd.DataFrame): DataFrame containing data for plotting.
    - add_sep_point (pd.DataFrame): DataFrame containing additional data for point showed on the fig.
    - name_sep_point (pd.DataFrame): Name of additional data for point showed on the fig.

    Returns:
    - None
    """

    # Extract data
    x = list(df_init_out_gr.total_ghg / 1e3)
    y = list(df_init_out_gr.total_costs / 1e3)

    # Reference values
    ref_ghgs = df_init_out_gr.ghg_init.max() / 1e3
    ref_costs = df_init_out_gr.cost_init.max() / 1e3

    # Create a plot with specified style
    rc = {
        'axes.facecolor': 'white',
        'axes.grid': True,
        'grid.color': '.8',
        'font.family': 'Arial',
        'font.size': 70,
        'xtick.bottom': False,
        'xtick.top': False,
        'ytick.left': False,
        'ytick.right': False,
        'patch.facecolor': 'b',
        'patch.linewidth': 3,
        'patch.force_edgecolor': True,
        'axes.linewidth': 3
    }
    plt.rcParams.update(rc)
    fig, ax = plt.subplots(figsize=(20, 10))

    loose_dashed = ((0, (5, 10)))

    # Plot trade-off curve
    x2 = np.linspace(x[0], x[-1], 11)
    y2 = pchip_interpolate(x, y, x2)
    ax.plot(x2, y2, linewidth=4, marker="o", markersize=15, color='darkblue', fillstyle="none", zorder=10, mew=1.3)

    # Set plot limits and tick parameters
    ax.set_ylim(0, max(max(y),ref_costs)*1.2)
    ax.set_xlim(min(min(x),ref_ghgs)*1.2, max(max(x),ref_ghgs)*1.2)
    font = 26
    ax.tick_params(axis='both', which='major', labelsize=font, rotation=0)

    if min(min(x),ref_ghgs)<0:
        # add zero point for ghgs on the x-axis.
        ax.vlines(x=0, ymin=0, ymax=max(max(y),ref_costs)*1.2, colors='black', linestyle="--", lw=1.2)

    # Add additional point if provided:
    if not isinstance(add_sep_points, bool):
        for i, add_sep_point in enumerate(add_sep_points):
            ax.scatter(add_sep_point.total_ghg.item()/1e3, add_sep_point.total_costs.item()/1e3, marker='o', s=210, color='darkgreen', edgecolor='darkblue', lw=1.5)
            ax.annotate(name_sep_points[i], (add_sep_point.total_ghg.item()/1e3 - 0.02, add_sep_point.total_costs.item()/1e3 + 0.05), fontsize=font - 5)
            if os.path.exists(r"figs\Results_trade_offs_{}.png".format(name_sep_points[i])):
                image = plt.imread(r"figs\Results_trade_offs_{}.png".format(name_sep_points[i]) )
                plot_image(add_sep_point.total_ghg.item()/1e3, add_sep_point.total_costs.item()/1e3, image, ax=ax)

    # Annotate reference points
    ax.scatter(ref_ghgs, ref_costs, marker='o', s=210, color='red', edgecolor='darkblue', lw=1.5)
    ax.annotate("BAU", (ref_ghgs - 0.02, ref_costs + 0.05), fontsize=font - 7)

    # Annotate trade-off points
    texts = ["GHG-Min", "GHG-80", "GHG-50", "Cost-Min"]
    ax.scatter(x, y, marker='o', s=100, color='red')
    for i, txt in enumerate(texts):
        ax.annotate(txt, (x[i] - 0.02, y[i] + 0.05), fontsize=font - 5)

    # Annotate percentage changes with BAU
    for i in range(0,len(x)):
        ax.vlines(x=ref_ghgs, ymin=ref_costs, ymax=y[i], colors='black', linestyle=loose_dashed, lw=0.8)
        if abs(((y[i] - ref_costs) / ref_costs * 100))>=1:
            ax.annotate("{}%".format(int(round(((y[i] - ref_costs) / ref_costs * 100), 0))),
                        (ref_ghgs - 0.05, y[i] - 0.05), fontsize=font - 7)
        ax.hlines(y=y[i], xmin=ref_ghgs, xmax=x[i], colors='black', linestyle=loose_dashed, lw=0.8)
        if abs(((x[i] - ref_ghgs) / ref_ghgs * 100))>=1:
            ax.annotate("{}%".format(int(round(((x[i] - ref_ghgs) / ref_ghgs * 100), 0))),
                        (((x[i] + ref_ghgs) / 1.4) - 0.05, y[i] - 0.05), fontsize=font - 7)

    # Set axis labels
    ax.set_ylabel("Annual costs [M€/year]", fontsize=font)
    ax.set_xlabel("Annual life cycle GHG emissions [kt CO$_{2}$-eq./year]", fontsize=font)

    paths = [
        r"figs\Results_trade_offs_GHG-Min.png",
        r"figs\Results_trade_offs_Cost-Min.png",
    ]

    x_wo = list(df_init_out_gr.total_ghg/1e3)[:1] + list(df_init_out_gr.total_ghg/1e3)[3:]
    y_wo = list(df_init_out_gr.total_costs/1e3)[:1] + list(df_init_out_gr.total_costs/1e3)[3:]

    for i, pt in enumerate(paths):
        if (pt!="") and (os.path.exists(pt)):
            image = plt.imread(pt)
            plot_image(x_wo[i], y_wo[i], image, ax=ax)
        
    if os.path.exists(r"figs\Results_trade_offs_bau.png"):
        image = plt.imread(r"figs\Results_trade_offs_bau.png")
        plot_image(ref_ghgs, ref_costs, image, ax=ax)

    plt.savefig("figs/Pareto_{}.png".format(ASSESS_COUNTRY), dpi = 200, bbox_inches = 'tight')

def plot_sensitivity_results_param_inc(df_info_increased, parameter_name, technologies, max_nr=30, step_inc=400):
    """
    Plot sensitivity results for cost increase.

    Parameters:
    - df_info_increased (pd.DataFrame): Results of sensitivity analysis for increased parameter.
    - parameter_name (str): Name of the parameter being increased (e.g., 'grid_capex').
    - technologies (list): List of technologies to be assessed (e.g., ['cap_wind_on', 'cap_pv', 'p_grid_connection', 'cap_bat_en']).
    - max_nr (int): Maximum number of sensitivity iterations.
    - step_inc (int): Step increase of the parameter per iteration.

    This function generates a plot displaying the sensitivity results for cost increase.
    It visualizes the installed capacity of different energy sources, the share of the specified parameter,
    and the assumed value of the parameter.

    """
    rc = {
        'figure.figsize': (10, 4),
        'axes.facecolor': 'white',
        'axes.grid': True,
        'grid.color': '.9',
        'font.family': 'Arial',
        'font.size': 12,
        'xtick.bottom': True,
        'xtick.top': True,
        'ytick.left': True,
        'ytick.right': True,
        'patch.facecolor': 'b',
        'patch.linewidth': 3,
        'patch.force_edgecolor': False
    }
    plt.rcParams.update(rc)
    plt.rcParams['axes.axisbelow'] = True
    plt.rcParams.update({'mathtext.default': 'regular'})
    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 1

    # Plot
    fig, ax = plt.subplots(figsize=(7, 4))

    colors = ["#495100", "#f1c232", "#600000", "red", "#014940", "darkgrey", "darkorange", "silver", "#38628c", "#ff9992"]
    font = 14
    lw = 2

    x2 = np.linspace(list(df_info_increased.index)[0], list(df_info_increased.index)[-1], 400)

    for tech, color in zip(technologies, colors):
        y2 = pchip_interpolate(list(df_info_increased.index), list(df_info_increased[tech]), x2)
        ax.plot(x2, y2, lw=lw, color=color, label=tech)

    ax.grid(False)

    ax.set_xlabel(f"{parameter_name.capitalize()} [euro/kW]".replace("_"," "), size=font)
    ax.set_ylabel("Installed capacity [MW(h)]", size=font)
    ax.xaxis.grid(False, which='both')
    ax.yaxis.grid(True, which='both')
    ax.set_xlim(0, (max_nr - 1) * step_inc)
    ax.set_ylim(0, round(1.15*np.ceil(df_info_increased[technologies].max().max())))

    ax.tick_params(axis='both', which='major', labelsize=font, rotation=360)
    ax.tick_params(axis='both', which='minor', labelsize=font, rotation=47)

    ax2 = ax.twinx()
    x2 = np.linspace(list(df_info_increased.index)[0], list(df_info_increased.index)[-1], 250)
    y2 = pchip_interpolate(list(df_info_increased.index), list(df_info_increased['share_grid_inv']), x2)
    ax2.plot(x2, y2, lw=lw, label=f"Share {parameter_name}".replace("_"," "), color='grey')

    ax2.fill_between(x2, y2, 0, color='grey', alpha=0.3, label='Positive area')
    ax2.set_ylabel(f"Share overall {parameter_name} [%]".replace("_"," "), size=font, color='darkgrey')

    ax2.tick_params(axis='both', which='major', labelsize=font, rotation=360)
    ax2.tick_params(axis='both', which='minor', labelsize=font, rotation=47)

    ax2.xaxis.grid(False, which='both')
    ax2.yaxis.grid(False, which='both')
    ax2.set_xlim(0, (max_nr - 1) * step_inc)
    ax2.set_ylim(0, round(1.2*np.ceil(df_info_increased[technologies].max().max())))
    
    y2 = pchip_interpolate(list(df_info_increased.index), list(df_info_increased['total_costs']/1e3), x2) 
    ax.plot(x2, y2, lw=lw, color='darkblue', label="total_costs")

    handles, labels = ax.get_legend_handles_labels()
    new_labels = [cap_mapping[label] for label in labels]
    print(new_labels)

    # To represent the assumption used:
    ax.vlines(x=COST_DATA[NAME_REF_DB].grid_capex, ymin=0, ymax=round(1.15*np.ceil(df_info_increased[technologies].max().max())), colors='black', lw=1, linestyle="--", \
        label=f"{parameter_name.capitalize()} assumed".replace("_"," "))

    ax.legend(handles, new_labels, fontsize=font - 3, ncol=2, frameon=False, loc='best')

    plt.savefig(f"figs/sens_{parameter_name}_increase.png", dpi=250, bbox_inches='tight')

# Assuming df_sens_op_grid_capex_gr is your DataFrame with sensitivity results
# and 'grid_capex' is the parameter being increased
# and ['cap_wind_on', 'cap_pv', 'p_grid_connection', 'cap_bat_en'] are the technologies assessed.
#plot_sensitivity_results(df_sens_op_grid_capex_gr, 'grid_capex', ['cap_wind_on', 'cap_pv', 'p_grid_connection', 'cap_bat_en'])

def plot_sensitivity_results(df_sens_min_gr, df_sens_out_gr, figsize=(5, 2.5), color=None):
    """
    Plot sensitivity analysis results for annual system costs.

    Parameters:
    - df_sens_min_gr (pd.DataFrame): DataFrame containing the results of the sensitivity analysis for minimum values.
    - df_sens_out_gr (pd.DataFrame): DataFrame containing the results of the sensitivity analysis for output values.
    - figsize (tuple, optional): Size of the resulting plot. Default is (5, 2.5).
    - color (list, optional): List of colors for the bar plot. If None, a default color scheme is used.

    Returns:
    - None: The function generates and saves the plot as 'figs/sens_parameters.png'.

    Note:
    - The input DataFrames should have columns like 'Costs [$\Delta$%]' for plotting.

    Example:
    plot_sensitivity_results(df_sens_min_gr, df_sens_out_gr)

    The resulting plot illustrates the sensitivity analysis for annual system costs, showing the impact of different
    parameters on the overall system costs.

    """
    rc = {
        'axes.facecolor': 'white',
        'axes.grid': False,
        'grid.color': '.8',
        'font.family': 'Arial',
        'font.size': 20,
        'xtick.bottom': False,
        'xtick.top': False,
        'ytick.left': False,
        'ytick.right': False,
        'patch.facecolor': 'white',
        'patch.linewidth': 1,
        'patch.force_edgecolor': True,
        'axes.linewidth': 1
    }

    fontiz = 15
    plt.rcParams.update(rc)

    fig, ax = plt.subplots(figsize=figsize)

    if color is None:
        color = ["indigo", "#f1c232", "red", "green", "azure", "darkgrey", "darkorange", "yellow", "blue",
                 "#ff9992", "lightgrey", "black", "blue"]

    df_sens_min_gr["Costs [$\Delta$%]"].plot(kind='barh', ax=ax, position=0.5, width=1, zorder=9, legend=False,
                                             color=color)
    # Plot mins
    df_sens_out_gr["Costs [$\Delta$%]"].plot(kind='barh', ax=ax, position=0.5, width=1, zorder=9, legend=False,
                                             color=color, hatch='//---')

    ax.grid(True)
    ax.set_xlabel("Annual system costs [$\Delta$%]", size=fontiz - 2)
    ax.xaxis.grid(False, which='both')
    ax.yaxis.grid(False, which='both')

    # Dynamically set x and y-axis limits based on min and max values in the DataFrame
    min_value = min(df_sens_min_gr["Costs [$\Delta$%]"].min(), df_sens_out_gr["Costs [$\Delta$%]"].min())
    max_value = max(df_sens_min_gr["Costs [$\Delta$%]"].max(), df_sens_out_gr["Costs [$\Delta$%]"].max())
    ax.set_xlim(min_value - 1, max_value + 1)

    min_value_y = 0
    max_value_y = max(len(df_sens_min_gr.index), len(df_sens_out_gr.index))
    ax.set_ylim(min_value_y - 0.6, max_value_y + 0.6-1 )

    ax.tick_params(axis='both', which='major', labelsize=fontiz - 3, rotation=360)
    ax.tick_params(axis='both', which='minor', labelsize=fontiz - 3, rotation=47)

    plt.savefig("figs/sens_parameters.png", dpi=250, bbox_inches='tight')

def create_heatmap(df):
    """
    Create a heatmap from a DataFrame and save it as an image.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data to be visualized.

    Returns:
    - None

    Example:
    create_heatmap(my_dataframe)
    """

    ref_name = "opt_results_{}_1_0".format(ASSESSMENT_YEAR)

    # Define colors for the custom colormap
    colors = [(0, 'darkgreen'), (0.5, 'white'), (1, 'darkred')]

    # Create the custom colormap
    custom_cmap = LinearSegmentedColormap.from_list('custom_RdYlGn', colors)

    # Data preprocessing
    df_init = df.copy()
    df = df.div(df[ref_name], axis=0).dropna()

    # This is to get the initial capacities
    df_init = df_init.T[list(df.index)].T

    # Rename columns for later
    df.rename(index=cap_mapping, inplace=True)
    df.columns = ['Off-grid (ref.)'] + list(df.columns[1:])

    # Plotting setup
    vmax_set=2 #df.max().max()
    vmin_set=0 #df.min().min()
    fig, ax = plt.subplots(figsize=(10, 8))
    norm = SymLogNorm(linthresh=2, linscale=2, vmin=vmin_set, vmax=vmax_set)
    cax = ax.matshow(df, cmap=custom_cmap, norm=norm)
    cbar = fig.colorbar(cax, format="%.2f", extend='both', shrink=0.5)
    cbar.set_label('Relative change [-]', fontsize=15)

    # Customize colorbar ticks and labels
    cbar_ticks = np.linspace(vmin_set, vmax_set, num=3)
    cbar_ticks = np.insert(cbar_ticks, 1, 1)
    cbar.set_ticks(cbar_ticks)
    decimal_places = 1
    cbar.set_ticklabels([f'{tick:.{decimal_places}f}' for tick in cbar_ticks], fontsize=15)

    # Annotate the plot with real values
    for i in range(len(df_init.index)):
        for j in range(len(df_init.columns)):
            value = df_init.iloc[i, j] 
            if df_init.index[i] == 'total_costs':
                value = value/1e3 # Round annual costs from thousands to Meuro.
            ax.text(j, i, "{}".format(round(value,1) if value>10 else round(value,2)), ha='center', va='center', color='black', fontsize=8)

    # Set labels and title
    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels(df.columns, rotation=45, ha='left', fontsize=15)
    ax.set_yticks(range(len(df.index)))
    ax.set_yticklabels(df.index, fontsize=15)

    plt.xlabel('Scenarios (-50% Capex)', fontsize=15)
    plt.ylabel('Capacities of energy technologies [MW(h)],\n cost [M€/a], and curtailment [–]', fontsize=15)

    # Save the plot as an image
    plt.savefig("figs/output_ext.png", dpi=250, bbox_inches='tight')

# Example usage:
# create_heatmap(your_dataframe)
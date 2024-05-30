import gurobipy as gp
import brightway2 as bw
import pandas as pd 
import numpy as np
import uuid
from functools import partial
from bw2io.strategies import add_database_name, csv_restore_tuples
from bw2io.importers.base_lci import LCIImporter
from collections import defaultdict
import time

import energy_data_processor as ep
from config import (CC_METHOD, NAME_REF_DB, CALC_ALL_LCA_IMPACTS, DB_NAME, MJ_KG_H2, ASSESSMENT_YEAR, DB_NAME_CONS,
                    PROJECT_NAME, DAYS, LATITUDE, LONGITUDE, FIT, LOC_ELECT, ASSESS_COUNTRY, COST_DATA, TOTAL_HH, AV_COP_ASHP) # import global vars
from mapping import my_methods, contribution_mapping_system  #import mappings

# Name for the MES to be created
IMPORTER = LCIImporter(DB_NAME) #

def optimize_mes(df_data, w_cost, w_env, parm, dict_ghg, dict_limits, sec_db = NAME_REF_DB,
                        credit_env_export=True, autonomous_elect=False, autonomous_gas=False, autonomous_balanced=False, export_results=True,
                        eps_ghg_constraint = False, euro_ton_co2 = 0,
                        consider_down_times = True, biomass_constrained = True,
                        high_t_heat = True, road_transport = True, micro_wt = False,
                        ghg_init=0, cost_init=0, heuristics=False, export_alias = "", res_heat_stor = False,
                        warm_start=False, save_for_warm_start=False,
                        mip_gap=0.005, int_feas_tol = 1e-7, 
                        time_limit = 10*3600, calc_all_lca_impacts=CALC_ALL_LCA_IMPACTS):
    
    """
    Solves a (single or multi-objective) mixed integer linear optimization problem for a multi-energy system (MES).

    Args:
        df_data (DataFrame): DataFrame with hourly input data.
        w_cost (float): Weight of cost objective (between 0 and 1).
        w_env (float): Weight of environmental objective (between 0 and 1).
        parm (dict): Dictionary with techno-economic data.
        dict_ghg (dict): Dictionary with environmental emission data.
        dict_limits (dict): Dictionary with limits of technologies and sizes.

        sec_db (str, optional): Database for calculating environmental results. Default is NAME_REF_DB.
        credit_env_export (bool, optional): Use environmental credit for exporting an energy commodity (electricity + hydrogen).
            Default is False.
        autonomous_elect (bool, optional): Optimize as an autonomous energy system without a power grid connection.
            Default is False.
        autonomous_gas (bool, optional): Optimize as an autonomous energy system without a gas grid connection.
            Default is False.
        autonomous_balanced (bool, optional): Optimize as an autonomous energy system with balanced energy autonomy, 
                meaning own renewable power generation >= than power demand.
            Default is False.
        export_results (bool, optional): Export results to a file (Excel). Default is True.
        eps_ghg_constraint (bool, optional): Add an epsilon constraint (float) on life cycle GHG emission objective.
            Default is False.
        euro_ton_co2 (float, optional): CO2 price in Euro per ton. Default is 0.
        consider_down_times (bool, optional): Consider up- and downtimes of specified technologies. Default is True.
        biomass_constrained (bool, optional): Constrain biomass availability. Default is True.
        high_t_heat (bool, optional): Consider high-temperature heat requirement in the MES. Default is True.
        road_transport (bool, optional): Consider residential road transportation in the MES. Default is True.
        micro_wt (bool, optional): Use hourly renewable energy generation with smaller wind turbines (lower capacity factor).
            Default is False.
        ghg_init (float, optional): Reference GHG emissions for comparison. Default is 0.
        cost_init (float, optional): Reference cost for comparison. Default is 0.
        heuristics (bool, optional): Apply other heuristics method if true. Default is False.
        export_alias (str, optional): Export alias for Excel export. Default is an empty string.
        res_heat_stor (bool, optional): Include residential heat storage unit if true. Default is False.
        warm_start (bool, optional): Use warm start from a previous optimization. Default is False.
        save_for_warm_start (bool, optional): Save the solution for warm start to disk. Default is False.
        time_limit (int, optional): Time limit of optimization in seconds. Default is 36000 (10 hours).
        calc_all_lca_impacts (bool, optional): whether to calculate other environmental impact categories. Default is True.

    Returns:
        overview_totals (DataFrame): DataFrame with one-dimensional key results.
        lca_results (DataFrame): DataFrame with LCA results.
        df_out (DataFrame): DataFrame with hourly decisions of technologies.
    """
    
    start = time.time()
    delta_t = 1
    
    # If no credit assumed for injection of grid electrciy, then replace that with zeros
    if credit_env_export == False:
        df_data.ghg_impact_cons = 0
        dict_ghg['ghg_imp_smr'] = 0
        
    T=len(df_data)    
    
    # If residential road transport, overwrite
    total_hh_r = TOTAL_HH if road_transport else 0
    
    # Decide to install micro wind turbines or not:
    wind_MW_array_on = df_data.wind_MW_array_on_micro if micro_wt else df_data.wind_MW_array_on
    
    """
    Step 1: Create a model and specify parameters of the model
    """

    m = gp.Model("obj")
    
    # Set parameters of model, tighten or relax
    m.setParam('MIPGap', mip_gap )
    m.setParam('IntFeasTol', int_feas_tol)
    m.setParam('TimeLimit', time_limit )
    m.setParam('Presolve', 2)
    m.setParam('MIPFocus', 1)
    
    if heuristics:
        m.setParam('NormAdjust', 2)
        m.setParam('Heuristics', 0.1)

    """
    Step 2: Define variables
    """
    # Grid variables
    p_grid_inj = m.addVars(T, name='p_grid_inj', ub = dict_limits['max_grid_cap'])
    p_grid_abs = m.addVars(T, name='p_grid_abs', ub = dict_limits['max_grid_cap'])
    bin_grid = m.addVars(T, vtype=gp.GRB.BINARY, name="bin_grid")

    # Battery variables
    p_battdis = m.addVars(T, name='p_battdis', ub=dict_limits['max_bat']) 
    p_battch = m.addVars(T, name='p_battch', ub=dict_limits['max_bat']) 
    E_batt = m.addVars(T, name="E_batt", ub=dict_limits['max_bat'])
    bin_bat = m.addVars(T, vtype=gp.GRB.BINARY, name='bin_bat') 

    # battery electric vehicle
    f_bev_1 = m.addVars(T, name='f_bev_1', ub=dict_limits['st_max'])
    f_bev_2 = m.addVars(T, name='f_bev_2', ub=dict_limits['st_max'])
    f_bev_3 = m.addVars(T, name='f_bev_3', ub=dict_limits['st_max'])

    # Hydrogen modelling
    E_h2_ves = m.addVars(T, name='E_h2_ves', ub=dict_limits['max_h2_storage'])
    p_h2_ves = m.addVars(T, name='p_h2_ves', lb=-dict_limits['max_h2_storage'], ub=dict_limits['max_h2_storage'])
    f_elect = m.addVars(T, name='f_elect', ub=dict_limits['st_max']) 
    p_elect = m.addVars(T, name="p_elect", ub=dict_limits['st_max'])
    p_fc = m.addVars(T, name='p_fc', ub=dict_limits['st_max'])
    p_fc_heat = m.addVars(T, name='p_fc_heat', ub=dict_limits['st_max'])
    f_fc = m.addVars(T, name='f_fc', ub=dict_limits['st_max'])

    # Wind electricity production
    pv_total = m.addVars(T, name="pv_total", ub=dict_limits['st_max'])
    wind_total_off = m.addVars(T, name="wind_total_off", ub=dict_limits['st_max'])
    wind_total_on = m.addVars(T, name="wind_total_on", ub=dict_limits['st_max'])

    #HP variables
    p_th_hp = m.addVars(T, name='p_th_hp', ub=dict_limits['st_max'])
    f_hp = m.addVars(T, name='f_hp', ub=dict_limits['st_max'])

    #LNG boiler
    p_d_boiler = m.addVars(T, name='p_d_boiler', ub=dict_limits['st_max'])
    f_d_boiler = m.addVars(T, name='f_d_boiler', ub=dict_limits['st_max'])

    # Heat storage medium
    if res_heat_stor:
        E_heat_stor = m.addVars(T, name='E_heat_stor', ub=dict_limits['st_max'])
        heat_stor_ch= m.addVars(T, name='heat_stor_ch', ub=dict_limits['st_max'])
        heat_stor_dis= m.addVars(T, name='heat_stor_dis', ub=dict_limits['st_max'])
        bin_h_stor = m.addVars(T, vtype=gp.GRB.BINARY, name='bin_h_stor') 

    # TS
    p_solar_heat = m.addVars(T, name='p_solar_heat', ub=dict_limits['st_max'])

    # innovative chp
    f_chp_h2 = m.addVars(T, name='f_chp_h2', ub=dict_limits['st_max'])
    p_th_chp_mix = m.addVars(T, name='p_th_chp_mix', ub=dict_limits['st_max'])
    p_chp_mix = m.addVars(T, name='p_chp_mix', ub=dict_limits['st_max'])
    
    if consider_down_times:
        #variables to model start-up and down times of mixed chp
        x_chp_mix = m.addVars(T, vtype=gp.GRB.BINARY, name="x_chp_mix")
        y_chp_mix = m.addVars(T, vtype=gp.GRB.BINARY, name="y_chp_mix")
        z_chp_mix = m.addVars(T, vtype=gp.GRB.BINARY, name="z_chp_mix")
        aux_chp_mix = m.addVars(T, name="aux_chp_mix", ub=dict_limits['st_max'])
    f_chp_mix = m.addVars(T, name='f_chp_mix', ub=dict_limits['st_max'])

    # binaries mixer
    f_chp_syngas = m.addVars(T, name="f_chp_syngas", ub=dict_limits['st_max']) 
    p_gasif = m.addVars(T, name="p_gasif", ub=dict_limits['st_max']) 
    f_gasif_wood = m.addVars(T, name="f_gasif_wood", ub=dict_limits['st_max']) 
    
    if consider_down_times:
        #variables to model start-up and down times of mixed wg
        x_wg = m.addVars(T, vtype=gp.GRB.BINARY, name="x_wg")
        y_wg = m.addVars(T, vtype=gp.GRB.BINARY, name="y_wg")
        z_wg = m.addVars(T, vtype=gp.GRB.BINARY, name="z_wg")
        aux_gasif = m.addVars(T, name="aux_gasif", ub=dict_limits['st_max'])

    # Export of hydrogen
    export_h2 = m.addVars(T, name="export_h2", ub=0 if (autonomous_balanced or autonomous_gas)
                                                        else dict_limits['st_max']) # no H2 export if energy autonomy

    #Residential NG boilers
    p_gb_res = m.addVars(T, name='p_gb_res', ub=dict_limits['st_max'])
    f_gb_res = m.addVars(T, name='f_gb_res', ub=dict_limits['st_max'])

    # Biogas from waste treatment plant
    f_biogas = m.addVars(T, name='f_biogas', ub=dict_limits['st_max'])
    f_chp_biogas = m.addVars(T, name='f_chp_biogas', ub=dict_limits['st_max'])
    
    # Biogas from waste treatment plant
    p_e_heating = m.addVars(T, name='p_e_heating', ub=dict_limits['st_max'])
    f_e_heating = m.addVars(T, name='f_e_heating', ub=dict_limits['st_max'])
    
    # e boiler
    p_e_boiler = m.addVars(T, name='p_e_boiler', ub=dict_limits['st_max'])
    f_e_boiler = m.addVars(T, name='f_e_boiler', ub=dict_limits['st_max'])

    ###### SINGLE VARS ######
    cap_wind_off = m.addVar(name="cap_wind_off", ub = dict_limits['max_wind_off'])
    cap_wind_on = m.addVar(name="cap_wind_on", ub = dict_limits['max_wind_on'])
    cap_pv = m.addVar(name="cap_pv", ub = dict_limits['max_pv'])
    cap_bat_en = m.addVar(name="cap_bat_en", ub = dict_limits['max_bat'])
    cap_bat_p = m.addVar(name="cap_bat_p", ub = dict_limits['max_bat'])
    cap_h2_ves = m.addVar(name="cap_h2_ves", ub = dict_limits['max_h2_storage']) 
    cap_electrolyzer = m.addVar(name="cap_electrolyzer", ub = dict_limits['max_electrolyzer']) 
    cap_fc = m.addVar(name="cap_fc", ub = dict_limits['max_fc']) 
    cap_hp = m.addVar(name="cap_hp", ub = dict_limits['max_hp'])
    cap_d_boiler = m.addVar(name="cap_d_boiler", ub = 0 if autonomous_gas else dict_limits['max_boiler'])
    cap_heat_storage = m.addVar(name="cap_heat_storage", ub = dict_limits['max_heat_storage'] if res_heat_stor else 0)
    cap_solar_heat = m.addVar(name="cap_solar_heat", ub = dict_limits['max_solar_heat_ha'])
    cap_chp_mix = m.addVar(name="cap_chp_mix", ub = dict_limits['max_chp_mix'])
    cap_gasif = m.addVar(name="cap_gasif", ub = dict_limits['max_gasif'])
    cap_gb_res = m.addVar(name="cap_gb_res", ub = 0 if autonomous_gas else dict_limits['max_boiler'])
    cap_e_heating = m.addVar(name="cap_e_heating", ub = dict_limits['st_max'])
    cap_e_boiler = m.addVar(name="cap_e_boiler", ub=dict_limits['st_max'])

    # to determine the selection of technologies
    bin_wind_off = m.addVar(vtype=gp.GRB.BINARY, name="bin_wind_off")
    bin_wind_on = m.addVar(vtype=gp.GRB.BINARY, name="bin_wind_on")
    bin_pv = m.addVar(vtype=gp.GRB.BINARY, name="bin_pv")
    bin_bat_en = m.addVar(vtype=gp.GRB.BINARY, name="bin_bat_en")
    bin_h2_ves = m.addVar(vtype=gp.GRB.BINARY, name="bin_h2_ves") 
    bin_electrolyzer = m.addVar(vtype=gp.GRB.BINARY, name="bin_electrolyzer") 
    bin_fc = m.addVar(vtype=gp.GRB.BINARY, name="bin_fc") 
    bin_hp = m.addVar(vtype=gp.GRB.BINARY, name="bin_hp")
    bin_boiler = m.addVar(vtype=gp.GRB.BINARY, name="bin_boiler")
    bin_heat_storage = m.addVar(vtype=gp.GRB.BINARY, name="bin_heat_storage")
    bin_solar_heat = m.addVar(vtype=gp.GRB.BINARY, name="bin_solar_heat")
    bin_chp_mix = m.addVar(vtype=gp.GRB.BINARY, name="bin_chp_mix")
    bin_gasif = m.addVar(vtype=gp.GRB.BINARY, name="bin_gasif")
    bin_gb_res = m.addVar(vtype=gp.GRB.BINARY, name="bin_gb_res")
    bin_e_boiler = m.addVar(vtype=gp.GRB.BINARY, name="bin_e_boiler")
    bin_e_heating = m.addVar(vtype=gp.GRB.BINARY, name="bin_e_heating")
    p_grid_connection = m.addVar(name="p_grid_connection", ub=0 if autonomous_elect else dict_limits['max_grid_cap']) # The grid connection is zero capacity in case there cannot per a grid power exchange
    share_bev = m.addVar(name="share_bev", lb=0 , ub=1)
    share_gasoline = m.addVar(name="share_gasoline", lb=0 , ub=1)
    
    # If warm start, read the solution from the file and set it as the initial solution for warm start
    if warm_start:
        m.read(warm_start)

    """
    Step 3: Add constraints
    """
    
    # Get total shares of residential road transporation
    m.addConstr(share_gasoline == 1 - share_bev)
    
    "Minimum island specific capacity constraints considering potential of renewables"
    m.addConstr(cap_pv >= dict_limits['min_cap'] * bin_pv)
    m.addConstr(cap_wind_on >= dict_limits['min_cap'] * bin_wind_on)
    m.addConstr(cap_wind_off >= dict_limits['min_cap_wind_off'] * bin_wind_off)
    
    """Maximum capacity constraints"""
    #Available biomass constraint
    if biomass_constrained:
        m.addConstr( gp.quicksum(f_gasif_wood[t]*delta_t for t in range(T)) <= dict_limits['max_bm_annu'] ) 

    #Set constraint on maximum hydrogen which can be exported
    m.addConstr( gp.quicksum(export_h2[t]*delta_t for t in range(T)) <= ((MJ_KG_H2/3.6) * dict_limits['max_h2_export']) ) # max_h2_export is in tonnes
    
    # Biogas
    for t in range(0, T, 24):
        m.addConstr((gp.quicksum(f_biogas[i] * delta_t for i in range(t, (t + 24)))) <= (dict_limits['biogas_potential_MWh'] / 365))
    
    for t in range(T):
        m.addConstr( f_chp_biogas[t] == f_biogas[t])# * biogas_eff)    

    """Other maximum capacity constraints"""
    if dict_limits['min_electrolyzer']>0: 
        #electrolyzer
        m.addConstr(cap_electrolyzer <= dict_limits['max_electrolyzer'] * bin_electrolyzer)
        m.addConstr(cap_electrolyzer >= dict_limits['min_electrolyzer'] * bin_electrolyzer)
    
    if dict_limits['min_cap']>0:
        # Maximum battery capacity installed
        m.addConstr(cap_bat_en <= dict_limits['max_bat'] * bin_bat_en)
        m.addConstr(cap_bat_en >= dict_limits['min_cap'] * bin_bat_en)

        # hydrogen
        m.addConstr(cap_h2_ves >= dict_limits['min_cap'] * bin_h2_ves)
        m.addConstr(cap_h2_ves <= dict_limits['max_h2_storage'] * bin_h2_ves)

        # fuel cell
        m.addConstr(cap_fc <= dict_limits['max_fc'] * bin_fc)
        m.addConstr(cap_fc >= dict_limits['min_cap'] * bin_fc)

        # Gasification 
        m.addConstr(cap_gasif <= dict_limits['max_gasif'] * bin_gasif)
        m.addConstr(cap_gasif >= dict_limits['min_cap'] * bin_gasif)
        
    if dict_limits['min_cap_ind']>0:
        # chp mix
        m.addConstr(cap_chp_mix <= dict_limits['max_chp_mix'] * bin_chp_mix)
        m.addConstr(cap_chp_mix >= dict_limits['min_cap_ind'] * bin_chp_mix)
        
        # gas boiler
        m.addConstr(cap_d_boiler <= dict_limits['max_boiler'] * bin_boiler)
        m.addConstr(cap_d_boiler >= dict_limits['min_cap_ind'] * bin_boiler)
        
        # electric boiler
        m.addConstr(cap_e_boiler <= dict_limits['max_boiler'] * bin_e_boiler)
        m.addConstr(cap_e_boiler >= dict_limits['min_cap_ind'] * bin_e_boiler)
        
        # Limit the amount of technologies used for high T heat, to reduce system complexity and integration issues.
        m.addConstr(bin_chp_mix + bin_boiler + bin_e_boiler <= 2)        
        
    # hp
    if dict_limits['min_cap_res']>0:
        m.addConstr(cap_hp <= dict_limits['max_hp'] * bin_hp)
        m.addConstr(cap_hp >= dict_limits['min_cap_res'] * bin_hp)

        # gas boiler - residential
        m.addConstr(cap_gb_res <= dict_limits['max_gb_res'] * bin_gb_res)
        m.addConstr(cap_gb_res >= dict_limits['min_cap_res'] * bin_gb_res)

        # Capacity storage medium
        if res_heat_stor:
            m.addConstr(cap_heat_storage <= dict_limits['max_heat_storage'] * bin_heat_storage)
            m.addConstr(cap_heat_storage >= dict_limits['min_cap_res'] * bin_heat_storage)

        # e heating
        m.addConstr(cap_e_heating <= bin_e_heating * dict_limits['st_max'])
        m.addConstr(cap_e_heating >= dict_limits['min_cap_res'] * bin_e_heating)

    # cap TS
    if dict_limits['min_solar_heat_ha']>0:
        m.addConstr(cap_solar_heat <= dict_limits['max_solar_heat_ha'] * bin_solar_heat)
        m.addConstr(cap_solar_heat >= dict_limits['min_solar_heat_ha'] * bin_solar_heat)

    """3.1. Local balance electricity"""
    for t in range(T):
        m.addConstr( (p_grid_abs[t] - p_grid_inj[t]) + 
                      (p_battdis[t] - p_battch[t]) + pv_total[t] + wind_total_off[t] 
                     + wind_total_on[t] + p_fc[t] + p_chp_mix[t]
                     == f_hp[t] + df_data.electricity_demand[t] + f_bev_1[t] + f_bev_2[t] + f_bev_3[t] + f_elect[t] + f_e_boiler[t] + f_e_heating[t])
    
    # Constraint(s) to ensure balanced energy autonomy
    if autonomous_balanced:
        # Locally generated renewable energy sources (wind, solar, and biomass) >= power demand, including (some) demand for heat provision
        m.addConstr( gp.quicksum( (pv_total[t] + wind_total_off[t] 
                     + wind_total_on[t] + p_chp_mix[t]) * delta_t for t in range(T)) >= gp.quicksum( (f_hp[t] + df_data.electricity_demand[t] + 
                                                                                                   f_bev_1[t] + f_bev_2[t] + f_bev_3[t] + f_elect[t] + 
                                                                                                   f_e_boiler[t] + f_e_heating[t]) * delta_t for t in range(T) ) ) 

    """ Balances energy carriers"""
    #For the analysis without high T heat
    for t in range(T):
        """ Heat balance - high T"""  
        if high_t_heat:
            m.addConstr(p_th_chp_mix[t] + p_d_boiler[t] + p_e_boiler[t] == df_data.ind_heat_demand[t])
        else:
            m.addConstr(p_th_chp_mix[t] + p_d_boiler[t] + p_e_boiler[t] == 0)  
            
        """ Residential heat balance - low T"""   
        # if residential heat storage is included, we have to consider charging and discharging of it
        if res_heat_stor:
            m.addConstr(p_gb_res[t] + p_solar_heat[t] + p_fc_heat[t] + p_th_hp[t] + p_e_heating[t] + (heat_stor_dis[t] - heat_stor_ch[t]) == df_data.bld_heat_demand[t])
        else:
            m.addConstr(p_gb_res[t] + p_solar_heat[t] + p_fc_heat[t] + p_th_hp[t] + p_e_heating[t] == df_data.bld_heat_demand[t])            
            
        """ Hydrogen balance"""
        m.addConstr(p_elect[t] == f_chp_h2[t] + f_fc[t] + export_h2[t] + p_h2_ves[t])

        """ Syngas balance"""
        m.addConstr(p_gasif[t] == f_chp_syngas[t])
        
    "Electric high T boiler - Industrial heat"
    for t in range(T): 
        m.addConstr(p_e_boiler[t] == f_e_boiler[t] * parm['e_boiler_eff'])
        m.addConstr(f_e_boiler[t] <= cap_e_boiler) 
        
    """3.2. Power boundaries"""
    for t in range(T):
        m.addConstr(p_grid_abs[t] <= dict_limits['max_grid_cap'] * bin_grid[t])
        # Use GUROBI indicator constraint instead
        m.addGenConstrIndicator(bin_grid[t], True, p_grid_inj[t], gp.GRB.EQUAL, 0)

    # Get the maximum abs or inj peak and store in var to be used for demand charge
    for t in range(T):
        m.addConstr(p_grid_inj[t] <= p_grid_connection)
        m.addConstr(p_grid_abs[t] <= p_grid_connection)

    """3.2. Wood gasification"""
    for t in range(T):
        m.addConstr(p_gasif[t] == (f_gasif_wood[t] * parm['gasif_eff'])) 
        m.addConstr(f_gasif_wood[t] <= cap_gasif)         
            
    if consider_down_times:
        wg_TU = int(parm['wg_TU'])
        wg_TD = int(parm['wg_TD'])
        
        # Add linearization constraints aux_gasif <= cap_gasif * x_wg[t]
        # to avoid m.addConstrs(f_gasif_wood[t] <= cap_gasif * x_wg[t] for t in range(T))
        for t in range(T):  
            m.addConstr(f_gasif_wood[t] <= aux_gasif[t]) 
            m.addConstr(f_gasif_wood[t] >= parm['gasif_pl_min'] * aux_gasif[t])
            
            # Add linearization constraints for f_gasif_wood[t] <= cap_gasif * x_wg[t]:
            m.addConstr(aux_gasif[t] <= dict_limits['max_gasif'] * x_wg[t])
            m.addConstr(dict_limits['min_cap'] * x_wg[t] <= aux_gasif[t])

            m.addConstr(cap_gasif - dict_limits['max_gasif'] * (1-x_wg[t]) <= aux_gasif[t])
            m.addConstr(aux_gasif[t] <= cap_gasif)

        for t in range(0,1):       
            # For timestep 0, we assume same or less output as t+1
            #2b
            m.addConstr(f_gasif_wood[0] <= f_gasif_wood[t+1])

        for t in range(1,T):       
            # Skip timestep 0, otherwise keyError issue
            #2b
            if parm['wg_ramp']>1:
                m.addConstr(f_gasif_wood[t] - f_gasif_wood[t-1] <= (cap_gasif/parm['wg_ramp']) )
                #2c
                m.addConstr(f_gasif_wood[t-1] - f_gasif_wood[t] <= (cap_gasif/parm['wg_ramp']) )

            # Logical constraint, 2d:
            m.addConstr(x_wg[t] - x_wg[t-1] == y_wg[t] - z_wg[t])

        if wg_TU>1:
            # 2e Minimum up- and downtimes, constraint 6 for uptimes:
            for t in range(wg_TU, T):   
                m.addConstr(gp.quicksum(y_wg[i] for i in range(t-wg_TU+1,t+1)
                                       ) <= x_wg[t])

        if wg_TD>1:
            # 2f Minimum up- and downtimes, constraint 6 for uptimes:
            for t in range(wg_TD, T):
                m.addConstr(gp.quicksum(z_wg[i] for i in range(t-wg_TD+1,t+1)
                                    ) <= 1 - x_wg[t])   
    """3.3. Battery model"""
    # Battery dynamics
    for t in range(0,1):
        m.addConstr(E_batt[t] == E_batt[0] * (1-parm['bat_dis_loss']*delta_t) + 
                    (parm['bat_eff_ch']*p_battch[t]*delta_t) - ((p_battdis[t]*delta_t)/parm['bat_eff_dis']))

    for t in range(1,T):
        m.addConstr(E_batt[t] == E_batt[t-1] * (1-parm['bat_dis_loss']*delta_t) + 
                    (parm['bat_eff_ch']*p_battch[t]*delta_t) - ((p_battdis[t]*delta_t)/parm['bat_eff_dis']))

    # Periodicity constraint
    m.addConstr(E_batt[T-1] == E_batt[0]) 
    
    for t in range(T):
        # Use generator constraint instead, to improve performance https://www.gurobi.com/documentation/9.5/refman/py_model_agc_xxx.html
        m.addConstr(p_battch[t] <= dict_limits['max_bat'] * bin_bat[t] )
        m.addGenConstrIndicator(bin_bat[t], True, p_battdis[t], gp.GRB.EQUAL, 0)
        
        m.addConstr(p_battch[t] <= cap_bat_p )
        m.addConstr(p_battdis[t] <= cap_bat_p )
        
        # respect SoC
        m.addConstr(E_batt[t] >= cap_bat_en * parm['bat_soc_min'])
        m.addConstr(E_batt[t] <= cap_bat_en * parm['bat_soc_max'])

    """3.4. Hydrogen storage, electrolyzer and fuel cell"""
    for t in range(T):
        m.addConstr(p_elect[t] == f_elect[t] * parm['electr_eff'])  
        m.addConstr(f_elect[t] <= cap_electrolyzer)  
    
    # Periodicity constraint
    m.addConstr(E_h2_ves[T-1] == E_h2_ves[0]) 
    
    # Within this balance
    for t in range(0,1):
        m.addConstr(E_h2_ves[t] == E_h2_ves[0] + p_h2_ves[t])

    for t in range(1,T):
        m.addConstr(E_h2_ves[t] == E_h2_ves[t-1] + p_h2_ves[t])
        
    for t in range(T):
        m.addConstr(E_h2_ves[t] <= cap_h2_ves) 

    # consider ramp rates h2_stor_ramp
    for t in range(T):
        m.addConstr(p_h2_ves[t] <= (cap_h2_ves/parm['h2_ramp_ves']) )
        m.addConstr( -(cap_h2_ves/parm['h2_ramp_ves']) <= p_h2_ves[t] )
        
    # FCCHP
    for t in range(T):
        # FCCHP
        m.addConstr(p_fc[t] == f_fc[t] * parm['fc_eff_e'])
        m.addConstr(p_fc_heat[t] == f_fc[t] * parm['fc_eff_th'])
        m.addConstr(f_fc[t] <= cap_fc)

    """Heat storage dynamics"""
    # Steen et al.: https://www.sciencedirect.com/science/article/pii/S0306261914007181
    # Using the approach in Steen et al. and final approach from Gabrielli et al.:
    if res_heat_stor:
        for t in range(0,1):    
            m.addConstr(E_heat_stor[t] == E_heat_stor[0] * (1-parm['heat_stor_dis_loss']*delta_t) - parm['heat_stor_loss_coeff'] * cap_heat_storage * (
                (parm['heat_stor_t_min']-df_data.temperature[t])/(parm['heat_stor_t_max']-parm['heat_stor_t_min'])) + 
                        (heat_stor_ch[t] * delta_t * parm['heat_stor_eff_ch']) - ((heat_stor_dis[t]*delta_t)/parm['heat_stor_eff_dis']))          

        for t in range(1,T):
            m.addConstr(E_heat_stor[t] == E_heat_stor[t-1] * (1-parm['heat_stor_dis_loss']*delta_t) - parm['heat_stor_loss_coeff'] * cap_heat_storage * (
                (parm['heat_stor_t_min']-df_data.temperature[t])/(parm['heat_stor_t_max']-parm['heat_stor_t_min'])) + 
                        (heat_stor_ch[t] * delta_t * parm['heat_stor_eff_ch']) - ((heat_stor_dis[t]*delta_t)/parm['heat_stor_eff_dis'])) 

        # Periodicity constraint
        m.addConstr(E_heat_stor[T-1] == E_heat_stor[0]) 

        for t in range(T):  
            # For the heat storage capacity determination size 
            m.addConstr(E_heat_stor[t] <= cap_heat_storage) 

            m.addConstr(heat_stor_ch[t] <= (cap_heat_storage/parm['heat_stor_ramp']) )
            m.addConstr(heat_stor_dis[t] <= (cap_heat_storage/parm['heat_stor_ramp']) )
            m.addConstr( -(cap_heat_storage/parm['heat_stor_ramp']) <= heat_stor_ch[t] )
            m.addConstr( -(cap_heat_storage/parm['heat_stor_ramp']) <= heat_stor_dis[t] ) 

            # Use generator constraint instead, to improve performance https://www.gurobi.com/documentation/9.5/refman/py_model_agc_xxx.html
            m.addConstr(heat_stor_ch[t] <= dict_limits['max_heat_storage'] * bin_h_stor[t] )
            m.addGenConstrIndicator(bin_h_stor[t], True, heat_stor_dis[t], gp.GRB.EQUAL, 0)

    """BEV"""
    # Charging schedules
    schedules = ['TP_1', "TP_2", "TP_3"]
    char_schedules = len(schedules) #assumed to be evenly distributed
    
    # Define common charging constraints
    common_charge_constraint = (share_bev * parm['km_per_day'] * parm['bev_MWh_km'] * total_hh_r / char_schedules)
    
    for t in range(0, T, 24):
        m.addConstr((gp.quicksum(f_bev_1[i] for i in range(t, (t+24)))) == common_charge_constraint)
        m.addConstr((gp.quicksum(f_bev_2[i] for i in range(t, (t+24)))) == common_charge_constraint)
        m.addConstr((gp.quicksum(f_bev_3[i] for i in range(t, (t+24)))) == common_charge_constraint)

    for t in range(T):
        if parm['bev_min_p_ch']>0:
            m.addConstr(f_bev_1[t] >= ( df_data.TP_1[t] * parm['bev_min_p_ch'] * total_hh_r * share_bev * (1/char_schedules) ))  # Charging schedule 1
            m.addConstr(f_bev_2[t] >= ( df_data.TP_2[t] * parm['bev_min_p_ch'] * total_hh_r * share_bev * (1/char_schedules) ))  # Charging schedule 2
            m.addConstr(f_bev_3[t] >= ( df_data.TP_3[t] * parm['bev_min_p_ch'] * total_hh_r * share_bev * (1/char_schedules) ))  # Charging schedule 3
            
        m.addConstr(f_bev_1[t] <= ( df_data.TP_1[t] * parm['bev_max_p_ch'] * total_hh_r * share_bev * (1/char_schedules) ))  # Charging schedule 1
        m.addConstr(f_bev_2[t] <= ( df_data.TP_2[t] * parm['bev_max_p_ch'] * total_hh_r * share_bev * (1/char_schedules) ))  # Charging schedule 2
        m.addConstr(f_bev_3[t] <= ( df_data.TP_3[t] * parm['bev_max_p_ch'] * total_hh_r * share_bev * (1/char_schedules) ))  # Charging schedule 3

    """CHP: mixed fuel"""
    for t in range(T):   
        m.addConstr(f_chp_mix[t] == (f_chp_h2[t] + f_chp_syngas[t] + f_chp_biogas[t]) )
        m.addConstr(p_chp_mix[t] == f_chp_mix[t] * parm['chp_mix_eff_e'])
        m.addConstr(p_th_chp_mix[t] == f_chp_mix[t] * parm['chp_mix_eff_th'])
        m.addConstr(f_chp_h2[t] <= parm['chp_mix_max_h2_share'] * f_chp_mix[t] ) 
        m.addConstr(f_chp_mix[t] <= cap_chp_mix) 
    
    if consider_down_times:
        chp_mix_TU = int(parm['chp_mix_TU'])
        chp_mix_TD = int(parm['chp_mix_TD'])
        
        # Add linearization constraints aux_chp_mix <= cap_chp_mix * x_chp_mix[t]
        # to avoid m.addConstrs(f_chp_mix[t] <= cap_chp_mix * x_chp_mix[t] for t in range(T))
        for t in range(T):  
            # Consider part load ratio
            m.addConstr(f_chp_mix[t] <= aux_chp_mix[t]) 
            m.addConstr(f_chp_mix[t] >= parm['chp_mix_pl_min'] * aux_chp_mix[t]) 

            # Add linearization constraints for f_chp_mix[t] <= cap_chp_mix * x_chp_mix[t]:
            m.addConstr(aux_chp_mix[t] <= dict_limits['max_chp_mix'] * x_chp_mix[t])
            m.addConstr(dict_limits['min_chp_mix'] * x_chp_mix[t] <= aux_chp_mix[t])

            m.addConstr(cap_chp_mix - dict_limits['max_chp_mix'] * (1-x_chp_mix[t]) <= aux_chp_mix[t])
            m.addConstr(aux_chp_mix[t] <= cap_chp_mix)

        for t in range(0,1):       
            # For timestep 0, we assume same or less output as t+1
            #2b
            m.addConstr(f_chp_mix[0] <= f_chp_mix[t+1])

        for t in range(1,T):       
            # Skip timestep 0, otherwise keyError issue
            #2b
            if parm['chp_mix_ramp']>1:
                m.addConstr(f_chp_mix[t] - f_chp_mix[t-1] <= (cap_chp_mix/parm['chp_mix_ramp']) )
                #2c
                m.addConstr(f_chp_mix[t-1] - f_chp_mix[t] <= (cap_chp_mix/parm['chp_mix_ramp']) )

            # Logical constraint, 2d:
            m.addConstr(x_chp_mix[t] - x_chp_mix[t-1] == y_chp_mix[t] - z_chp_mix[t])

        if chp_mix_TU>1:
            # 2e Minimum up- and downtimes, constraint 6 for uptimes:
            for t in range(chp_mix_TU, T):   
                m.addConstr(gp.quicksum(y_chp_mix[i] for i in range(t-chp_mix_TU+1,t+1)
                                       ) <= x_chp_mix[t])

        if chp_mix_TD>1:
            # 2f Minimum up- and downtimes, constraint 6 for uptimes:
            for t in range(chp_mix_TD, T):
                m.addConstr(gp.quicksum(z_chp_mix[i] for i in range(t-chp_mix_TD+1,t+1)
                                    ) <= 1 - x_chp_mix[t])   

    "Diesel boiler - Industrial heat"
    for t in range(T): 
        m.addConstr(p_d_boiler[t] == f_d_boiler[t] * parm['d_boiler_eff'])
        m.addConstr(f_d_boiler[t] <= cap_d_boiler) 

    "NG boiler - Residential heat"
    for t in range(T): 
        m.addConstr(p_gb_res[t] == f_gb_res[t] * parm['gb_eff_res'])
        m.addConstr(f_gb_res[t] <= cap_gb_res) 
        
    "Electric - Residential heat"
    for t in range(T): 
        m.addConstr(p_e_heating[t] == f_e_heating[t] * parm['e_heating_eff'])
        m.addConstr(f_e_heating[t] <= cap_e_heating) 

    "Thermal solar (TS)"
    for t in range(T): 
        m.addConstr(p_solar_heat[t] <= (10000/1000) * df_data.irradiance[t] * parm['solar_th_eff'] * cap_solar_heat) # Conversion factor included from ha and MW
        
    "HP"
    for t in range(T): 
        m.addConstr(p_th_hp[t] == f_hp[t] * df_data.cop_array[t]) 
        m.addConstr(f_hp[t] <= cap_hp) # take the power capacity of the heat pump.

    """3.2. Wind and PV supply as well as curtailment"""
    for t in range(T): 
        m.addConstr(wind_total_off[t] <= cap_wind_off * df_data.wind_MW_array_off[t])
        m.addConstr(wind_total_on[t] <= cap_wind_on * wind_MW_array_on[t])
        m.addConstr(pv_total[t] <= cap_pv * df_data.pv_MW_array[t])

    """
    Step 4: Set objective one of multi-objective function
    """
    petrol_car_kms = share_gasoline * DAYS * parm['km_per_day'] * total_hh_r
    bev_car_kms = share_bev * DAYS * parm['km_per_day'] * total_hh_r

    """Objective 1: Costs"""
    """Operation (variable OPEX)"""
    an_op_wood = gp.quicksum(f_gasif_wood[t] * parm['wood_price'] for t in range(T))
    an_op_d_boiler = gp.quicksum(f_d_boiler[t] * parm['diesel_price'] for t in range(T))
    an_op_ng_res = gp.quicksum(f_gb_res[t] * parm['ng_price'] for t in range(T))
    an_op_grid_abs = gp.quicksum(p_grid_abs[t] * df_data.grid_abs_price[t] for t in range(T))
    an_op_grid_inj = gp.quicksum(p_grid_inj[t] * df_data.rev_inj[t] for t in range(T))
    an_op_export_h2 = gp.quicksum(export_h2[t] * parm['h2_sell_price'] for t in range(T))
    an_op_biogas = gp.quicksum(f_biogas[t] * parm['biogas_price'] for t in range(T))
    
    # Add fuel costs for gasoline cars, if not based on BEVs only
    an_op_gasoline_cars = petrol_car_kms * parm['price_km_gasoline']

    an_op = (an_op_wood +
                        an_op_d_boiler + an_op_ng_res + an_op_biogas +
                        an_op_grid_abs + an_op_gasoline_cars - an_op_grid_inj -
                        an_op_export_h2
                       )

    """Investments (CAPEX)"""
    an_capex_boiler = calc_crf(parm['dr'],parm['project_lt']) * (cap_d_boiler * parm['d_boiler_capex'])   
    an_capex_gb_res = calc_crf(parm['dr'],parm['project_lt']) * (cap_gb_res * parm['gb_capex_res']) 
    an_capex_e_heating = calc_crf(parm['dr'],parm['project_lt']) * (cap_e_heating * parm['e_heating_capex']) 
    
    an_capex_hp = calc_crf(parm['dr'],parm['project_lt']) * (cap_hp * parm['hp_capex'])
    an_capex_electrolyzer = calc_crf(parm['dr'],parm['project_lt']) * (cap_electrolyzer * parm['electr_capex'])
    an_capex_fc = calc_crf(parm['dr'],parm['project_lt']) * (cap_fc * parm['fc_capex'])   
    an_capex_h2_ves = calc_crf(parm['dr'],parm['project_lt']) * (cap_h2_ves * parm['h2_ves_capex']) 
    an_capex_pv = calc_crf(parm['dr'],parm['project_lt']) * (cap_pv * parm['pv_capex']) 
    an_capex_wind_off = calc_crf(parm['dr'],parm['project_lt']) *(cap_wind_off * parm['wind_off_capex'])
    an_capex_wind_on = calc_crf(parm['dr'],parm['project_lt']) * (cap_wind_on * parm['wind_on_capex'])
    an_capex_bat_en = calc_crf(parm['dr'],parm['project_lt']) * (cap_bat_en * parm['bat_en_capex']) 
    an_capex_bat_p = calc_crf(parm['dr'],parm['project_lt']) * (cap_bat_p * parm['bat_p_capex'])
    an_capex_heat_storage = calc_crf(parm['dr'],parm['project_lt']) * (cap_heat_storage * parm['heat_stor_capex'])
    an_capex_p_solar_heat = calc_crf(parm['dr'],parm['project_lt']) * (cap_solar_heat * parm['solar_th_capex'])
    an_capex_chp_mixer = calc_crf(parm['dr'],parm['project_lt']) * (cap_chp_mix * parm['chp_mix_capex'] * parm['chp_mix_eff_e']) # The costs for the CHP unit is provided per MWe power delivered
    an_capex_gasif = calc_crf(parm['dr'],parm['project_lt']) * (cap_gasif * parm['gasif_capex'])
    an_capex_grid_ins = calc_crf(parm['dr'],parm['project_lt']) * (p_grid_connection * parm['grid_capex'])
    an_capex_e_boiler = calc_crf(parm['dr'],parm['project_lt']) * (cap_e_boiler * parm['e_boiler_capex'])   
    
    an_capex_trans_gas = share_gasoline * total_hh_r * calc_crf(parm['dr'],parm['project_lt']) * parm['gas_car_capex']
    an_capex_trans_bev = share_bev * total_hh_r * calc_crf(parm['dr'],parm['project_lt']) * parm['bev_car_capex']

    an_capex  = (an_capex_boiler + an_capex_gb_res + an_capex_e_heating + an_capex_hp + an_capex_electrolyzer +
                an_capex_fc + an_capex_h2_ves + 
                an_capex_pv + an_capex_wind_off + an_capex_wind_on +
                an_capex_bat_en + an_capex_bat_p + an_capex_heat_storage +
                an_capex_p_solar_heat + an_capex_chp_mixer +
                an_capex_gasif + an_capex_grid_ins + an_capex_e_boiler +
                an_capex_trans_gas + an_capex_trans_bev
               )

    """Replacements"""
    an_rep = (
          rep_annual_int((cap_d_boiler * parm['d_boiler_capex']),parm['project_lt'],parm['dr'], parm['d_boiler_lt'])   
        + rep_annual_int((cap_gb_res * parm['gb_capex_res']),parm['project_lt'],parm['dr'], parm['gb_lt_res'])  
        + rep_annual_int((cap_e_heating * parm['e_heating_capex']),parm['project_lt'],parm['dr'], parm['e_heating_lt'])  
        + rep_annual_int((cap_hp * parm['hp_capex']),parm['project_lt'],parm['dr'],parm['hp_lt'])
        + rep_annual_int((cap_h2_ves * parm['h2_ves_capex']),parm['project_lt'],parm['dr'],parm['h2_ves_lt'])
        + rep_annual_int((cap_electrolyzer * parm['electr_capex']),parm['project_lt'],parm['dr'],parm['electr_lt']) 
        + rep_annual_int((cap_fc * parm['fc_capex']),parm['project_lt'],parm['dr'],parm['fc_lt'])
        + rep_annual_int((cap_pv * parm['pv_capex']),parm['project_lt'],parm['dr'],parm['pv_lt'])
        + rep_annual_int((cap_wind_off * parm['wind_off_capex']),parm['project_lt'],parm['dr'],parm['wind_off_lt'])
        + rep_annual_int((cap_wind_on * parm['wind_on_capex']),parm['project_lt'],parm['dr'],parm['wind_on_lt'])
        + rep_annual_int((cap_bat_en * parm['bat_en_capex']),parm['project_lt'],parm['dr'],parm['bat_en_lt'])
        + rep_annual_int((cap_bat_p * parm['bat_p_capex']),parm['project_lt'],parm['dr'],parm['bat_p_lt'])
        + rep_annual_int((cap_heat_storage * parm['heat_stor_capex']),parm['project_lt'],parm['dr'],parm['heat_stor_lt'])
        + rep_annual_int((cap_solar_heat * parm['solar_th_capex']),parm['project_lt'],parm['dr'],parm['solar_th_lt'])
        + rep_annual_int((cap_chp_mix * parm['chp_mix_capex'] * parm['chp_mix_eff_e']),parm['project_lt'],parm['dr'],parm['chp_mix_lt'])# The costs for the CHP unit is provided per MWe power delivered
        + rep_annual_int((cap_gasif * parm['gasif_capex']),parm['project_lt'],parm['dr'],parm['gasif_lt'])
        + rep_annual_int((p_grid_connection * parm['grid_capex']),parm['project_lt'],parm['dr'],parm['grid_lt'])
        + rep_annual_int((cap_e_boiler * parm['e_boiler_capex']),parm['project_lt'],parm['dr'],parm['e_boiler_lt'])
        
        + rep_annual_int((share_gasoline * total_hh_r * parm['gas_car_capex']),parm['project_lt'],parm['dr'],parm['gas_car_lt'])
        + rep_annual_int((share_bev * total_hh_r * parm['bev_car_capex']),parm['project_lt'],parm['dr'],parm['bev_car_lt'])        
             )

    """Fixed O&M (without replacements)"""
    an_om = ((cap_d_boiler * parm['d_boiler_capex'] * parm['d_boiler_om']) 
            + (cap_gb_res * parm['gb_capex_res'] * parm['gb_om_res']) 
            + (cap_e_heating * parm['e_heating_capex'] * parm['e_heating_om']) 
            + (cap_hp * parm['hp_capex'] * parm['hp_om']) 
            + (cap_h2_ves * parm['h2_ves_capex'] * parm['h2_ves_om'])          
            + (cap_electrolyzer * parm['electr_capex'] * parm['electr_om'])
            + (cap_fc * parm['fc_capex'] * parm['fc_om']) 
            + (cap_pv * parm['pv_capex'] * parm['pv_om']) 
            + (cap_wind_off * parm['wind_off_capex'] * parm['wind_off_om']) 
            + (cap_wind_on * parm['wind_on_capex'] * parm['wind_on_om']) 
            + (cap_bat_en * parm['bat_en_capex'] * parm['bat_en_om']) 
            + (cap_bat_p * parm['bat_p_capex'] * parm['bat_p_om'])
            + (cap_heat_storage * parm['heat_stor_capex'] * parm['heat_stor_om']) 
            + (cap_solar_heat * parm['solar_th_capex'] * parm['solar_th_om']) 
            + (cap_chp_mix * parm['chp_mix_capex'] * parm['chp_mix_om'] * parm['chp_mix_eff_e']) # The costs for the CHP unit is provided per MWe power delivered
            + (cap_gasif * parm['gasif_capex'] * parm['gasif_om']) 
            + (p_grid_connection * parm['grid_capex'] * parm['grid_om']) 
            + (cap_e_boiler * parm['e_boiler_capex'] * parm['e_boiler_om']) 
             
            + (share_gasoline * parm['gas_car_capex'] * total_hh_r * parm['gas_car_om']) 
            + (share_bev * parm['bev_car_capex'] * total_hh_r * parm['bev_car_om']) 
            )
    

    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    """Objective 2: life cycle GHG emissions"""
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    """Operation"""
    # Note: biogas is assumed to be burden-free
    an_ghg_op_wood = gp.quicksum(p_gasif[t] * dict_ghg['ghg_imp_syngas'] for t in range(T))
    an_ghg_op_d_boiler = gp.quicksum(p_d_boiler[t] * dict_ghg['ghg_imp_diesel_ht'] for t in range(T)) 
    an_ghg_op_ng_res = gp.quicksum(p_gb_res[t] * dict_ghg['ghg_imp_ng'] for t in range(T))
    an_ghg_op_h2 = gp.quicksum(export_h2[t] * dict_ghg['ghg_imp_smr'] for t in range(T))
    an_ghg_op_chp = gp.quicksum(f_chp_mix[t] * dict_ghg['ghg_imp_chp_fuel'] for t in range(T))
    an_ghg_op_grid_abs = gp.quicksum( (p_grid_abs[t] * df_data.ghg_impact[t]) for t in range(T)) 
    an_ghg_op_grid_inj = gp.quicksum( (p_grid_inj[t] * df_data.ghg_impact_cons[t]) for t in range(T))
    
    an_ghg_op_gasoline_cars = petrol_car_kms * dict_ghg['ghg_gas_car_km']

    # No purchase of heat as it will be balanced in the optimal design
    an_op_ghg = (an_ghg_op_wood + an_ghg_op_d_boiler + an_ghg_op_chp + an_ghg_op_grid_abs + 
                 an_ghg_op_ng_res + an_ghg_op_gasoline_cars - an_ghg_op_grid_inj - an_ghg_op_h2)

    """Production and replacements"""    
    an_ghg_boiler = (cap_d_boiler * dict_ghg['ghg_imp_gb'] * (parm['project_lt']/parm['d_boiler_lt'])) / parm['project_lt']  
    an_ghg_gb_res = (cap_gb_res * dict_ghg['ghg_imp_gb'] * (parm['project_lt']/parm['gb_lt_res'])) / parm['project_lt']
    an_ghg_e_heating = (cap_e_heating * dict_ghg['ghg_imp_e_heating'] * (parm['project_lt']/parm['e_heating_lt'])) / parm['project_lt']  
    an_ghg_hp = (cap_hp * dict_ghg['ghg_imp_hp'] * (parm['project_lt']/parm['hp_lt'])) / parm['project_lt']  
    an_ghg_electrolyzer = (cap_electrolyzer * dict_ghg['ghg_imp_electr'] * (parm['project_lt']/parm['electr_lt'])) / parm['project_lt']  
    an_ghg_fc = parm['fc_eff_e'] * (cap_fc * dict_ghg['ghg_imp_fc'] * (parm['project_lt']/parm['fc_lt'])) / parm['project_lt']     
    an_ghg_h2_ves = (cap_h2_ves * dict_ghg['ghg_imp_h2_ves'] * (parm['project_lt']/parm['h2_ves_lt'])) / parm['project_lt']     
    an_ghg_pv = (cap_pv * dict_ghg['ghg_imp_pv'] * (parm['project_lt']/parm['pv_lt'])) / parm['project_lt']  
    an_ghg_wind_off = (cap_wind_off * dict_ghg['ghg_imp_wind_off'] * (parm['project_lt']/parm['wind_off_lt'])) / parm['project_lt']  
    an_ghg_wind_on = (cap_wind_on * dict_ghg['ghg_imp_wind_on'] * (parm['project_lt']/parm['wind_on_lt'])) / parm['project_lt']   
    an_ghg_bat_en = (cap_bat_en * dict_ghg['ghg_imp_bat_cap'] * (parm['project_lt']/parm['bat_en_lt'])) / parm['project_lt']  
    an_ghg_bat_p = (cap_bat_p * dict_ghg['ghg_imp_bat_p'] * (parm['project_lt']/parm['bat_p_lt'])) / parm['project_lt']  
    an_ghg_heat_storage = (cap_heat_storage * dict_ghg['ghg_imp_heat_storage'] * (parm['project_lt']/parm['heat_stor_lt'])) / parm['project_lt']  
    an_ghg_p_solar_heat = (cap_solar_heat * dict_ghg['ghg_imp_solar_heat_ha'] * (parm['project_lt']/parm['solar_th_lt'])) / parm['project_lt']  
    an_ghg_chp_mixer = (cap_chp_mix * dict_ghg['ghg_imp_chp_mix'] * (parm['project_lt']/parm['chp_mix_lt'])) / parm['project_lt']   # The costs for the CHP unit is provided per MWe power delivered
    an_ghg_gasif = (cap_gasif * dict_ghg['ghg_imp_gasif'] * (parm['project_lt']/parm['gasif_lt'])) / parm['project_lt']  
    an_ghg_grid_ins = (p_grid_connection * dict_ghg['ghg_impact_grid_network'] * (parm['project_lt']/parm['grid_lt'])) / parm['project_lt']
    an_ghg_e_boiler = (cap_e_boiler * dict_ghg['ghg_imp_e_heating'] * (parm['project_lt']/parm['e_boiler_lt'])) / parm['project_lt']  
    
    an_ghg_bev_cars = bev_car_kms * dict_ghg['ghg_bev_car_km'] 

    an_op_ghg_inv = (an_ghg_boiler + 
                    an_ghg_gb_res + an_ghg_e_heating + an_ghg_hp + an_ghg_electrolyzer +
                    an_ghg_fc + an_ghg_h2_ves +
                    an_ghg_pv +  an_ghg_wind_off + an_ghg_wind_on +
                    an_ghg_bat_en + an_ghg_bat_p + an_ghg_heat_storage +
                    an_ghg_p_solar_heat + an_ghg_chp_mixer + 
                    an_ghg_gasif + an_ghg_grid_ins + an_ghg_e_boiler + an_ghg_bev_cars
               )

    obj2 = an_op_ghg + an_op_ghg_inv
  

    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    """Total costs, including possible costs for CO2"""
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    an_op_co2 = obj2 * (euro_ton_co2/1000) if euro_ton_co2!=0 else 0  
    obj1 = an_op + an_capex + an_rep + an_om + an_op_co2
        
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    """Set objectives"""
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    # add constraint to have x % less GHG emissions, implemented as inequality constraints to have a faster acceptable solution
    if eps_ghg_constraint:
        m.addConstr(obj2 <= eps_ghg_constraint)
        m.setObjective(obj1)       
    elif w_cost==1 and w_env==0:
        m.setObjective(obj1) 
    elif w_cost==0 and w_env==1:
        m.setObjective(obj2)         
    else:
        # Multi-objective problem with weights
        m.setObjectiveN(obj1, index = 0, weight = w_cost)
        m.setObjectiveN(obj2, index = 1, weight = w_env)

    """
    Step 5: Solve model
    """

    # Add error handling for optimization
    try:
        m.optimize()

        # Check optimization status
        if m.status == gp.GRB.OPTIMAL:
            print("Optimization successful. Objective value:", m.objVal)
            if save_for_warm_start:
                initial_solution = m.getAttr('X', m.getVars())
                # Export the solution to a file
                m.write('results\initial_solution_{}.sol'.format(int(obj1.getValue())))
        elif m.status == gp.GRB.INFEASIBLE:
            print("Optimization failed: The model is infeasible.")
        elif m.status == gp.GRB.UNBOUNDED:
            print("Optimization failed: The model is unbounded.")
        else:
            print("Optimization failed with status:", m.status)

    except gp.GurobiError as e:
        print("Gurobi Error:", e)

    """
    Step 6: Store variables values from optimal solution
    """ 

    """Get and store results in df"""     
    df_out = pd.DataFrame({
                            "time": df_data.index,
                            "electricity_demand":df_data.electricity_demand,
                            "ind_heat_demand":df_data.ind_heat_demand,
                            "bld_heat_demand":df_data.bld_heat_demand,
                            "p_grid_inj": m.getAttr('x',p_grid_inj).values(),
                            "p_grid_abs": m.getAttr('x',p_grid_abs).values(),
                            "bin_grid": m.getAttr('x',bin_grid).values(),
                            "p_battdis": m.getAttr('x',p_battdis).values(),
                            "p_battch": m.getAttr('x',p_battch).values(),
                            "E_batt": m.getAttr('x',E_batt).values(),
                            "bin_bat": m.getAttr('x',bin_bat).values(),
                            "f_bev_1": m.getAttr('x',f_bev_1).values(),
                            "f_bev_2": m.getAttr('x',f_bev_2).values(),
                            "f_bev_3": m.getAttr('x',f_bev_3).values(),
                            "E_h2_ves": m.getAttr('x',E_h2_ves).values(),
                            "p_h2_ves": m.getAttr('x',p_h2_ves).values(),
                            "f_elect": m.getAttr('x',f_elect).values(),
                            "p_elect": m.getAttr('x',p_elect).values(),
                            "fcchp_out": m.getAttr('x',p_fc).values(),
                            "p_fc_heat": m.getAttr('x',p_fc_heat).values(),
                            "f_fc": m.getAttr('x',f_fc).values(),
                            "pv_total": m.getAttr('x',pv_total).values(),
                            "wind_total_off": m.getAttr('x',wind_total_off).values(),
                            "wind_total_on": m.getAttr('x',wind_total_on).values(),
                            "p_th_hp": m.getAttr('x',p_th_hp).values(),
                            "f_hp": m.getAttr('x',f_hp).values(),
                            "p_d_boiler": m.getAttr('x',p_d_boiler).values(),
                            "f_d_boiler": m.getAttr('x',f_d_boiler).values(),
                            "E_heat_stor": m.getAttr('x',E_heat_stor).values() if res_heat_stor else [0] * len(df_data),
                            "heat_stor_ch": m.getAttr('x',heat_stor_ch).values() if res_heat_stor else [0] * len(df_data),
                            "heat_stor_dis": m.getAttr('x',heat_stor_dis).values() if res_heat_stor else [0] * len(df_data),
                            "bin_h_stor": m.getAttr('x',bin_h_stor).values() if res_heat_stor else [0] * len(df_data), 
                            "p_solar_heat": m.getAttr('x',p_solar_heat).values(),
                            "f_chp_h2": m.getAttr('x',f_chp_h2).values(),
                            "p_th_chp_mix": m.getAttr('x',p_th_chp_mix).values(),
                            "p_chp_mix": m.getAttr('x',p_chp_mix).values(),
                            "x_chp_mix": m.getAttr('x',x_chp_mix).values() if consider_down_times else [0] * len(df_data),
                            "y_chp_mix": m.getAttr('x',y_chp_mix).values() if consider_down_times else [0] * len(df_data),
                            "z_chp_mix": m.getAttr('x',z_chp_mix).values() if consider_down_times else [0] * len(df_data),
                            "aux_chp_mix": m.getAttr('x',aux_chp_mix).values() if consider_down_times else [0] * len(df_data),  
                            "f_chp_mix": m.getAttr('x',f_chp_mix).values(),
                            "f_chp_syngas": m.getAttr('x',f_chp_syngas).values(),
                            "p_gasif": m.getAttr('x',p_gasif).values(),
                            "f_gasif_wood": m.getAttr('x',f_gasif_wood).values(),
                            "x_wg": m.getAttr('x',x_wg).values() if consider_down_times else [0] * len(df_data),
                            "y_wg": m.getAttr('x',y_wg).values() if consider_down_times else [0] * len(df_data),
                            "z_wg": m.getAttr('x',z_wg).values() if consider_down_times else [0] * len(df_data),
                            "aux_gasif": m.getAttr('x',aux_gasif).values() if consider_down_times else [0] * len(df_data),
                            "export_h2": m.getAttr('x',export_h2).values(),
                            "p_gb_res": m.getAttr('x',p_gb_res).values(),
                            "f_gb_res": m.getAttr('x',f_gb_res).values(),
                            "f_biogas": m.getAttr('x',f_biogas).values(),
                            "f_chp_biogas": m.getAttr('x',f_chp_biogas).values(),
                            "p_e_heating": m.getAttr('x',p_e_heating).values(),
                            "f_e_heating": m.getAttr('x',f_e_heating).values(),  
                            "p_e_boiler": m.getAttr('x',p_e_boiler).values(),
                            "f_e_boiler": m.getAttr('x',f_e_boiler).values(),   
                                                       }).set_index("time")

    # Get single variables values
    cap_wind_off = cap_wind_off.x
    cap_wind_on = cap_wind_on.x
    cap_pv = cap_pv.x
    cap_bat_en = cap_bat_en.x
    cap_bat_p = cap_bat_p.x
    cap_h2_ves = cap_h2_ves.x
    cap_electrolyzer = cap_electrolyzer.x
    cap_fc = cap_fc.x
    cap_hp = cap_hp.x
    cap_d_boiler = cap_d_boiler.x
    cap_heat_storage = cap_heat_storage.x
    cap_solar_heat = cap_solar_heat.x
    cap_chp_mix = cap_chp_mix.x
    cap_gasif = cap_gasif.x
    cap_gb_res = cap_gb_res.x
    cap_e_heating = cap_e_heating.x
    cap_e_boiler = cap_e_boiler.x
    bin_wind_off = bin_wind_off.x
    bin_wind_on = bin_wind_on.x
    bin_pv = bin_pv.x
    bin_bat_en = bin_bat_en.x
    bin_h2_ves = bin_h2_ves.x
    bin_electrolyzer = bin_electrolyzer.x
    bin_fc = bin_fc.x
    bin_hp = bin_hp.x
    bin_boiler = bin_boiler.x
    bin_heat_storage = bin_heat_storage.x
    bin_solar_heat = bin_solar_heat.x
    bin_chp_mix = bin_chp_mix.x
    bin_gasif = bin_gasif.x
    bin_gb_res = bin_gb_res.x
    bin_e_boiler = bin_e_boiler.x
    bin_e_heating = bin_e_heating.x
    p_grid_connection = p_grid_connection.x
    share_bev = share_bev.x
    share_gasoline = share_gasoline.x

    # Check whether obj function is same as expectations in formula
    obj_value = m.objVal
    end = time.time()
    
    # Get calculation time
    total_time = (end-start)/3600 #hours
    
    # Calculate percentage change with initial results
    percent_c_costs = round((100*(obj1.getValue()-cost_init)/cost_init),2)
    percent_c_ghgs = round((100*(obj2.getValue()-ghg_init)/ghg_init),2)

    # Curtailment of renewables
    if ( sum(cap_wind_off * df_data.wind_MW_array_off) - sum(df_out.wind_total_off) ) > 0:
        curtailed_wind_off = (cap_wind_off * df_data.wind_MW_array_off) - df_out.wind_total_off
        ratio_wind_curtailed_off = curtailed_wind_off.sum() / sum(cap_wind_off * df_data.wind_MW_array_off)
    else:
        curtailed_wind_off = [0] * len(df_data)
        ratio_wind_curtailed_off = 0
        
    if ( sum(cap_wind_on *  wind_MW_array_on) - sum(df_out.wind_total_on) ) > 0:
        curtailed_wind_on = (cap_wind_on * wind_MW_array_on) - df_out.wind_total_on
        ratio_wind_curtailed_on = curtailed_wind_on.sum() / sum(cap_wind_on * wind_MW_array_on)
    else:
        curtailed_wind_on =[0] * len(df_data)
        ratio_wind_curtailed_on = 0
        
    if ( sum(cap_pv * df_data.pv_MW_array) - sum(df_out.pv_total) ) > 0:   
        curtailed_pv = (cap_pv * df_data.pv_MW_array) - df_out.pv_total
        ratio_pv_curtailed = curtailed_pv.sum() / sum(cap_pv * df_data.pv_MW_array)
    else:
        curtailed_pv = [0] * len(df_data)
        ratio_pv_curtailed = 0     
    
    # Only get the optimized CO2 costs when required
    an_op_co2 = an_op_co2.getValue() if euro_ton_co2!=0 else 0
    
    #Create overview of electricity loads
    overview_totals = pd.DataFrame({
                # Capacities
                 "cap_pv": cap_pv,
                 "cap_bat_en": cap_bat_en,
                 "cap_bat_p": cap_bat_p,   
                 "cap_wind_off": cap_wind_off,
                 "cap_wind_on": cap_wind_on,   
                 "cap_h2_ves": cap_h2_ves,
                 "cap_electrolyzer": cap_electrolyzer, 
                 "cap_fc": cap_fc,
                 "cap_solar_heat": cap_solar_heat,   
                 "cap_chp_mix": cap_chp_mix,
                 "cap_gasif": cap_gasif,   
                 "cap_gb_res": cap_gb_res,   
                 "cap_d_boiler": cap_d_boiler,
                 "cap_hp": cap_hp,
                 "cap_heat_storage":cap_heat_storage,
                 "cap_e_heating": cap_e_heating,
                 "cap_e_boiler":cap_e_boiler,
                 "p_grid_connection":p_grid_connection,
                "ghg_init":ghg_init,
                "cost_init":cost_init,
                "cost_reduction":percent_c_costs,
                "ghg_reduction":percent_c_ghgs,
                "total_time_h":total_time, #hours
                "mip_gap":m.MIPGap*100, # in %
                 "total_costs": obj1.getValue(), 
                 #Operation
                 "operation_costs": an_op.getValue(), 
                 "an_costs_op_wood": an_op_wood.getValue(), 
                 "an_costs_op_d_boiler": an_op_d_boiler.getValue(),
                 "an_costs_op_ng_res": an_op_ng_res.getValue(),
                 "an_costs_op_grid_abs": an_op_grid_abs.getValue(),  
                 "an_costs_op_grid_inj": -an_op_grid_inj.getValue(),
                 "an_costs_op_export_h2": -an_op_export_h2.getValue(),
                 "an_costs_op_co2": an_op_co2,
                 "an_costs_op_gasoline_cars":an_op_gasoline_cars.getValue(),
                 "an_costs_op_biogas": an_op_biogas.getValue(),

                 # Investment
                 "investment_costs": an_capex.getValue(), 
                 "an_costs_capex_boiler": an_capex_boiler.getValue(),
                 "an_costs_capex_hp": an_capex_hp.getValue(),
                 "an_costs_capex_pv": an_capex_pv.getValue(), 
                 "an_costs_capex_bat_en": an_capex_bat_en.getValue(),
                 "an_costs_capex_bat_p": an_capex_bat_p.getValue(),
                 "an_costs_capex_heat_storage": an_capex_heat_storage.getValue(),
                 "an_costs_capex_e_heating":an_capex_e_heating.getValue(),

                 "an_costs_capex_chp_mixer": an_capex_chp_mixer.getValue(),
                 "an_costs_capex_gasif": an_capex_gasif.getValue(),
                 "an_costs_capex_grid": an_capex_grid_ins.getValue(), 
                 "an_costs_capex_p_solar_heat": an_capex_p_solar_heat.getValue(),
                 "an_costs_capex_wind_on": an_capex_wind_on.getValue(),
                 "an_costs_capex_wind_off": an_capex_wind_off.getValue(),

                 "an_costs_capex_electrolyzer": an_capex_electrolyzer.getValue(),
                 "an_costs_capex_fc": an_capex_fc.getValue(),
        
                 "an_costs_capex_h2_ves": an_capex_h2_ves.getValue(), 
        
                 "an_costs_capex_gb_res": an_capex_gb_res.getValue(),  
                 "an_costs_capex_e_boiler": an_capex_e_boiler.getValue(),  
        
                 "an_costs_capex_trans_gas": an_capex_trans_gas.getValue(),  
                 "an_costs_capex_trans_bev":an_capex_trans_bev.getValue(),

                 # Replacement and O&M
                 "an_costs_rep": an_rep.getValue(), 
                 "an_costs_om": an_om.getValue(),

                 #GHG emissions
                 "total_ghg":obj2.getValue(), 
                 "operational_ghg": an_op_ghg.getValue(), 
                 "an_ghg_op_wood":an_ghg_op_wood.getValue(), 
                 "an_ghg_op_d_boiler":an_ghg_op_d_boiler.getValue(),         
                 "an_ghg_op_ng_res":an_ghg_op_ng_res.getValue(),
                 "an_ghg_op_grid_abs":an_ghg_op_grid_abs.getValue(),
                 "an_ghg_op_grid_inj": -an_ghg_op_grid_inj.getValue(),
                 "an_ghg_op_h2":-an_ghg_op_h2.getValue(),
                 "an_ghg_op_gasoline_cars": an_ghg_op_gasoline_cars.getValue(),
                 "an_ghg_op_chp":an_ghg_op_chp.getValue(),
        
                 "an_ghg_boiler": an_ghg_boiler.getValue(),
                 "an_ghg_hp": an_ghg_hp.getValue(),
                 "an_ghg_pv": an_ghg_pv.getValue(), 
                 "an_ghg_bat_en": an_ghg_bat_en.getValue(),
                 "an_ghg_bat_p": an_ghg_bat_p.getValue(),
                 "an_ghg_heat_storage": an_ghg_heat_storage.getValue(),
                 "an_ghg_chp_mixer": an_ghg_chp_mixer.getValue(),
                 "an_ghg_gasif": an_ghg_gasif.getValue(),
                 "an_ghg_grid_ins": an_ghg_grid_ins.getValue(), 
                 "an_ghg_p_solar_heat": an_ghg_p_solar_heat.getValue(),
                 "an_ghg_wind_on": an_ghg_wind_on.getValue(),
                 "an_ghg_wind_off": an_ghg_wind_off.getValue(),
                 "an_ghg_e_heating": an_ghg_e_heating.getValue(),
                 "an_ghg_electrolyzer": an_ghg_electrolyzer.getValue(),
                 "an_ghg_fc": an_ghg_fc.getValue(),
                 "an_ghg_h2_ves": an_ghg_h2_ves.getValue(),      
                 "an_ghg_e_boiler": an_ghg_e_boiler.getValue(),
                 "an_ghg_gb_res": an_ghg_gb_res.getValue(), 
                 "an_ghg_bev_cars":an_ghg_bev_cars.getValue(),
                    
                 "share_ghg_prod": an_op_ghg_inv.getValue()/obj2.getValue(),
                 "share_inv_total": an_capex.getValue()/obj1.getValue(),
        
                "ratio_pv_curtailed": ratio_pv_curtailed,
                "ratio_wind_curtailed_off": ratio_wind_curtailed_off,      
                "ratio_wind_curtailed_on": ratio_wind_curtailed_on,    
        
                 "share_bev": share_bev,
                 "share_gasoline":share_gasoline,
        
                 "ghg_grid_mean": df_data.ghg_impact.mean(),
                 "cost_grid_mean":df_data.grid_abs_price.mean(),   
                 "grid_elect_demand": sum(df_out.p_grid_abs) - sum(df_out.p_grid_inj), 
                 "f_biogas_waste_use":sum(df_out.f_chp_biogas),
                 "exported_h2":sum(df_out.export_h2),
                                   }, index=[ "opt_results_{}_{}_{}".format(ASSESSMENT_YEAR,w_cost,
                                                                                              round(eps_ghg_constraint,2))])
    if export_results:
        overview_totals.T.to_excel(r"results\opt_results_{}_{}_{}_{}_{}_{}_{}_{}.xlsx".format(ASSESSMENT_YEAR,w_cost,autonomous_elect,
                                                                                                                                              autonomous_gas,round(eps_ghg_constraint,2),euro_ton_co2,road_transport,
                                                                                                                                              export_alias))       
        df_out.to_excel(r"weekly_plots\result_{}_{}_{}_{}_{}_{}_{}_{}.xlsx".format(ASSESSMENT_YEAR,w_cost,autonomous_elect,
                                                                                                                                              autonomous_gas,round(eps_ghg_constraint,2),euro_ton_co2, road_transport,
                                                                                                                                                    export_alias))
    # Checks for simultaneous operations and raises an error if found.
    check_simultaneous_operations(df_out, "Simultaneous operation detected.")
    
    # Define input variables for LCA
    input_vars = {
                "parm": parm,
                "loc_elect": LOC_ELECT,
                "cap_wind_off": cap_wind_off,
                "cap_wind_on": cap_wind_on,
                "cap_pv": cap_pv,
                "cap_bat_en": cap_bat_en,
                "cap_bat_p": cap_bat_p,
                "cap_h2_ves": cap_h2_ves,
                "cap_electrolyzer": cap_electrolyzer,
                "cap_fc": cap_fc,
                "cap_hp": cap_hp,
                "cap_boiler": cap_d_boiler,
                "cap_boiler_res": cap_gb_res,
                "cap_heat_storage": cap_heat_storage,
                "cap_solar_heat": cap_solar_heat,
                "cap_chp_mix": cap_chp_mix,
                "cap_gasif": cap_gasif,
                "p_grid_connection": p_grid_connection,
                "cap_e_boiler": cap_e_boiler,
                "cap_e_heating": cap_e_heating,
                "summed_grid_abs": sum(df_out.p_grid_abs),
                "summed_grid_inj": sum(df_out.p_grid_inj),
                "summed_wood_syngas": sum(df_out.p_gasif),
                "summed_heat_boiler_ind": sum(df_out.p_d_boiler),
                "summed_heat_boiler_res": sum(df_out.p_gb_res),
                "summed_chp_fuel": sum(df_out.f_chp_mix),
                "ghgs_opt": obj2.getValue() * 1000, # it is in tonnes, convert to kg
                "w_cost": w_cost,
                "bev_car_kms": bev_car_kms.getValue(),
                "sec_db": sec_db,
                "petrol_car_kms": petrol_car_kms.getValue(),
                "h2_export": sum(df_out.export_h2),
                "credit_env_export": credit_env_export,
                "lcia_method": CC_METHOD,
                "epsilon_constraint": eps_ghg_constraint,
                "boiler_oil_res": False # not included in optimization problem.
            }

    # Call environmental_lca function with input variables
    if calc_all_lca_impacts:
        lca_results = environmental_lca(**input_vars)
    else:
        lca_results = ""       

    # Printing Cost Weight
    print(f"Cost weight '{w_cost}'")

    # Printing Annual Costs and GHG Emissions
    annual_costs_million_euro = round(obj1.getValue() / 1e3, 2)
    annual_ghg_emissions_kt_co2_eq = round(obj2.getValue() / 1e3, 2)
    print(f"Annual costs are '{annual_costs_million_euro}' M€")
    print(f"Annual GHG emissions are '{annual_ghg_emissions_kt_co2_eq}' kt CO2-eq.")

    # Separator Line
    print("***************************************************************")
    # Printing Percentage Changes and Reductions
    print(f"Percentual change costs = '{percent_c_costs}%'")
    print(f"Percentual change GHGs = '{percent_c_ghgs}%'")
    # Printing Calculation Time and MIP Gap
    calculation_time_hours = round(total_time, 2)
    mip_gap_percentage = round(m.MIPGap * 100, 2)
    print(f"Calculation time is '{calculation_time_hours}' hours")
    print(f"Final MIP gap value: '{mip_gap_percentage}%'")
    
    return overview_totals, lca_results, df_out

def calc_energy_system_results(df_data, parm, dict_ghg, loc_elect, total_hh, cap_wind_off, cap_wind_on, 
                      cap_pv, cap_bat_en, cap_bat_p, 
                      cap_h2_ves, cap_electrolyzer, cap_fc, 
                      cap_hp, cap_d_boiler, cap_gb_res, 
                      cap_heat_storage, cap_solar_heat, cap_chp_mix, 
                      cap_gasif, p_grid_connection, cap_e_boiler, cap_e_heating,
                     cap_boiler_oil_res,
                      summed_grid_abs, summed_grid_inj, summed_syngas, 
                      summed_heat_boiler_ind, summed_heat_boiler_res, summed_heat_oil_boiler_res,
                      summed_chp_fuel, summed_biogas, 
                      bev_car_kms, petrol_car_kms,
                      sec_db = NAME_REF_DB,
                      h2_export = 0, credit_env_export=True,
                               euro_ton_co2=0,
                      lcia_method=CC_METHOD, calc_all_lca_impacts=CALC_ALL_LCA_IMPACTS
                              ):
    
    """
    Calculates cost and environmental burdens of a non-optimized Multi-Energy System (MES) considering selected environmental impact categories.
    
    Args:
        df_data (pd.DataFrame): DataFrame with hourly input data.
        parm (dict): Dictionary of techno-economic cost data and assumptions [-].
        dict_ghg (dict): Dictionary with environmental emission data [kg CO2-eq./unit].
        loc_elect (str): Ecoinvent location of MES, abbreviation string (e.g., "GR") [str].
        total_hh (int): Total number of households in the system [int].
        cap_wind_off (float): Capacity of offshore wind [MW].
        cap_wind_on (float): Capacity of onshore wind [MW].
        cap_pv (float): Capacity of solar PV [MW].
        cap_bat_en (float): Energy capacity of battery electricity storage system [MWh].        
        cap_bat_p (float): Power capacity of battery electricity storage system [MW].
        cap_h2_ves (float): Energy capacity of hydrogen storage system [MWh].    
        cap_electrolyzer (float): Electrical capacity of electrolyzer [MWe].    
        cap_fc (float): Capacity of fuel cell [MW].
        cap_hp (float): Thermal output capacity of residential heat pump [MWth].
        cap_d_boiler (float): Fuel capacity of industrial boiler considered [MW].
        cap_gb_res (float): Fuel capacity of residential gas boiler considered [MW].
        cap_heat_storage (float): Energy capacity of residential heat storage system [MWh].
        cap_solar_heat (float): Installed area of residential thermal solar heat [ha].
        cap_chp_mix (float): Fuel capacity of advanced CHP unit [MW].
        cap_gasif (float): Fuel capacity of wood gasifier [MW].
        p_grid_connection (float): (Max) connection capacity of power grid [MW].
        cap_e_boiler (float): Fuel capacity of electric boiler for high temperature industrial heat [MW].
        cap_e_heating (float): Fuel capacity of residential electric heating [MW].
        cap_boiler_oil_res (float): Fuel capacity of an oil boiler for residential heating.
        summed_grid_abs (float): Grid absorption from power grid [MWh].
        summed_grid_inj (float): Grid injection to power grid [MWh].
        summed_syngas (float): Syngas produced [MWh].
        summed_heat_boiler_ind (float): Heat produced with industrial boiler [MWh].
        summed_heat_boiler_res (float): Heat produced with residential boiler [MWh].
        summed_heat_oil_boiler_res (float): Heat produced with residential oil boiler [MWh].
        summed_chp_fuel (float): Energy fuel input of advanced CHP unit [MWh].
        summed_biogas (float): Biogas produced [MWh].
        bev_car_kms (float): Kilometers driven in MES with Battery Electric Vehicles (BEVs) [kilometers].
        petrol_car_kms (float): Kilometers driven in MES with gasoline vehicles [kilometers].
        sec_db (str): Ecoinvent database used for calculation of LCA impacts [-].
        h2_export (float): Amount of exported hydrogen [MWh].
        credit_env_export (bool): If True, environmental credit is given for power injection and hydrogen export [-].
        euro_ton_co2 (float): CO2 price in Euro per ton [float].
        lcia_method (tuple): Standard climate change impact category used (e.g., CC_METHOD).
        calc_all_lca_impacts (bool, optional): whether to calculate other environmental impact categories. Default is True.
    Returns:
        lca_results (pd.DataFrame): DataFrame with environmental burdens of all selected environmental impact categories.
        overview_totals (pd.DataFrame): DataFrame with techno-economic information and other results.
    """

    """Objective 1: Costs"""
    """Operation (variable OPEX)"""
    an_op_wood = (summed_syngas/parm['gasif_eff']) * parm['wood_price'] # assumed with wood, but in actual optimization problem it can also be supplied with H2 and biogas.
    an_op_d_boiler = (summed_heat_boiler_ind/parm['d_boiler_eff']) * parm['diesel_price'] 
    an_op_ng_res = (summed_heat_boiler_res/parm['gb_eff_res']) * parm['ng_price'] 
    an_op_grid_abs = ( (df_data.electricity_demand * df_data.grid_abs_price).sum() # electricity
                    + (summed_grid_abs - df_data.electricity_demand.sum() ) * df_data.grid_abs_price.mean() # rest electricity for heating and electric vehicles
                     ) 
    an_op_grid_inj = (summed_grid_inj * df_data.rev_inj).sum() 
    an_op_export_h2 = h2_export * parm['h2_sell_price'] 
    an_op_biogas = summed_biogas * parm['biogas_price'] 
    an_op_oil_res =  (summed_heat_oil_boiler_res / parm['o_boiler_eff']) * parm['oil_price']

    share_gasoline = petrol_car_kms/(petrol_car_kms+bev_car_kms) if (petrol_car_kms+bev_car_kms)>0 else 0
    share_bev = 1-share_gasoline if (petrol_car_kms+bev_car_kms)>0 else 0
    
    # Add fuel costs for gasoline cars, if not based on BEVs only
    an_op_gasoline_cars = petrol_car_kms * parm['price_km_gasoline']

    an_op = (an_op_wood + an_op_d_boiler + an_op_ng_res + an_op_biogas +
                        an_op_grid_abs + an_op_gasoline_cars + an_op_oil_res - an_op_grid_inj -
                        an_op_export_h2
                       )

    """Investments (CAPEX)"""
    an_capex_boiler = calc_crf(parm['dr'],parm['project_lt']) * (cap_d_boiler * parm['d_boiler_capex'])   
    an_capex_gb_res = calc_crf(parm['dr'],parm['project_lt']) * (cap_gb_res * parm['gb_capex_res']) 
    an_capex_e_heating = calc_crf(parm['dr'],parm['project_lt']) * (cap_e_heating * parm['e_heating_capex']) 
    an_capex_boiler_oil_res = calc_crf(parm['dr'],parm['project_lt']) * (cap_boiler_oil_res * parm['o_boiler_capex']) 
    an_capex_hp = calc_crf(parm['dr'],parm['project_lt']) * (cap_hp * parm['hp_capex'])
    an_capex_electrolyzer = calc_crf(parm['dr'],parm['project_lt']) * (cap_electrolyzer * parm['electr_capex'])
    an_capex_fc = calc_crf(parm['dr'],parm['project_lt']) * (cap_fc * parm['fc_capex'])   
    an_capex_h2_ves = calc_crf(parm['dr'],parm['project_lt']) * (cap_h2_ves * parm['h2_ves_capex']) 
    an_capex_pv = calc_crf(parm['dr'],parm['project_lt']) * (cap_pv * parm['pv_capex']) 
    an_capex_wind_off = calc_crf(parm['dr'],parm['project_lt']) *(cap_wind_off * parm['wind_off_capex'])
    an_capex_wind_on = calc_crf(parm['dr'],parm['project_lt']) * (cap_wind_on * parm['wind_on_capex'])
    an_capex_bat_en = calc_crf(parm['dr'],parm['project_lt']) * (cap_bat_en * parm['bat_en_capex']) 
    an_capex_bat_p = calc_crf(parm['dr'],parm['project_lt']) * (cap_bat_p * parm['bat_p_capex'])
    an_capex_heat_storage = calc_crf(parm['dr'],parm['project_lt']) * (cap_heat_storage * parm['heat_stor_capex'])
    an_capex_p_solar_heat = calc_crf(parm['dr'],parm['project_lt']) * (cap_solar_heat * parm['solar_th_capex'])
    an_capex_chp_mixer = calc_crf(parm['dr'],parm['project_lt']) * (cap_chp_mix * parm['chp_mix_capex'] * parm['chp_mix_eff_e']) # The costs for the CHP unit is provided per MWe power delivered
    an_capex_gasif = calc_crf(parm['dr'],parm['project_lt']) * (cap_gasif * parm['gasif_capex'])
    an_capex_grid_ins = calc_crf(parm['dr'],parm['project_lt']) * (p_grid_connection * parm['grid_capex'])
    an_capex_e_boiler = calc_crf(parm['dr'],parm['project_lt']) * (cap_e_boiler * parm['e_boiler_capex'])   
    
    an_capex_trans_gas = share_gasoline * total_hh * calc_crf(parm['dr'],parm['project_lt']) * parm['gas_car_capex']
    an_capex_trans_bev = share_bev * total_hh * calc_crf(parm['dr'],parm['project_lt']) * parm['bev_car_capex']

    an_capex  = (an_capex_boiler + an_capex_gb_res + an_capex_e_heating + an_capex_boiler_oil_res 
                 + an_capex_hp + an_capex_electrolyzer +
                an_capex_fc + an_capex_h2_ves + 
                an_capex_pv + an_capex_wind_off + an_capex_wind_on +
                an_capex_bat_en + an_capex_bat_p + an_capex_heat_storage +
                an_capex_p_solar_heat + an_capex_chp_mixer +
                an_capex_gasif + an_capex_grid_ins + an_capex_e_boiler +
                an_capex_trans_gas + an_capex_trans_bev
               )

    """Replacements"""
    an_rep = (
          rep_annual_int((cap_d_boiler * parm['d_boiler_capex']),parm['project_lt'],parm['dr'], parm['d_boiler_lt']) 
        + rep_annual_int((cap_boiler_oil_res * parm['o_boiler_capex']),parm['project_lt'],parm['dr'], parm['o_boiler_lt']) 
        + rep_annual_int((cap_gb_res * parm['gb_capex_res']),parm['project_lt'],parm['dr'], parm['gb_lt_res'])  
        + rep_annual_int((cap_e_heating * parm['e_heating_capex']),parm['project_lt'],parm['dr'], parm['e_heating_lt'])  
        + rep_annual_int((cap_hp * parm['hp_capex']),parm['project_lt'],parm['dr'],parm['hp_lt'])
        + rep_annual_int((cap_h2_ves * parm['h2_ves_capex']),parm['project_lt'],parm['dr'],parm['h2_ves_lt'])
        + rep_annual_int((cap_electrolyzer * parm['electr_capex']),parm['project_lt'],parm['dr'],parm['electr_lt']) 
        + rep_annual_int((cap_fc * parm['fc_capex']),parm['project_lt'],parm['dr'],parm['fc_lt'])
        + rep_annual_int((cap_pv * parm['pv_capex']),parm['project_lt'],parm['dr'],parm['pv_lt'])
        + rep_annual_int((cap_wind_off * parm['wind_off_capex']),parm['project_lt'],parm['dr'],parm['wind_off_lt'])
        + rep_annual_int((cap_wind_on * parm['wind_on_capex']),parm['project_lt'],parm['dr'],parm['wind_on_lt'])
        + rep_annual_int((cap_bat_en * parm['bat_en_capex']),parm['project_lt'],parm['dr'],parm['bat_en_lt'])
        + rep_annual_int((cap_bat_p * parm['bat_p_capex']),parm['project_lt'],parm['dr'],parm['bat_p_lt'])
        + rep_annual_int((cap_heat_storage * parm['heat_stor_capex']),parm['project_lt'],parm['dr'],parm['heat_stor_lt'])
        + rep_annual_int((cap_solar_heat * parm['solar_th_capex']),parm['project_lt'],parm['dr'],parm['solar_th_lt'])
        + rep_annual_int((cap_chp_mix * parm['chp_mix_capex'] * parm['chp_mix_eff_e']),parm['project_lt'],parm['dr'],parm['chp_mix_lt'])# The costs for the CHP unit is provided per MWe power delivered
        + rep_annual_int((cap_gasif * parm['gasif_capex']),parm['project_lt'],parm['dr'],parm['gasif_lt'])
        + rep_annual_int((p_grid_connection * parm['grid_capex']),parm['project_lt'],parm['dr'],parm['grid_lt'])
        + rep_annual_int((cap_e_boiler * parm['e_boiler_capex']),parm['project_lt'],parm['dr'],parm['e_boiler_lt'])
        + rep_annual_int((share_gasoline * total_hh * parm['gas_car_capex']),parm['project_lt'],parm['dr'],parm['gas_car_lt'])
        + rep_annual_int((share_bev * total_hh * parm['bev_car_capex']),parm['project_lt'],parm['dr'],parm['bev_car_lt'])        
             )

    """Fixed O&M (without replacements)"""
    an_om = ((cap_d_boiler * parm['d_boiler_capex'] * parm['d_boiler_om']) 
            + (cap_boiler_oil_res * parm['o_boiler_capex'] * parm['o_boiler_om']) 
            + (cap_gb_res * parm['gb_capex_res'] * parm['gb_om_res']) 
            + (cap_e_heating * parm['e_heating_capex'] * parm['e_heating_om']) 
            + (cap_hp * parm['hp_capex'] * parm['hp_om']) 
            + (cap_h2_ves * parm['h2_ves_capex'] * parm['h2_ves_om'])          
            + (cap_electrolyzer * parm['electr_capex'] * parm['electr_om'])
            + (cap_fc * parm['fc_capex'] * parm['fc_om']) 
            + (cap_pv * parm['pv_capex'] * parm['pv_om']) 
            + (cap_wind_off * parm['wind_off_capex'] * parm['wind_off_om']) 
            + (cap_wind_on * parm['wind_on_capex'] * parm['wind_on_om']) 
            + (cap_bat_en * parm['bat_en_capex'] * parm['bat_en_om']) 
            + (cap_bat_p * parm['bat_p_capex'] * parm['bat_p_om'])
            + (cap_heat_storage * parm['heat_stor_capex'] * parm['heat_stor_om']) 
            + (cap_solar_heat * parm['solar_th_capex'] * parm['solar_th_om']) 
            + (cap_chp_mix * parm['chp_mix_capex'] * parm['chp_mix_om'] * parm['chp_mix_eff_e']) # The costs for the CHP unit is provided per MWe power delivered
            + (cap_gasif * parm['gasif_capex'] * parm['gasif_om']) 
            + (p_grid_connection * parm['grid_capex'] * parm['grid_om']) 
            + (cap_e_boiler * parm['e_boiler_capex'] * parm['e_boiler_om']) 
             
            + (share_gasoline * parm['gas_car_capex'] * total_hh * parm['gas_car_om']) 
            + (share_bev * parm['bev_car_capex'] * total_hh * parm['bev_car_om']) 
            )
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    """Total costs, including possible costs for CO2"""
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
    an_op_co2 = obj2 * (euro_ton_co2/1000) if euro_ton_co2!=0 else 0  
    obj1 = an_op + an_capex + an_rep + an_om + an_op_co2
    
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    """Objective 2: life cycle GHG emissions"""
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    """Operation"""
    
    # Note: biogas is assumed to be burden-free
    an_ghg_op_wood = summed_syngas * dict_ghg['ghg_imp_syngas'] 
    #dataset of heat relates to fuel input, so divide by efficiency
    an_ghg_op_d_boiler = summed_heat_boiler_ind * dict_ghg['ghg_imp_diesel_ht']
    an_ghg_op_oil_res = summed_heat_oil_boiler_res * dict_ghg['ghg_imp_oil_heating'] # Construction is included for this unit.
    
    an_ghg_op_ng_res = summed_heat_boiler_res * dict_ghg['ghg_imp_ng'] 
    an_ghg_op_h2 = h2_export * dict_ghg['ghg_imp_smr'] if credit_env_export else 0
    an_ghg_op_chp = summed_chp_fuel * dict_ghg['ghg_imp_chp_fuel'] 
    an_ghg_op_grid_abs = summed_grid_abs * df_data.ghg_impact.mean() # grid GHG intensity is here assumed stable. 
    an_ghg_op_grid_inj = (summed_grid_inj * df_data.ghg_impact_cons).sum() if credit_env_export else 0
    an_ghg_op_gasoline_cars = petrol_car_kms * dict_ghg['ghg_gas_car_km']
    
    an_op_ghg = (an_ghg_op_wood +
                 an_ghg_op_d_boiler +
                 an_ghg_op_oil_res +
                 an_ghg_op_ng_res +
                 an_ghg_op_h2 +
                 an_ghg_op_chp +
                 an_ghg_op_grid_abs -
                 an_ghg_op_grid_inj +
                 an_ghg_op_gasoline_cars)

    """Production and replacements"""    
    an_ghg_boiler = (cap_d_boiler * dict_ghg['ghg_imp_gb'] * (parm['project_lt']/parm['d_boiler_lt'])) / parm['project_lt']  
    an_ghg_gb_res = (cap_gb_res * dict_ghg['ghg_imp_gb'] * (parm['project_lt']/parm['gb_lt_res'])) / parm['project_lt']
    an_ghg_e_heating = (cap_e_heating * dict_ghg['ghg_imp_e_heating'] * (parm['project_lt']/parm['e_heating_lt'])) / parm['project_lt']  
    an_ghg_hp = (cap_hp * dict_ghg['ghg_imp_hp'] * (parm['project_lt']/parm['hp_lt'])) / parm['project_lt']  
    an_ghg_electrolyzer = (cap_electrolyzer * dict_ghg['ghg_imp_electr'] * (parm['project_lt']/parm['electr_lt'])) / parm['project_lt']  
    an_ghg_fc = parm['fc_eff_e'] * (cap_fc * dict_ghg['ghg_imp_fc'] * (parm['project_lt']/parm['fc_lt'])) / parm['project_lt']     
    an_ghg_h2_ves = (cap_h2_ves * dict_ghg['ghg_imp_h2_ves'] * (parm['project_lt']/parm['h2_ves_lt'])) / parm['project_lt']     
    an_ghg_pv = (cap_pv * dict_ghg['ghg_imp_pv'] * (parm['project_lt']/parm['pv_lt'])) / parm['project_lt']  
    an_ghg_wind_off = (cap_wind_off * dict_ghg['ghg_imp_wind_off'] * (parm['project_lt']/parm['wind_off_lt'])) / parm['project_lt']  
    an_ghg_wind_on = (cap_wind_on * dict_ghg['ghg_imp_wind_on'] * (parm['project_lt']/parm['wind_on_lt'])) / parm['project_lt']   
    an_ghg_bat_en = (cap_bat_en * dict_ghg['ghg_imp_bat_cap'] * (parm['project_lt']/parm['bat_en_lt'])) / parm['project_lt']  
    an_ghg_bat_p = (cap_bat_p * dict_ghg['ghg_imp_bat_p'] * (parm['project_lt']/parm['bat_p_lt'])) / parm['project_lt']  
    an_ghg_heat_storage = (cap_heat_storage * dict_ghg['ghg_imp_heat_storage'] * (parm['project_lt']/parm['heat_stor_lt'])) / parm['project_lt']  
    an_ghg_p_solar_heat = (cap_solar_heat * dict_ghg['ghg_imp_solar_heat_ha'] * (parm['project_lt']/parm['solar_th_lt'])) / parm['project_lt']  
    an_ghg_chp_mixer = (cap_chp_mix * dict_ghg['ghg_imp_chp_mix'] * (parm['project_lt']/parm['chp_mix_lt'])) / parm['project_lt']   # The costs for the CHP unit is provided per MWe power delivered
    an_ghg_gasif = (cap_gasif * dict_ghg['ghg_imp_gasif'] * (parm['project_lt']/parm['gasif_lt'])) / parm['project_lt']  
    an_ghg_grid_ins = (p_grid_connection * dict_ghg['ghg_impact_grid_network'] * (parm['project_lt']/parm['grid_lt'])) / parm['project_lt']
    an_ghg_e_boiler = (cap_e_boiler * dict_ghg['ghg_imp_e_heating'] * (parm['project_lt']/parm['e_boiler_lt'])) / parm['project_lt']  
    
    an_ghg_bev_cars = bev_car_kms * dict_ghg['ghg_bev_car_km'] 
    
    an_op_ghg_inv = (an_ghg_boiler + 
                    an_ghg_gb_res + an_ghg_e_heating + an_ghg_hp + an_ghg_electrolyzer +
                    an_ghg_fc + an_ghg_h2_ves +
                    an_ghg_pv +  an_ghg_wind_off + an_ghg_wind_on +
                    an_ghg_bat_en + an_ghg_bat_p + an_ghg_heat_storage +
                    an_ghg_p_solar_heat + an_ghg_chp_mixer + 
                    an_ghg_gasif + an_ghg_grid_ins + an_ghg_e_boiler + an_ghg_bev_cars
               )

    obj2 = an_op_ghg + an_op_ghg_inv
    
    #Create overview of electricity loads
    overview_totals = pd.DataFrame({
                # Capacities
                 "cap_pv": cap_pv,
                 "cap_bat_en": cap_bat_en,
                 "cap_bat_p": cap_bat_p,   
                 "cap_wind_off": cap_wind_off,
                 "cap_wind_on": cap_wind_on,   
                 "cap_h2_ves": cap_h2_ves,
                 "cap_electrolyzer": cap_electrolyzer, 
                 "cap_fc": cap_fc,
                 "cap_solar_heat": cap_solar_heat,   
                 "cap_chp_mix": cap_chp_mix,
                 "cap_gasif": cap_gasif,   
                 "cap_gb_res": cap_gb_res,   
                 "cap_d_boiler": cap_d_boiler,
                 "cap_hp": cap_hp,
                 "cap_heat_storage":cap_heat_storage,
                 "cap_e_heating": cap_e_heating,
                 "cap_e_boiler":cap_e_boiler,
                 "p_grid_connection":p_grid_connection,
                 'cap_boiler_oil_res': cap_boiler_oil_res,
                 "total_costs": obj1, 
                 
                 #Operation
                 "operation_costs": an_op, 
                 "an_costs_op_wood": an_op_wood, 
                 "an_costs_op_d_boiler": an_op_d_boiler,
                 "an_costs_op_ng_res": an_op_ng_res,
                 "an_costs_op_grid_abs": an_op_grid_abs,  
                 "an_costs_op_grid_inj": -an_op_grid_inj,
                 "an_costs_op_export_h2": -an_op_export_h2,
                 "an_costs_op_co2": an_op_co2,
                 "an_costs_op_gasoline_cars":an_op_gasoline_cars,
                 "an_costs_op_biogas": an_op_biogas,

                 # Investment
                 "investment_costs": an_capex, 
                 "an_costs_capex_boiler": an_capex_boiler,
                 "an_costs_capex_hp": an_capex_hp,
                 "an_costs_capex_pv": an_capex_pv, 
                 "an_costs_capex_bat_en": an_capex_bat_en,
                 "an_costs_capex_bat_p": an_capex_bat_p,
                 "an_costs_capex_heat_storage": an_capex_heat_storage,
                 "an_costs_capex_e_heating":an_capex_e_heating,

                 "an_costs_capex_chp_mixer": an_capex_chp_mixer,
                 "an_costs_capex_gasif": an_capex_gasif,
                 "an_costs_capex_grid": an_capex_grid_ins, 
                 "an_costs_capex_p_solar_heat": an_capex_p_solar_heat,
                 "an_costs_capex_wind_on": an_capex_wind_on,
                 "an_costs_capex_wind_off": an_capex_wind_off,

                 "an_costs_capex_electrolyzer": an_capex_electrolyzer,
                 "an_costs_capex_fc": an_capex_fc,
        
                 "an_costs_capex_h2_ves": an_capex_h2_ves, 
        
                 "an_costs_capex_gb_res": an_capex_gb_res,  
                 "an_costs_capex_e_boiler": an_capex_e_boiler,  
        
                 "an_costs_capex_trans_gas": an_capex_trans_gas,  
                 "an_costs_capex_trans_bev":an_capex_trans_bev,

                 # Replacement and O&M
                 "an_costs_rep": an_rep, 
                 "an_costs_om": an_om,

                 #GHG emissions
                 "total_ghg":obj2, 
                 "operational_ghg": an_op_ghg, 
                 "an_ghg_op_wood":an_ghg_op_wood, 
                 "an_ghg_op_d_boiler":an_ghg_op_d_boiler,         
                 "an_ghg_op_ng_res":an_ghg_op_ng_res,
                 "an_ghg_op_oil_res":an_ghg_op_oil_res,
                 "an_ghg_op_grid_abs":an_ghg_op_grid_abs,
                 "an_ghg_op_grid_inj": -an_ghg_op_grid_inj,
                 "an_ghg_op_h2":-an_ghg_op_h2,
                 "an_ghg_op_gasoline_cars": an_ghg_op_gasoline_cars,
                 "an_ghg_op_chp":an_ghg_op_chp,
        
                 "an_ghg_boiler": an_ghg_boiler,
                 "an_ghg_hp": an_ghg_hp,
                 "an_ghg_pv": an_ghg_pv, 
                 "an_ghg_bat_en": an_ghg_bat_en,
                 "an_ghg_bat_p": an_ghg_bat_p,
                 "an_ghg_heat_storage": an_ghg_heat_storage,

                 "an_ghg_chp_mixer": an_ghg_chp_mixer,
                 "an_ghg_gasif": an_ghg_gasif,
                 "an_ghg_grid_ins": an_ghg_grid_ins, 
                 "an_ghg_p_solar_heat": an_ghg_p_solar_heat,
                 "an_ghg_wind_on": an_ghg_wind_on,
                 "an_ghg_wind_off": an_ghg_wind_off,
                 "an_ghg_e_heating": an_ghg_e_heating,

                 "an_ghg_electrolyzer": an_ghg_electrolyzer,
                 "an_ghg_fc": an_ghg_fc,
        
                 "an_ghg_h2_ves": an_ghg_h2_ves,      
                 "an_ghg_e_boiler": an_ghg_e_boiler,
        
                 "an_ghg_gb_res": an_ghg_gb_res, 
                 "an_ghg_bev_cars":an_ghg_bev_cars,
        
                 "share_bev": share_bev,
                 "share_gasoline":share_gasoline,
                                   }, index=[ "opt_results_{}_{}".format(loc_elect,sec_db)])
    
    if calc_all_lca_impacts:
        lca_results = environmental_lca(
                        parm=parm,
                        loc_elect=loc_elect,
                        cap_wind_off=cap_wind_off,
                        cap_wind_on=cap_wind_on,
                        cap_pv=cap_pv,
                        cap_bat_en=cap_bat_en,
                        cap_bat_p=cap_bat_p,
                        cap_h2_ves=cap_h2_ves,
                        cap_electrolyzer=cap_electrolyzer,
                        cap_fc=cap_fc,
                        cap_hp=cap_hp,
                        cap_boiler=cap_d_boiler,
                        cap_boiler_res=cap_gb_res,
                        cap_heat_storage=cap_heat_storage,
                        cap_solar_heat=cap_solar_heat,
                        cap_chp_mix=cap_chp_mix,
                        cap_gasif=cap_gasif,
                        p_grid_connection=p_grid_connection,
                        cap_e_boiler=cap_e_boiler,
                        cap_e_heating=cap_e_heating,
                        summed_grid_abs=summed_grid_abs,
                        summed_grid_inj=summed_grid_inj,
                        summed_wood_syngas=summed_syngas,
                        summed_heat_boiler_ind=summed_heat_boiler_ind,
                        summed_heat_boiler_res=summed_heat_boiler_res,
                        summed_chp_fuel=summed_chp_fuel,
                        ghgs_opt=obj2 * 1000,
                        w_cost=1,
                        bev_car_kms=bev_car_kms,
                        petrol_car_kms=petrol_car_kms,
                        h2_export=h2_export,
                        boiler_oil_res=summed_heat_oil_boiler_res
                    )
    else:
        lca_results = ""    
    
    return overview_totals, lca_results

# # This equation is needed to calculate environmental burdens, other than climate change
def environmental_lca(parm, loc_elect, cap_wind_off, cap_wind_on, 
                      cap_pv, cap_bat_en, cap_bat_p, 
                      cap_h2_ves, cap_electrolyzer, cap_fc, 
                      cap_hp, cap_boiler, cap_boiler_res, 
                      cap_heat_storage, cap_solar_heat, cap_chp_mix, 
                      cap_gasif, p_grid_connection, cap_e_boiler, cap_e_heating,
                      summed_grid_abs, summed_grid_inj, summed_wood_syngas, 
                      summed_heat_boiler_ind, summed_heat_boiler_res, 
                      summed_chp_fuel, ghgs_opt, w_cost, bev_car_kms,
                      sec_db = NAME_REF_DB,
                      petrol_car_kms = 0, h2_export = 0, credit_env_export=True,
                      lcia_method=CC_METHOD, 
                      epsilon_constraint = False, boiler_oil_res=False):
    
    """
    Calculates environmental burdens of a MES, including all selected environmental impact categories. Returns a dataframe with all this information.

    Args:
        parm (dict): dictionary of techno-ecomic cost data and assumptions [-].    
        loc_elect (str): ecoinvent location of MES, abbreviation string, e.g., "GR" [str].
        cap_wind_off (float): capacity of offshore wind [MW].
        cap_wind_on (float): capacity of onshore wind [MW].
        cap_pv (float): capacity of solar PV [MW].
        cap_bat_en (float): energy capacity of battery electricity storage system [MWh].        
        cap_bat_p (float): power capacity of battery electricity storage system [MW].
        cap_h2_ves (float): energy capacity of hydrogen storage system [MWh].    
        cap_electrolyzer (float): electrical capacity electrolyzer [MWe].    
        cap_fc (float): capacity of fuel cell [MW].
        cap_hp (float): thermal output capacity of residential heat pump [MWth].
        cap_boiler (float): fuel capacity of boiler considered [MW].
        cap_boiler_res (float): fuel capacity of residential boiler considered [MW].
        cap_heat_storage (float): energy capacity of residential heat storage system [MWh].
        cap_solar_heat (float): installed area of residential thermal solar heat [ha].
        cap_chp_mix (float): fuel capacity of advanced CHP unit [MW].
        cap_gasif (float): fuel capacity of wood gasifier [MW].
        p_grid_connection (float): (max) connection of power grid [MW].
        cap_e_boiler (float): fuel capacity of electric boiler for high temperature industrial heat [MW].
        cap_e_heating (float): fuel capacity of residential electric heating [MW].
        summed_grid_abs (float): grid absorption from power grid [MWh].
        summed_grid_inj (float): grid injection to power grid [MWh].
        summed_wood_syngas (float): syngas produced [MWh].
        summed_heat_boiler_ind (float): heat produced with industrial boiler [MWh].
        summed_heat_boiler_res (float): heat produced with residential boiler [MWh].
        summed_chp_fuel (float): energy fuel input of advanced CHP unit [MWh].
        ghgs_opt (float): GHG emissions as a results from the optimization, this serves as input here to check whether outcomes are the identical [kg CO2-eq.].
        w_cost (float): weight of cost objective (between 0 and 1) [-].
        bev_car_kms (float): kilometers driven in MES with BEVs [kilometers].
        sec_db (str): ecoinvent database used for calculation of LCA impacts [-].
        petrol_car_kms (float): kilometers driven in MES with gasoline vehicles [kilometers].
        h2_export (float): amount of exported hydrogen [MWh].
        credit_env_export (bool): if True, environmental credit is given for power injection and hydrogen export [-].
        lcia_method (str): standard climate change impact category used [-].
        epsilon_constraint (str): If true, then there is a constraint on GHG emissions.
        boiler_oil_res (bool or float): If true, then a oil boiler is used for residential heating.
    Returns:
        lca_results (pd.DataFrame): dataframe with environmental burdens of all selected environmental impact categories.
    """
    
    # Make activity:
    db_name = 'island_energy_system'
    
    exchanges = []
    
    # Offshore wind, per MW
    if cap_wind_off>0:
        exch_wind_off_1 = {'name': "market for wind power plant, 2MW, offshore, fixed parts",
                'reference product': "wind power plant, 2MW, offshore, fixed parts",
                'database': sec_db,
                'location': "GLO",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': (cap_wind_off/2 * (parm['project_lt']/parm['wind_off_lt']))/parm['project_lt']}
        exchanges.append(exch_wind_off_1)

        exch_wind_off_2 = {'name': "market for wind power plant, 2MW, offshore, moving parts",
                'reference product': "wind power plant, 2MW, offshore, moving parts",
                'database': sec_db,
                'location': "GLO",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': (cap_wind_off/2 * (parm['project_lt']/parm['wind_off_lt']))/parm['project_lt']}
        exchanges.append(exch_wind_off_2)
    
    # Onshore wind, per MW
    if cap_wind_on > 0:
        exch_wind_on = {'name': "market for wind turbine, 2MW, onshore",
                'reference product': "wind turbine, 2MW, onshore",
                'database': sec_db,
                'location': "GLO",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': (cap_wind_on/2 * (parm['project_lt']/parm['wind_on_lt']))/parm['project_lt']}
        exchanges.append(exch_wind_on)

        exch_wind_on_network = {'name': "market for wind turbine network connection, 2MW, onshore",
                'reference product': "wind turbine network connection, 2MW, onshore",
                'database': sec_db,
                'location': "GLO",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': (cap_wind_on/2 * (parm['project_lt']/parm['wind_on_lt']))/parm['project_lt']}
        exchanges.append(exch_wind_on_network)
        
    # Solar PV
    if cap_pv>0:    
        exch_pv = {'name': "photovoltaic open ground installation, 570 kWp, multi-Si, on open ground",
                   "reference product": "photovoltaic open ground installation, 570 kWp, multi-Si, on open ground",
                'database': sec_db,
                'location': "RER",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': (cap_pv/(0.570*0.895) * (parm['project_lt']/parm['pv_lt']))/parm['project_lt']}
        exchanges.append(exch_pv)
    
    # Battery - Energy unit
    if cap_bat_en>0:
        exch_bat_bms = {'name': "battery management system, kWh",
                'reference product': "battery management system, kWh",
                'database': sec_db,
                'location': "GLO",
                'type': 'technosphere',
                'unit': 'kilowatt hour',
                'amount': (cap_bat_en * 1e3 * (parm['project_lt']/parm['bat_en_lt']))/parm['project_lt']}
        exchanges.append(exch_bat_bms)

        exch_bat_ems = {'name': "energy management system, kWh",
                'reference product': "energy management system, kWh",
                'database': sec_db,
                'location': "GLO",
                'type': 'technosphere',
                'unit': 'kilowatt hour',
                'amount': (cap_bat_en * 1e3 * (parm['project_lt']/parm['bat_en_lt']))/parm['project_lt']}
        exchanges.append(exch_bat_ems)

        exch_bat_cap = {'name': "li-ion (NMC)",
                'reference product': "li-ion (NMC)",
                'database': sec_db,
                'location': "GLO",
                'type': 'technosphere',
                'unit': 'kilowatt hour',
                'amount': (cap_bat_en * 1e3 * (parm['project_lt']/parm['bat_en_lt']))/parm['project_lt']}
        exchanges.append(exch_bat_cap)

    # Battery - Power unit
    if cap_bat_p>0:
        exch_bat_pcs = {'name': "power conditioning system, container system",
                'reference product': "power conditioning system, container system",
                'database': sec_db,
                'location': "GLO",
                'type': 'technosphere',
                'unit': 'kW',
                'amount': (cap_bat_p * 1e3 * (parm['project_lt']/parm['bat_p_lt']))/parm['project_lt']}
        exchanges.append(exch_bat_pcs)
    
    # Hydrogen storage
    if cap_h2_ves>0:
        exch_stor_ves = {'name': "high pressure hydrogen storage tank",
            'reference product': "high pressure hydrogen storage tank",
            'database': sec_db,
            'location': "GLO",
            'type': 'technosphere',
            'unit': 'kilogram',
            'amount': ( 1e3 * (cap_h2_ves/(MJ_KG_H2/3.6)) * (parm['project_lt']/parm['h2_ves_lt']))/parm['project_lt']}
        exchanges.append(exch_stor_ves)
    
    # Electrolyzer
    if cap_electrolyzer>0:
        # replacements are already accounted for in these LCIs, so there are set to BoS lifetime (20 years)
        exch_electr = {'name': "electrolyzer production, 1MWe, PEM, Stack",
                'reference product': "electrolyzer, 1MWe, PEM, Stack",
                'database': sec_db,
                'location': "RER",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': ( cap_electrolyzer * (parm['project_lt']/parm['electr_lt']))/parm['project_lt']}
        exchanges.append(exch_electr)
        
        exch_electr_bop = {'name': "electrolyzer production, 1MWe, PEM, Balance of Plant",
                'reference product': "electrolyzer, 1MWe, PEM, Balance of Plant",
                'database': sec_db,
                'location': "RER",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': ( cap_electrolyzer * (parm['project_lt']/parm['electr_bos_lt']))/parm['project_lt']}
        exchanges.append(exch_electr_bop)
    
    # Fuel cell
    if cap_fc:
        exch_fc = {'name': "fuel cell system assembly, 1 kWe, proton exchange membrane (PEM)",
                'reference product': "fuel cell system, 1 kWe, proton exchange membrane (PEM)",
                'database': sec_db,
                'location': "GLO",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': (parm['fc_eff_e']*(cap_fc/0.001) * (parm['project_lt']/parm['fc_lt']))/parm['project_lt']}
        exchanges.append(exch_fc)    
    
    # Heat pump
    if cap_hp>0:
        exch_hp = {'name': "heat pump, air-water, 7kW",
                'reference product': "heat pump, air-water, 7kW",
                'database': sec_db,
                'location': "RER",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': ((cap_hp/(0.007*AV_COP_ASHP) ) * (parm['project_lt']/parm['hp_lt']))/parm['project_lt']}
        exchanges.append(exch_hp)  
    
    # Gas boiler
    if cap_boiler>0:
        exch_gb = {'name': "gas boiler production",
                'reference product': "gas boiler",
                'database': sec_db,
                'location': "RoW",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': ( (cap_boiler/0.01) * (parm['project_lt']/parm['d_boiler_lt']))/parm['project_lt']}
        exchanges.append(exch_gb) 
    
    # Gas boiler - Residential
    if cap_boiler_res>0:
        exch_gb_res = {'name': "gas boiler production",
                'reference product': "gas boiler",
                'database': sec_db,
                'location': "RoW",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': ((cap_boiler_res/0.01) * (parm['project_lt']/parm['gb_lt_res']))/parm['project_lt']}
        exchanges.append(exch_gb_res) 
        
    # Electric boiler - Industry
    if cap_e_boiler>0:      
        exch_cap_e_boiler = {'name': "electric storage heater, 5kW",
                'reference product': "electric storage heater, 5kW",
                'database': sec_db,
                'location': "CH",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': ( (cap_e_boiler/0.005) * (parm['project_lt']/parm['e_boiler_lt']))/parm['project_lt'] }#per 5 kW
        exchanges.append(exch_cap_e_boiler)    
    
    # Heat storage
    MWh_water = (4.2 * (parm['heat_stor_t_max'] - parm['heat_stor_t_min']) * 2000)/ 3600 / 1e3 #https://www.engineeringtoolbox.com/energy-storage-water-d_1463.html
    
    if cap_heat_storage>0:
        exch_hs = {'name': "heat storage production, 2000l",
                'reference product': "heat storage, 2000l",
                'database': sec_db,
                'location': "RoW",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': ((cap_heat_storage/MWh_water) * (parm['project_lt']/parm['heat_stor_lt']))/parm['project_lt']}
        exchanges.append(exch_hs) 
    
    # Solar heat
    if cap_solar_heat>0:
        exch_sh = {'name': "solar collector system installation, Cu flat plate collector, one-family house, combined system, w\o storage",
                'reference product': "solar collector system, Cu flat plate collector, one-family house, combined system",
                'database': sec_db,
                'location': "CH",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': ( 10000 * (cap_solar_heat/12.3) * (parm['project_lt']/parm['solar_th_lt']))/parm['project_lt']} # provided in hectares, so multiply with factor 100m*100m here.
        exchanges.append(exch_sh)  
    
    # CHP mixer
    if cap_chp_mix>0:
        exch_chp_mix = {'name': "advanced gas turbine, 400kWe",
                'reference product': "advanced gas turbine, 400kWe",
                'database': sec_db,
                'location': "FI",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': parm['chp_mix_eff_e'] * (cap_chp_mix/0.400 * (parm['project_lt']/parm['chp_mix_lt']))/parm['project_lt']}
        exchanges.append(exch_chp_mix)  
        
        exch_chp_fuel = {'name': "mixed CHP unit, operation",
                'reference product': "mixed CHP unit, operation",
                'database': sec_db,
                'location': "RER",
                'type': 'technosphere',
                'unit': 'kilowatt hour',
                'amount': summed_chp_fuel * 1e3}
        exchanges.append(exch_chp_fuel) 
    
    # Gasifier
    if cap_gasif>0:
        exch_gasif = {'name': "synthetic gas factory construction",
                'reference product': "synthetic gas factory",
                'database': sec_db,
                'location': "RoW",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': (cap_gasif/7.4 * (parm['project_lt']/parm['gasif_lt']))/parm['project_lt']}
        exchanges.append(exch_gasif)  
        
    if summed_wood_syngas>0:
        exch_syngas = {'name': 'synthetic gas production, from wood, at fixed bed gasifier, w\o infrastructure',
                'reference product': 'synthetic gas',
                'database': sec_db,
                'location': loc_elect,
                'type': 'technosphere',
                'unit': 'cubic meter',
                'amount': 1e3 * summed_wood_syngas / (5.2 / 3.6)}
        exchanges.append(exch_syngas)
    
    # Grid electricity
    if summed_grid_abs>0:
        exch_grid_elect = {'name': "market for electricity, low voltage, Crete",
                'reference product': "electricity, low voltage",
                'database': sec_db,
                'location': loc_elect,
                'type': 'technosphere',
                'unit': 'kilowatt hour',
                'amount': 1e3 * summed_grid_abs}
        exchanges.append(exch_grid_elect)
    
    if credit_env_export:
        # Grid electricity - CONSEQUENTIAL FOR GRID INJECTION
        exch_grid_elect_cons = {'name': "market for electricity, low voltage",
                'reference product': "electricity, low voltage",
                'database': DB_NAME_CONS,
                'location': loc_elect,
                'type': 'technosphere',
                'unit': 'kilowatt hour',
                'amount': -summed_grid_inj * 1e3}
        exchanges.append(exch_grid_elect_cons)
    
    # Natural gas
    if summed_heat_boiler_ind>0:
        exch_db = {'name': 'heat production, light fuel oil, at industrial furnace 1MW',
                'reference product': 'heat, district or industrial, other than natural gas',
                'database': sec_db,
                'location': "Europe without Switzerland",
                'type': 'technosphere',
                'unit': 'megajoule',
                'amount': summed_heat_boiler_ind * 3.6 * 1e3}
        exchanges.append(exch_db)
    
    # Natural gas - Residential
    if summed_heat_boiler_res>0:
        exch_ng_res = {'name': 'heat production, natural gas, w\o infrastructure',
                'reference product': 'heat production, natural gas, w\o infrastructure',
                'database': sec_db,
                'location': "Europe without Switzerland",
                'type': 'technosphere',
                'unit': 'megajoule',
                'amount': summed_heat_boiler_res * 3.6 * 1e3}
        exchanges.append(exch_ng_res)
    
    if p_grid_connection>0:
        exch_impact_grid_network = {'name': "wind turbine network connection construction, 4.5MW, onshore",
                'reference product': "wind turbine network connection, 4.5MW, onshore",
                'database': sec_db,
                'location': "GLO",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': ( (p_grid_connection/4.5) * (parm['project_lt']/parm['grid_lt']))/parm['project_lt']}
        exchanges.append(exch_impact_grid_network)
    
    if (h2_export>0 and credit_env_export):
        exch_impact_h2_export = {'name': "hydrogen production, steam reforming",
                'reference product': "hydrogen, gaseous",
                'database': sec_db,
                'location': "RER",
                'type': 'technosphere',
                'unit': 'kilogram',
                'amount': - (h2_export / (MJ_KG_H2/3.6)) * 1e3} #kWh to kg
        exchanges.append(exch_impact_h2_export)   
    
    if petrol_car_kms>0:
        exch_impact_car = {'name': "market for transport, passenger car, medium size, petrol, EURO 5",
                'reference product': "transport, passenger car, medium size, petrol, EURO 5",
                'database': sec_db,
                'location': "GLO",
                'type': 'technosphere',
                'unit': 'kilometer',
                'amount': petrol_car_kms } #kilometers
        exchanges.append(exch_impact_car)    
    
    if bev_car_kms>0:
        exch_impact_bev_car = {'name': "transport, passenger car, electric, w\o fuel",
                'reference product': "transport, passenger car, electric, w\o fuel",
                'database': sec_db,
                'location': "GLO",
                'type': 'technosphere',
                'unit': 'kilometer',
                'amount': bev_car_kms } #kilometers
        exchanges.append(exch_impact_bev_car) 
        
    if cap_e_heating:
        exch_impact_e_heating = {'name': "electric storage heater, 5kW",
                'reference product': "electric storage heater, 5kW",
                'database': sec_db,
                'location': "CH",
                'type': 'technosphere',
                'unit': 'unit',
                'amount': ( (cap_e_heating/0.005) * (parm['project_lt']/parm['e_heating_lt']))/parm['project_lt'] }#per 5 kW
        exchanges.append(exch_impact_e_heating)    
        
    if boiler_oil_res>0:
        exch_impact_oil_heating = {'name': "heat production, light fuel oil, at boiler 10kW condensing, non-modulating",
                'reference product': "heat, central or small-scale, other than natural gas",
                'database': sec_db,
                'location': "Europe without Switzerland",
                'type': 'technosphere',
                'unit': 'megajoule',
                'amount': (boiler_oil_res * 3.6 * 1e3) }#per MWh to megajoule
        exchanges.append(exch_impact_oil_heating)    
    
    create_process_and_add(loc_elect, ASSESSMENT_YEAR, exchanges, sec_db)

    activities = list(bw.Database(DB_NAME))
    activities = [k for k in activities
                        if str('multi_energy_system') in str(k)]

    lca = bw.LCA({activities[0]: 1, activities[-1]: 1}, method=lcia_method)
    lca.lci()
    lca.lcia()
    
    if "climate change" in str(lcia_method):
        if abs((ghgs_opt-lca.score))>1:
            print("**************************************************")
            print(abs((ghgs_opt-lca.score)))
            print("Difference between scores, inititial calc score is '{}' and LCA score here is '{}'".format(ghgs_opt,lca.score))
            #report results for debugging, per exchange:
            mes = activities[0]
            lca.redo_lcia({mes: 1})
            for exc in mes.exchanges():
                if exc['type'] == 'technosphere':
                    lca.redo_lcia({exc.input: exc['amount']})
                    print("{}, amount: '{}', lca results: '{}'".format(exc['name'], 
                                                                       exc['amount'], lca.score))
                elif exc['type'] == 'biosphere':
                    # Need to multiple the amount times its CF
                    cf = lca.characterization_matrix[lca.biosphere_dict[exc.input], :].sum()
                    print("{}, amount: '{}', lca results: '{}'".format(exc['name'], 
                                                                       exc['amount'], cf * exc['amount']))        
            raise ValueError("ERROR: in calculation please check GHG calculation")

    # This part of the script is adopted from Antonini and Treyer et al. (2020)
    result_array = [[[] for _ in activities] for _ in my_methods]

    for i, method in enumerate(my_methods):
        lca.switch_method(method)
        for j, activity in enumerate(activities):
            lca.redo_lcia({activity: 1})
            # Add total to check afterwards that all exchanges add up to total
            result_array[i][j].append(("total", lca.score))
            for exc in activity.exchanges():
                if exc['type'] == 'technosphere':
                    lca.redo_lcia({exc.input: exc['amount']})
                    result_array[i][j].append((exc, lca.score))
                elif exc['type'] == 'biosphere':
                    # Need to multiple the amount times its CF
                    cf = lca.characterization_matrix[lca.biosphere_dict[exc.input], :].sum()
                    result_array[i][j].append((exc, cf * exc['amount']))
                    
    for method in result_array:
        for activity in method:
            if not np.allclose(activity[0][1], sum([o[1] for o in activity[1:]])):
                print("Mismatch")
                break

    grouped_array = [[group_exchange_scores(j) for j in i] for i in result_array]
    data_frames = []
    # Unpack in a dataframe:
    for i, group_data in enumerate(grouped_array):
        data_0 = dict(group_data[0])
        col_name = "multi_energy_system_{}_{}_eps:{}_credit:{}".format(loc_elect, w_cost, epsilon_constraint, credit_env_export)

        df_add = pd.DataFrame.from_dict(data_0, orient='index')
        df_add.rename(columns={"climate change total": col_name, 0: col_name}, inplace=True)
        df_add.index.names = ['contributor']
        df_add['category'] = str(my_methods[i][1])
        df_add['year'] = ASSESSMENT_YEAR
        df_add['db_name'] = sec_db
        data_frames.append(df_add)

    df_total = pd.concat(data_frames, axis=0)

    lca_results = df_total.reset_index().set_index(['category','contributor', 'year', 'db_name'])

    return lca_results

def get_tech_environmental_burdens(cost_dict, ei_loc, sec_db = NAME_REF_DB, lcia_method=CC_METHOD):
    """
    Gets the environmental impact from each activity used in the MES.

    Args:
        cost_dict (dic): dictionary with cost data [-].
        ei_loc (str): ecoinvent location [-].
        sec_db (str): ecoinvent database used [-].
        lcia_method (str): Standard LCIA method used, here CC_METHOD
        
    Returns:
        dict_env_impacts (dict): dictionary with environmental impacts for the 'lcia_method' specified.
        env_impact (float): environmental burden factor for grid absorption from the grid.
        env_impact (float): environmental burden credit for grid injection to the grid.
    """
    # Use BW project created in create_db_lca_functions
    # define name of project
    bw.projects.set_current(PROJECT_NAME) #Creating/accessing the project

    # Hydrogen storage vessel
    env_imp_h2_ves = get_activity_env("high pressure hydrogen storage tank", "GLO", 
                                      "high pressure hydrogen storage tank", sec_db, "", "", lcia_method=lcia_method) / (MJ_KG_H2/3.6) # per vessel, 1 kg H2 storage, convert to kWh

    # Ground-mounted solar PV panels     
    env_imp_pv = get_activity_env("photovoltaic open ground installation, 570 kWp, multi-Si, on open ground",
                                  "RER", "photovoltaic open ground installation, 570 kWp, multi-Si, on open ground", sec_db, "", "", lcia_method=lcia_method)/ (570*0.895) #per 570 kWp, however, consider degradation

    # Wind turbine:
    env_imp_wind_off = (get_activity_env("market for wind power plant, 2MW, offshore, fixed parts", "GLO", 
                                      "wind power plant, 2MW, offshore, fixed parts", sec_db, "", "", lcia_method=lcia_method) +
                       get_activity_env("market for wind power plant, 2MW, offshore, moving parts", "GLO", 
                                      "wind power plant, 2MW, offshore, moving parts", sec_db, "", "", lcia_method=lcia_method)) / 2000

    # onhore wind
    env_imp_wind_on = (get_activity_env("market for wind turbine, 2MW, onshore", "GLO", 
                                      "wind turbine, 2MW, onshore", sec_db, "", "", lcia_method=lcia_method) + get_activity_env("market for wind turbine network connection, 2MW, onshore", "GLO", 
                                      "wind turbine network connection, 2MW, onshore", sec_db, "", "", lcia_method=lcia_method)) / 2000
    # electrolyzer
    env_imp_electr = ( get_activity_env("electrolyzer production, 1MWe, PEM, Stack", "RER", 
                                      "electrolyzer, 1MWe, PEM, Stack", 
                                       sec_db, "", "", lcia_method=lcia_method) + (cost_dict['electr_lt']/cost_dict['electr_bos_lt']) * get_activity_env("electrolyzer production, 1MWe, PEM, Balance of Plant", "RER", 
                                      "electrolyzer, 1MWe, PEM, Balance of Plant", 
                                       sec_db, "", "", lcia_method=lcia_method) ) / 1000

    # battery, here NMC
    env_imp_bat_cap = get_activity_env("li-ion (NMC)", "GLO", "li-ion (NMC)", sec_db, "", "", lcia_method=lcia_method) + get_activity_env(
        "battery management system, kWh", "GLO", "battery management system, kWh", sec_db, "","", lcia_method=lcia_method) + get_activity_env(
        "energy management system, kWh", "GLO", "energy management system, kWh", sec_db, "", "", lcia_method=lcia_method)        

    env_imp_bat_p = get_activity_env("power conditioning system, container system", "GLO",
                                      "power conditioning system, container system", sec_db, "", "", lcia_method=lcia_method)

    # PEM Fuel cell
    env_imp_fc = get_activity_env("fuel cell system assembly, 1 kWe, proton exchange membrane (PEM)", "GLO", 
                                     "fuel cell system, 1 kWe, proton exchange membrane (PEM)", 
                                     sec_db, "", "", lcia_method=lcia_method) / 1

    # Residential heat pump
    env_imp_hp = get_activity_env("heat pump, air-water, 7kW", "RER", 
                                      "heat pump, air-water, 7kW", 
                                       sec_db, "", "", lcia_method=lcia_method) / (7 * AV_COP_ASHP)

    # CHP unit, assume micro gas turbine for now
    env_imp_chp = cost_dict['chp_mix_eff_e'] * get_activity_env("advanced gas turbine, 400kWe", "FI", 
                                      "advanced gas turbine, 400kWe", 
                                       sec_db, "", "", lcia_method=lcia_method) / 400
    env_imp_chp_mix = env_imp_chp

    # CHP unit, assume micro gas turbine for now
    env_imp_chp_fuel = get_activity_env("mixed CHP unit, operation", "RER", 
                                      "mixed CHP unit, operation", 
                                       sec_db, "", "", lcia_method=lcia_method)

    # gas boiler, assume oil boiler for now
    env_imp_gb = get_activity_env("gas boiler production", "RoW", 
                                      "gas boiler", 
                                       sec_db, "", "", lcia_method=lcia_method) / 10

    # natural gas impact per kWh heat
    env_imp_ng = get_activity_env('heat production, natural gas, w\o infrastructure', "Europe without Switzerland", 
                                      'heat production, natural gas, w\o infrastructure', 
                                       sec_db, "", "", lcia_method=lcia_method) * 3.6

    # Syngas
    syngas_co2_fu = get_activity_env('synthetic gas production, from wood, at fixed bed gasifier, w\o infrastructure', ei_loc, 
                                      'synthetic gas', sec_db, "", "", lcia_method=lcia_method) #m3
    energy_m3_syngas = 5.2 / 3.6 # MJ/m3 
    env_imp_syngas = syngas_co2_fu /energy_m3_syngas #kWh syngas

    # Gasifier
    env_imp_gasif = get_activity_env("synthetic gas factory construction", "RoW", 
                                     "synthetic gas factory", 
                                     sec_db, "", "", lcia_method=lcia_method) / 7400

    kwh_water = (4.2 * (cost_dict['heat_stor_t_max'] - cost_dict['heat_stor_t_min']) * 2000)/3600 #https://www.engineeringtoolbox.com/energy-storage-water-d_1463.html
    env_imp_heat_storage = get_activity_env("heat storage production, 2000l",  "RoW", "heat storage, 2000l",
                                      sec_db, "", "", lcia_method=lcia_method) / kwh_water

    env_imp_solar_heat_ha = (10000/1000) * get_activity_env("solar collector system installation, Cu flat plate collector, one-family house, combined system, w\o storage",  
                                             "CH", "solar collector system, Cu flat plate collector, one-family house, combined system",
                                      sec_db, "", "", lcia_method=lcia_method)/ 12.3 #m2

    env_impact_grid_network = get_activity_env("wind turbine network connection construction, 4.5MW, onshore", "GLO", 
                                              "wind turbine network connection, 4.5MW, onshore", sec_db, "", "", lcia_method=lcia_method) / 4500

    env_gas_car_km = get_activity_env("market for transport, passenger car, medium size, petrol, EURO 5", "GLO", 
                                              "transport, passenger car, medium size, petrol, EURO 5", 
                                               sec_db, "", "", lcia_method=lcia_method)

    env_bev_car_km = get_activity_env("transport, passenger car, electric, w\o fuel", "GLO", 
                                              "transport, passenger car, electric, w\o fuel", 
                                               sec_db, "", "", lcia_method=lcia_method)

    env_imp_smr = get_activity_env("hydrogen production, steam reforming", "RER", 
                                          "hydrogen, gaseous", 
                                           sec_db, "", "", lcia_method=lcia_method) / (MJ_KG_H2 / 3.6) #per kWh

    env_imp_e_heating = get_activity_env("electric storage heater, 5kW", "CH", 
                                          "electric storage heater, 5kW", 
                                           sec_db, "", "", lcia_method=lcia_method) / 5 #per kW

    env_imp_oil_heating = get_activity_env("heat production, light fuel oil, at boiler 10kW condensing, non-modulating", "Europe without Switzerland", 
                                          "heat, central or small-scale, other than natural gas", 
                                           sec_db, "", "", lcia_method=lcia_method) * 3.6 #per kWh
    
    env_imp_diesel_ht = get_activity_env("heat production, light fuel oil, at industrial furnace 1MW", "Europe without Switzerland", 
                                          "heat, district or industrial, other than natural gas", 
                                           sec_db, "", "", lcia_method=lcia_method) * 3.6 #per kWh

    env_impact, loc_spec = get_activity_env_elect("market for electricity, low voltage, Crete", ei_loc, 
                                      "electricity, low voltage", sec_db, "", "", lcia_method=lcia_method)

    env_impact_cons, loc_spec_cons = get_activity_env_elect("market for electricity, low voltage", ei_loc, 
                                  "electricity, low voltage", DB_NAME_CONS, "", "", lcia_method=lcia_method)
 
    ############################################## 
    ############################################## from kg CO2 to tonne is -1e3, from kWh to Mwh is factor 1e3
    ##############################################

    dict_env_impacts = {"ghg_imp_h2_ves":env_imp_h2_ves, #t/MWh
                    "ghg_imp_pv":env_imp_pv, #t/MWp
                    "ghg_imp_wind_off":env_imp_wind_off, #t/MWp
                    "ghg_imp_wind_on":env_imp_wind_on, #t/MWp 
                    "ghg_imp_electr": env_imp_electr, #t/MW
                    "ghg_imp_bat_cap": env_imp_bat_cap, #t/MWp
                    "ghg_imp_bat_p":env_imp_bat_p, #t/MW
                    "ghg_imp_fc": env_imp_fc, #t/MW
                    "ghg_imp_hp":env_imp_hp, #t/MWth
                    "ghg_imp_chp":env_imp_chp, #/MWe
                    "ghg_imp_chp_mix":env_imp_chp_mix,#t/MWe
                    "ghg_imp_gb":env_imp_gb,#t/MWth
                    "ghg_imp_ng":env_imp_ng,#t/MWth
                    "ghg_imp_syngas":env_imp_syngas,#t/MWh
                    "ghg_imp_gasif":env_imp_gasif,#t/MWth
                    "ghg_imp_heat_storage":env_imp_heat_storage,#t/MWh
                    "ghg_imp_solar_heat_ha":env_imp_solar_heat_ha,#t/ha
                    "ghg_impact_grid_network":env_impact_grid_network,#t/MW
                    "ghg_gas_car_km":env_gas_car_km/1e3,#t/km
                    "ghg_imp_smr": env_imp_smr, # t/MWh H2
                    "ghg_imp_chp_fuel":env_imp_chp_fuel,#t/MWh
                    "ghg_bev_car_km":env_bev_car_km/1e3,#t/km
                    "ghg_imp_e_heating":env_imp_e_heating, #t/MW
                    "ghg_imp_oil_heating":env_imp_oil_heating,#t/MWh
                    "ghg_imp_diesel_ht":env_imp_diesel_ht,#t/MWh
                    "ghg_imp_grid_power":env_impact,#t/MWh
                    "ghg_imp_grid_power_cons":env_impact_cons,#t/MWh
                   }
        
    if lcia_method != my_methods[0]:
         dict_env_impacts.keys.replace("ghg_imp","env_imp").replace("ghg","env")
    
    return dict_env_impacts


#ecoinvent_remind_SSP2-PkBudg1300_2030_all
def get_activity_env(name: str, location: str, ref_product: str, 
                     db: str, year: str, scenario: str, lcia_method = 
                     CC_METHOD) -> float:   
    """
    Gets the environmental impact of an activity from a specified ecoinvent database.

    Args:
        name (str): activity name [-].
        location (str): activity location [-].
        ref_product (str): reference product [-].
        year (str): year of database [-].
        scenario (str): IAM scenario used [-].
        lcia_method (str): Standard LCIA method used, here CC_METHOD
        
    Returns:
        float: environmental impact.
        string: location of activity found.
    """
    
    if db == "":
        db_name = "ecoinvent_remind_{}_{}_all".format(scenario, year)
    else:
        db_name = db
    
    # For PV db, we don't have a reference product
    if ref_product == "":
        activity = [x for x in bw.Database(db_name) if name == x['name'] and
               location == x['location'] ][0]
    else:
        activity = [x for x in bw.Database(db_name) if name == x['name'] and
               location == x['location'] and ref_product == x['reference product']
                   ][0]
    lca = bw.LCA({activity: 1}, method=lcia_method)
    lca.lci()
    lca.lcia()
    
    return lca.score

#ecoinvent_remind_SSP2-PkBudg1300_2030_all
def get_activity_env_elect(name: str, location: str, ref_product: str, 
                     db: str, year: str, scenario: str, lcia_method = 
                     CC_METHOD) -> float:   
    """
    Gets the environmental impact of an electricity activity from a specified ecoinvent database.

    Args:
        name (str): activity name [-].
        location (str): activity location [-].
        ref_product (str): reference product [-].
        year (str): year of database [-].
        scenario (str): IAM scenario used [-].
        lcia_method (str): Standard LCIA method used, here CC_METHOD
        
    Returns:
        float: environmental impact.
    """
    
    if db == "":
        db_name = "ecoinvent_remind_{}_{}_all".format(scenario, year)
    else:
        db_name = db
        
    activity = [x for x in bw.Database(db) if name == x['name'] and
           location == x['location'] and
           ref_product == x['reference product']
               ]
    
    if len(activity) < 1:
        # Check whether there is a market group activity for larger area
        activity = [x for x in bw.Database(db) if x['name'] == "market group for electricity, low voltage" and
               location == x['location'] and
               ref_product == x['reference product']
                   ]
        
        if len(activity) < 1:
            # Try to select GLO activity
            activity = [x for x in bw.Database(db) if x['name'] == "market group for electricity, low voltage" and
                   "GLO" == x['location'] and
                   ref_product == x['reference product']
                       ]

            if len(activity) < 1:
                # Try to select RoW activity
                activity = [x for x in bw.Database(db) if name == x['name'] and
                       "RoW" == x['location'] and
                       ref_product == x['reference product']]
                if len(activity) == 1:
                    # Select this activity
                    activity = activity[0]
                else:
                    print("ERROR: No activity found for '{}'".format(name))                
            else:
                # Select this activity
                activity = activity[0]                      

        elif len(activity) == 1:
            # Select this activity
            activity = activity[0]
        else:
            print("ERROR: More than 1 activity found for '{}'".format(name))
            
    elif len(activity) == 1:
        # Select this activity
        activity = activity[0]
    else:
        print("ERROR: More than 1 activity found for '{}'".format(name))
    
    lca = bw.LCA({activity: 1}, method=lcia_method)
    lca.lci()
    lca.lcia()
    
    return lca.score, activity['location']

def calc_crf(r: float, lt: float) -> float:
    """
    Calculates capital recovery factor.

    Args:
        r (float): dicount/interest rate/WACC [-].
        lt (float): lifetime of loan [years].
    Returns:
        float: capital recovery factor.
    """
    return (r * (1 + r) ** lt) / ((1 + r) ** lt - 1)

def rep_annual_int(unit_inv_cost: float, lt: int, r: float, ry: float, rep_factor=0.75) -> float:
    """
    Calculates replacement costs for a technology.

    Args:
        unit_inv_cost (float): Replacement cost [CHF/euro/...].
        lt (int): Lifetime of the project or system [years].
        r (float): Interest rate [%].
        ry (float): Replacement year [year].
        rep_factor (float): share of replacement costs compared to full investment [-].

    Returns:
        float: Total replacement cost.
    """
    c_rep_total = 0
    
    # 75% of initial costs for replacement, assuming some components can be used and due to technological improvement
    unit_inv_cost = rep_factor * unit_inv_cost
    
    # 1. If the lifetime of the system is bigger than the replacement year, there is a residual value
    if ry > lt:
        # Check the amount left over, i.e. how much of the component is left after the system lifetime
        rest = 1 - (lt/ry)
        crf = calc_crf(r, lt)
        c_rep_a = crf*((unit_inv_cost * rest) / ((1+r)**lt)) 
        # Reduce the replacement costs
        c_rep_total -= c_rep_a
        
    # 2. If the replacement year of a component is smaller than the system lifetime, then there are replacement costs
    elif ry < lt:
        replacement_times = (lt/ry)
        counter_ry = ry
        crf = calc_crf(r, lt)
        
        # Iterate through replacement times
        for _ in range(int(replacement_times)):
            # Calculate replacement costs
            c_rep_a = crf*(unit_inv_cost / ((1+r)**ry))
            
            # Add to the total replacement costs
            c_rep_total += c_rep_a
            
            # Update counters
            ry += counter_ry
                
        # Calculate residual value
        rest = 1 - (replacement_times % 1)
        c_rep_a = crf*((unit_inv_cost * rest) / ((1+r)**lt))
                
        # Reduce the replacement costs
        c_rep_total -= c_rep_a
    
    # 3. Lifetime of a system component is exactly the lifetime of the system   
    else:
        c_rep_total = 0
    
    return c_rep_total

def check_simultaneous_operations(df_out, error_message):
    """
    Checks for simultaneous operations and raises an error if found.

    Args:
        df_out (DataFrame): DataFrame containing the relevant columns for the checks (from the optimization).
        error_message (str): The error message to be raised if a simultaneous operation is found.
    """
    checks = [
        ('Battery', df_out['p_battch'], df_out['p_battdis']),
        ('Heat Storage Medium', df_out['heat_stor_ch'], df_out['heat_stor_dis']),
        ('Grid', df_out['p_grid_abs'], df_out['p_grid_inj'])
    ]

    for label, charge, discharge in checks:
        simultaneous = charge[(charge > 10**-4) & (discharge > 10**-4)]
        if simultaneous.sum() > 0:
            df_out.to_excel("error_file.xlsx")
            raise ValueError(f"WARNING: {label} is charging and discharging at the same time")

def create_process(location, year, exchanges):
    """
    Create a process for a multi-energy system.

    Parameters:
        location (str): The location of the multi-energy system.
        year (int): The year of the multi-energy system.
        exchanges (list): A list of exchanges associated with the process.

    Returns:
        dict: A dictionary representing the LCA process of the multi-energy system process.
    """
    
    name = "multi_energy_system_{}_{}".format(location, year)  
    return {
        'name': name,
         "code": str(uuid.uuid4().hex),
        'unit': 'unit',
        'reference product': name,
        'location' :location,
        'exchanges': exchanges, 
    }

def group_exchange_scores(lst):
    """
    Group exchange scores based on a list of exchanges and their scores.

    Parameters:
        lst (list): A list of tuples where each tuple contains an exchange and its associated score.

    Returns:
        dict: A dictionary with group labels as keys and corresponding LCIA scores as values.
    """
    
    # We will store results in a dict, with keys for group label and values of LCIA score.
    grouped_results = defaultdict(int)

    for exc, score in lst[1:]:
        grouped_results[contribution_mapping_system[exc.input['name']]] += score
    return grouped_results

def drop_empty_categories(db):
    """
    Drop categories with the value ('',) from the database.

    Parameters:
        db (list): A list of processes or datasets in a database.

    Returns:
        list: The modified database with empty categories removed.
    """

    DROP = ('',)
    for ds in db:
        if ds.get('categories') == DROP:
            del ds['categories']
        for exc in ds.get("exchanges", []):
            if exc.get('categories') == DROP:
                del exc['categories']
    return db

def strip_nonsense(db):
    """
    Strip leading and trailing spaces from strings in the database.

    Parameters:
        db (list): A list of processes or datasets in a database.

    Returns:
        list: The modified database with leading and trailing spaces removed from string values.
    """
    for ds in db:
        for key, value in ds.items():
            if isinstance(value, str):
                ds[key] = value.strip()
            for exc in ds.get('exchanges', []):
                for key, value in exc.items():
                    if isinstance(value, str):
                        exc[key] = value.strip()
    return db

def create_process_and_add(location, year, exchanges, sec_db):
    """Creates database with processes"""
    
    IMPORTER.data = [create_process(location, year, exchanges)]
    
    IMPORTER.strategies = [
        partial(add_database_name, name=DB_NAME),
        csv_restore_tuples,
        drop_empty_categories,
        strip_nonsense,
    ]
    IMPORTER.apply_strategies()
    IMPORTER.match_database(sec_db, fields=('name','unit','location','reference product', 'database'))
    IMPORTER.match_database(DB_NAME_CONS, fields=('name','unit','location','reference product', 'database'))   
    IMPORTER.match_database(fields = ('name',))
    IMPORTER.statistics()
    IMPORTER.write_excel(only_unlinked=True)
    IMPORTER.write_database()

def parameter_increase_sensitivity(boolean_ind, dict_ghg_impacts, param_change='grid_capex', mip_gap=0.001, 
        max_nr=30, step_inc=400, export_results=True):
    """
    Perform sensitivity analysis for a certain parameter increase, such as the grid network.

    Parameters:
    - boolean_ind (bool): Flag to indicate whether to calculate sensitivity.
    - dict_ghg_impacts (dict): Dictionary with emissio data.
    - param_change (str): Parameter to increase per iterative step. Default is 'grid_capex'.
    - mip_gap (float): MIP gap for optimization.
    - max_nr (int): Maximum number of sensitivity iterations.
    - step_inc (int): Step increase of the parameter per iteration.
    - export_results (bool): Flag to indicate whether to export results to Excel.

    Returns:
    - df_sens_op_param_gr (pd.DataFrame): Results of sensitivity analysis for the specified parameter.
    - df_sens_op_param_env_gr (pd.DataFrame): Environmental results of sensitivity analysis.
    """
    if boolean_ind:
        # Initialize and process data
        data_processor = ep.EnergyDataProcessor(
            ASSESS_COUNTRY, ASSESSMENT_YEAR, LATITUDE, LONGITUDE,
            fit=FIT, constrained=False, residential_sector=True
        )

        df_data, dict_limits = data_processor.process_data()
        df_data['ghg_impact'], df_data['ghg_impact_cons'] = dict_ghg_impacts["ghg_imp_grid_power"], dict_ghg_impacts["ghg_imp_grid_power_cons"]
        cost_dict = COST_DATA[NAME_REF_DB].to_dict()

        # Initialize results DataFrames outside the loop
        df_init_sens, df_init_sens_env = None, None

        # Loop for sensitivity analysis
        for i in range(max_nr):
            print(f"\nStatus [{i+1}/{max_nr}]\n")

            # Modify cost dict
            cost_dict[param_change] = step_inc * i

            total_result, lca_result, total_fig = optimize_mes(
                df_data, 1, 0, cost_dict, dict_ghg_impacts,
                dict_limits, ghg_init=1, cost_init=1,
                export_results=False, mip_gap=mip_gap,
            )
            total_result[param_change] = cost_dict[param_change]

            # Store results in DataFrames
            df_init_sens = pd.concat([df_init_sens, total_result], axis=0, ignore_index=True)
            df_init_sens_env = pd.concat([df_init_sens_env, lca_result], axis=1)

        # Calculate share of power grid inv compared to overall annual cost
        df_init_sens['share_grid_total_cost'] = 100 * df_init_sens.an_costs_capex_grid / df_init_sens.cost_init
        # Calculate share of power grid inv compared to overall annual cost
        df_init_sens['share_grid_inv'] = 100 * df_init_sens.an_costs_capex_grid / df_init_sens.filter(like='capex').T.sum()

        if export_results:
            df_init_sens.T.to_excel(f"opt_results\\df_sens_op_{param_change}_gr.xlsx")
            df_init_sens_env.to_excel(f"opt_results\\df_sens_op_{param_change}_env_gr.xlsx")

    else:
        df_init_sens = pd.read_excel(f"opt_results\\df_sens_op_{param_change}_gr.xlsx", index_col=[0]).T
        df_init_sens_env = pd.read_excel(f"opt_results\\df_sens_op_{param_change}_env_gr.xlsx", index_col=[0, 1, 2, 3])

    df_sens_op_param_gr = df_init_sens.set_index(param_change)
    return df_sens_op_param_gr, df_init_sens_env

def local_sens_analysis(dict_ghg_impacts, sens_factor, 
                        factor_change, cost_init, ghg_init, constrained=False, residential_sector=True,
                        autonomous_elect=True, autonomous_gas=True):
    """
    Perform local sensitivity analysis on off-grid systems.

    Parameters:
    - dict_ghg_impacts (dict): Dictionary containing environmental impacts for different technologies.
    - ghg_impact (float): Sensitivity factor for environmental impacts.
    - ghg_impact_cons (float): Sensitivity factor for environmental impacts of consumed electricity.
    - sens_factor (str): The factor to be sensitized (e.g., 'capex', 'wood_price', 'lt', 'om', 'irradiance', 'temp_air').
    - factor_change (float): Factor for changing the sensitized parameter. Default is -0.1 for a 10% decrease.
    - cost_init (float): Initial total cost for comparison.
    - ghg_init (float): Initial total GHG emissions for comparison.
    - autonomous_elect (bool, optional): Optimize as an autonomous energy system without a power grid connection.
            Default is True.
    - autonomous_gas (bool, optional): Optimize as an autonomous energy system without a gas grid connection.
            Default is False.

    Returns:
    - total_result (pd.DataFrame): DataFrame containing the results of the multi-objective optimization.
    - lca_result (pd.DataFrame): DataFrame containing the environmental impacts results.

    Note:
    - The input data paths, functions like 'data_processor.get_renewable_profiles', 'get_energy_demands_region', and other 
      external functions need to be defined or imported before calling this function.
    """

    # Increase factor
    factor = 1 + factor_change
    
    # Load data from Excel file for time preferences of BEVs.
    TPs = pd.read_excel(r"new_data\TP_BEVs.xlsx", index_col=[0])
    TP_1 = list(TPs.TP1) * 365
    TP_2 = list(TPs.TP2) * 365
    TP_3 = list(TPs.TP3) * 365
    
    data_processor = ep.EnergyDataProcessor(ASSESS_COUNTRY, ASSESSMENT_YEAR, LATITUDE, LONGITUDE, \
        fit=FIT, constrained=constrained, residential_sector=residential_sector)

    # ### 10.2.1 Generate some demand and generation profiles
    weather_data, el = data_processor.get_weather_data()

    # Check if the sensitivity factor is in the weather_data columns
    if sens_factor in list(weather_data.columns):
        print(weather_data[sens_factor].mean())
        if sens_factor == 'temp_air':
            weather_data.loc[weather_data[sens_factor] > 0] *= factor
            weather_data.loc[weather_data[sens_factor] < 0] *= (2 - factor)
        else:
            weather_data[sens_factor] *= factor
            
        print(weather_data[sens_factor].mean())
    
    # Get renewable electricity generation profiles
    pv_MW_array, wind_MW_array_off, wind_MW_array_on = data_processor.get_renewable_profiles(weather_data)
    
    # Copy the initial cost dictionary
    cost_dict = COST_DATA[NAME_REF_DB].to_dict()
    
    # Modify cost dictionary based on the sensitivity factor dynamically
    for k in cost_dict:
        if sens_factor == 'capex': 
            if 'capex' in str(k):
                cost_dict[k] *= factor 
                print("Dict changes {} {}".format(cost_dict[k], k))
        elif sens_factor == 'wood_price': 
            if (str(k) == 'wood_price') or (str(k) == 'wood_price_res'):
                cost_dict[k] *= factor
                print("Dict changes {} {}".format(cost_dict[k], k))
        elif sens_factor == "lt":                 
            if ('_lt' in str(k)) and (k != "project_lt"):
                cost_dict[k] *= factor     
                print("Dict changes {} {}".format(cost_dict[k], k))
        elif sens_factor == "om":                 
            if ('_om' in str(k)):
                cost_dict[k] *= factor
                print("Dict changes {} {}".format(cost_dict[k], k))
        else:
            if sens_factor == k: 
                cost_dict[k] *= factor     
                print("Dict changes {} {}".format(cost_dict[k], k))

    # Get energy demands for different sectors
    cop_array, elect_demand_ind, ind_heat_demand, bld_heat_demand, elect_demand_res = ep.get_energy_demands_region(weather_data.temp_air)

    # Generate wind profiles
    wind_MW_start_on_micro = ep.wind_profile_generator_tmy(weather_data, LATITUDE, LONGITUDE, turbine_spec="E48/800", hub_height=40)

    # Create a DataFrame with input data
    df_data = pd.DataFrame(data={'pv_MW_array': pv_MW_array, 'wind_MW_array_off': wind_MW_array_off,
                                 'wind_MW_array_on': wind_MW_array_on, 'wind_MW_array_on_micro': wind_MW_start_on_micro,
                                 "temperature": weather_data.temp_air, "irradiance": weather_data['ghi'] / 1000,
                                 "ind_heat_demand": np.array(ind_heat_demand), "bld_heat_demand": np.array(bld_heat_demand),
                                 "elect_demand_ind": np.array(elect_demand_ind),
                                 "cop_array": cop_array, "electricity_demand": np.array(
                                     elect_demand_res.values + elect_demand_ind.values),
                                 "elect_demand_res": np.array(elect_demand_res), "grid_abs_price": [0]*8760,
                                 "ghg_impact": dict_ghg_impacts["ghg_imp_grid_power"], "ghg_impact_cons": dict_ghg_impacts["ghg_imp_grid_power_cons"],
                                 "rev_inj": [0]*8760, "TP_1": TP_1, "TP_2": TP_2, "TP_3": TP_3,
                                 }, index=pd.date_range('1/1/{} 00:00'.format(ASSESSMENT_YEAR), periods=8760, freq='H'))
    
    # Check if the sensitivity factor is in the df_data columns and not in the weather_data columns
    if (sens_factor in list(df_data.columns)) and (sens_factor not in list(weather_data.columns)):
        print(df_data[sens_factor].mean())
        df_data[sens_factor] *= factor
        print(df_data[sens_factor].mean())

    # Get maximum capacity constraints for regions
    dict_limits = ep.get_max_caps_regions(constrained=constrained, residential_sector = residential_sector)

    # Perform multi-objective optimization
    total_result, lca_result, total_fig = optimize_mes(df_data, 1, 0, cost_dict, dict_ghg_impacts, 
                                                                dict_limits, export_results=False,
                                                                ghg_init=ghg_init,
                                                                cost_init=cost_init,
                                                                autonomous_elect=autonomous_elect, 
                                                                autonomous_gas=autonomous_gas,
                                                                mip_gap = 0.001 # smaller mip gap, as optimal installed capacities are quite sensitive to it.
                                                                )

    # Calculate and print cost and GHG changes
    cost_change = (total_result.total_costs.item() - cost_init) / cost_init
    total_result['Costs [$\Delta$%]'] = cost_change * 100

    ghg_change = (total_result.total_ghg.item() - ghg_init) / ghg_init
    total_result['GHG emissions [$\Delta$%]'] = ghg_change * 100

    print("Cost change is '{}'% and GHG change is '{}'%".format(round(cost_change * 100, 2), round(ghg_change * 100, 2)))
    
    return total_result, lca_result
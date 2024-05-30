import pandas as pd 
import matplotlib.pyplot as plt
import json

# import own Python files, vars, mappings, and functions
from config import (CC_METHOD, NAME_REF_DB, CALC_ALL_LCA_IMPACTS, ASSESSMENT_YEAR, DB_NAME_CONS, DB_NAME_INIT, LOCATION_DB,
                    LATITUDE, LONGITUDE, FIT, LOC_ELECT, ASSESS_COUNTRY, COST_DATA, GENERATE_NEW_LCA_DB, LOCATION_DB_CONS,
                    SHARE_E_RES, HP_SHARE, TOTAL_HH)

from mapping import (my_methods, list_prods_e, list_cons_e, list_prods_lt, list_cons_lt, list_prods_ht, list_cons_ht,
                    list_prods_h2, list_cons_h2, storage_h2, full_names_mapping_dict, storage_e)

import plotting_functions as pf
import energy_data_processor as ep
import opt_lca_functions as opt
import create_db_lca_functions as lcaf

# # SETTINGS
# All sector scenarios
COST_MIN_OPT = False
GHG_MIN_OPT = False
AUTO_ALL_MIN_OPT = False
COST_MIN_OPT_CONSTR= False
AUTO_BAL_MIN_OPT = False

# Bakery industry
COST_MIN_OPT_BI = False
COST_OPT_NON_CONSTR_BI = False

# Sensitivity
CALC_SENS_OP_GRID_CAPEX = False
CALC_SENS = False
CALC_SENS_EXTREMES = False

# Pareto
GHG_MIN_80 = False
CALC_PARETO = False
GEN_FIGS = True

if __name__ == "__main__":
    # ### Optimization model parameters
    cost_dict = COST_DATA[NAME_REF_DB].to_dict() # techno-economic data

    # If we want to use a full LCA approach, we have to set-up the LCA database:
    if GENERATE_NEW_LCA_DB:
        # Import LCIA methods, import ecoinvent cut-off and consequential dbs
        lcaf.import_additional_lcias()
        lcaf.import_ecoinvent_database(DB_NAME_INIT, LOCATION_DB)
        lcaf.import_ecoinvent_database(DB_NAME_CONS, LOCATION_DB_CONS)
        # Generate reference database with premise to import additional novel LCIs.
        lcaf.generate_reference_database()
        # Modify some datasets, needed for the location-specific MES.
        lcaf.generate_specialized_datasets()

    # Generate GHG emission data
    if CALC_ALL_LCA_IMPACTS:
        dict_ghg_impacts = opt.get_tech_environmental_burdens(cost_dict, ASSESS_COUNTRY) 
        with open("new_data\dict_ghg_impacts.txt", 'w') as file:
            json.dump(dict_ghg_impacts, file)
    else:
        #Otherwise, just used the stored data valid for ecoinvent 3.9.1 cut-off
        with open("new_data\dict_ghg_impacts.txt", 'r') as file:
            dict_ghg_impacts = json.load(file)

    # Initialize and process data
    data_processor = ep.EnergyDataProcessor(ASSESS_COUNTRY, ASSESSMENT_YEAR, LATITUDE, LONGITUDE, fit=FIT, constrained=False, residential_sector=True)
    df_data, dict_limits = data_processor.process_data()
    df_data['ghg_impact'],  df_data['ghg_impact_cons'] = dict_ghg_impacts["ghg_imp_grid_power"], dict_ghg_impacts["ghg_imp_grid_power_cons"]

    for col in df_data.columns:
        if 'demand' in str(col):
            print("Summed '{}' is '{}' GWh".format(col,round(df_data[col].sum()/1000,1)))

    # Initialize and process data
    data_processor_bi = ep.EnergyDataProcessor(ASSESS_COUNTRY, ASSESSMENT_YEAR, LATITUDE, LONGITUDE, fit=FIT, constrained=True, residential_sector=False)
    df_data_bi, dict_limits_bi = data_processor_bi.process_data()
    df_data_bi['ghg_impact'],  df_data_bi['ghg_impact_cons'] = dict_ghg_impacts["ghg_imp_grid_power"], dict_ghg_impacts["ghg_imp_grid_power_cons"]

    # Generate some demand and generation profiles
    weather_data, el = data_processor.get_weather_data()

    plot_info = {'temp_air': {'name': 'Temperature', 'y_name': 'Temperature [$^\circ$C]'},
                'ghi': {'name': 'Irradiance (GHI)', 'y_name': 'Irradiance [W/m$^2$]'},
                'wind_speed': {'name': 'Windspeed', 'y_name': 'Windspeed [m/s]'} }

    # Plot profiles for temperature, irradiance, and windspeed.
    pf.plot_weather_profiles(plot_info,weather_data)

    # Dictionary for mapping list_df_data values to names and y_names
    data_mapping = {
        #'ind_heat_demand': {'name': 'Industrial heat demand', 'y_name': 'Heat demand [MW]'}, # confidential
        'bld_heat_demand': {'name': 'Residential heat demand', 'y_name': 'Heat demand [MW]'},
        'electricity_demand': {'name': 'Electricity demand', 'y_name': 'Electricity demand [MW]'},
        'grid_abs_price': {'name': 'Grid electricity price', 'y_name': 'Price [k€/MWh]'},
        'ghg_impact': {'name': 'GHG intensity grid', 'y_name': 'GHG intensity [tCO$_2$-eq./MWh]'},
        'pv_MW_array': {'name': 'PV generation [MW$_p$]', 'y_name': 'Generation [MW/MW$_p$]'},
        'wind_MW_array_off': {'name': 'Offshore generation [MW$_p$]', 'y_name': 'Generation [MW/MW$_p$]'},
        'wind_MW_array_on': {'name': 'Onshore generation [MW$_p$]', 'y_name': 'Generation [MW/MW$_p$]'},
        'wind_MW_array_on_micro': {'name': 'Onshore generation (micro turbine) [MW$_p$]', 'y_name': 'Generation [MW/MW$_p$]'} }

    # Plot profiles.
    pf.plot_g_d_profiles(data_mapping, df_data)
    
    # Initialize bakery_industry_system
    bakery_industry_system = ep.ReferenceEnergySystemCapacities(df_data_bi, cost_dict, dict_limits_bi, total_hh = 0,
                                                         share_res_oil_boiler=0, share_res_hp=0, share_res_e_heating=0, 
                                                         share_res_ng_boiler=0, share_ind_boiler=1)
    
    bakery_industry_system.determine_capacities()
    bakery_industry_system.calculate_residential_heating_capacities()
    bakery_industry_system.calculate_grid_absorption()
    bakery_industry_system.calculate_industrial_boiler()
    bakery_industry_system.calculate_transportation()
    
    # calculate all impacts of the reference system
    overview_totals_zero_bi, lca_results_zero_bi = opt.calc_energy_system_results(df_data_bi, cost_dict, dict_ghg_impacts, LOC_ELECT, total_hh = 0,
                                cap_wind_off=0,
                                cap_wind_on=0,
                                cap_pv=0,
                                cap_bat_en=0,
                                cap_bat_p=0,
                                cap_h2_ves=0,
                                cap_electrolyzer=0,
                                cap_fc=0,
                                cap_hp=bakery_industry_system.cap_hp_res,
                                cap_d_boiler=bakery_industry_system.cap_d_boiler_ind,
                                cap_gb_res=bakery_industry_system.cap_boiler_ng_res,
                                cap_heat_storage=0,
                                cap_solar_heat=0,
                                cap_chp_mix=0,
                                cap_gasif=0,
                                p_grid_connection= bakery_industry_system.p_grid_connection,
                                cap_e_boiler= 0,
                                cap_e_heating = bakery_industry_system.cap_e_heating_res,
                                cap_boiler_oil_res = bakery_industry_system.cap_boiler_oil_res,
                                summed_grid_abs = bakery_industry_system.p_grid_abs,
                                summed_grid_inj = 0,
                                summed_syngas = 0,
                                summed_heat_boiler_ind = bakery_industry_system.p_d_boiler,
                                summed_heat_boiler_res = bakery_industry_system.p_boiler_gb_res,
                                summed_heat_oil_boiler_res = bakery_industry_system.p_boiler_oil_res,
                                summed_chp_fuel = 0,
                                summed_biogas= 0,
                                bev_car_kms = bakery_industry_system.bev_car_kms,
                                petrol_car_kms= bakery_industry_system.petrol_car_kms,
                                sec_db=NAME_REF_DB,
                                h2_export=0,
                                credit_env_export=False,
                                lcia_method=CC_METHOD)
    # Initialize entire_mes
    entire_mes = ep.ReferenceEnergySystemCapacities(df_data, cost_dict, dict_limits_bi, TOTAL_HH,
                                                 share_res_oil_boiler=(1-SHARE_E_RES-HP_SHARE), share_res_hp=HP_SHARE, share_res_e_heating=SHARE_E_RES, 
                                                 share_res_ng_boiler=0, share_ind_boiler=1)
    
    # Calculate capacities for entire_mes
    entire_mes.determine_capacities()
    entire_mes.calculate_residential_heating_capacities()
    entire_mes.calculate_grid_absorption()
    entire_mes.calculate_industrial_boiler()
    entire_mes.calculate_transportation()

    # Calculate impacts of the reference system for entire_mes
    overview_totals_zero, lca_results_zero = opt.calc_energy_system_results(df_data, cost_dict, dict_ghg_impacts, LOC_ELECT, total_hh = TOTAL_HH,
                                cap_wind_off=0,
                                cap_wind_on=0,
                                cap_pv=0,
                                cap_bat_en=0,
                                cap_bat_p=0,
                                cap_h2_ves=0,
                                cap_electrolyzer=0,
                                cap_fc=0,
                                cap_hp=entire_mes.cap_hp_res,
                                cap_d_boiler=entire_mes.cap_d_boiler_ind,
                                cap_gb_res=entire_mes.cap_boiler_ng_res,
                                cap_heat_storage=0,
                                cap_solar_heat=0,
                                cap_chp_mix=0,
                                cap_gasif=0,
                                p_grid_connection= entire_mes.p_grid_connection,
                                cap_e_boiler= 0,
                                cap_e_heating = entire_mes.cap_e_heating_res,
                                cap_boiler_oil_res = entire_mes.cap_boiler_oil_res,
                                summed_grid_abs = entire_mes.p_grid_abs,
                                summed_grid_inj = 0,
                                summed_syngas = 0,
                                summed_heat_boiler_ind = entire_mes.p_d_boiler,
                                summed_heat_boiler_res = entire_mes.p_boiler_gb_res,
                                summed_heat_oil_boiler_res = entire_mes.p_boiler_oil_res,
                                summed_chp_fuel = 0,
                                summed_biogas= 0,
                                bev_car_kms = entire_mes.bev_car_kms,
                                petrol_car_kms= entire_mes.petrol_car_kms,
                                sec_db=NAME_REF_DB,
                                h2_export=0,
                                credit_env_export=False,
                                lcia_method=CC_METHOD)   
    
    print("mes oil boiler {}, {}".format(entire_mes.cap_boiler_oil_res, entire_mes.p_boiler_gb_res))
    """
    Solve optimization scenarios.
    """
    # # 1.1 Scenario -- Cost Min (entire MES)
    ghg_init = overview_totals_zero.total_ghg.item()
    cost_init = overview_totals_zero.total_costs.item()

    if COST_MIN_OPT:
        totals_cost_min_gr, lca_cost_min_gr, output_cost_min_gr = opt.optimize_mes(df_data, 1, 0,
                                                            cost_dict, dict_ghg_impacts, dict_limits,
                                                            ghg_init = ghg_init,
                                                            cost_init = cost_init,
                                                            export_alias = 'cost_base'                           
                                                            )
        
        totals_cost_min_gr.T.to_excel(r"opt_results\totals_cost_min_gr.xlsx")
        lca_cost_min_gr.to_excel(r"opt_results\lca_cost_min_gr.xlsx")
    else:
        totals_cost_min_gr = pd.read_excel(r"opt_results\totals_cost_min_gr.xlsx", index_col=[0]).T
        lca_cost_min_gr = pd.read_excel(r"opt_results\lca_cost_min_gr.xlsx", index_col=[0,1,2,3])

    # # 1.2 Scenario -- Cost Min Constrained (entire MES)
    data_processor = ep.EnergyDataProcessor(ASSESS_COUNTRY, ASSESSMENT_YEAR, LATITUDE, LONGITUDE, fit=FIT, constrained=True, residential_sector=True)
    df_data_constr, dict_limits_constr = data_processor.process_data()

    if COST_MIN_OPT_CONSTR:
        totals_cost_min_gr_constr, lca_cost_min_gr_constr, output_cost_min_gr_constr = opt.optimize_mes(df_data, 1, 0,
                                                            cost_dict, dict_ghg_impacts, dict_limits_constr,
                                                            ghg_init = ghg_init,
                                                            cost_init = cost_init,
                                                            micro_wt = True
                                                            )
        
        totals_cost_min_gr_constr.T.to_excel(r"opt_results\totals_cost_min_gr_constr.xlsx")
        lca_cost_min_gr_constr.to_excel(r"opt_results\lca_cost_min_gr_constr.xlsx")
    else:
        totals_cost_min_gr_constr = pd.read_excel(r"opt_results\totals_cost_min_gr_constr.xlsx", index_col=[0]).T
        lca_cost_min_gr_constr = pd.read_excel(r"opt_results\lca_cost_min_gr_constr.xlsx", index_col=[0,1,2,3])

    # 1.3 Scenario -- GHG-Min (Entire MES)
    if GHG_MIN_OPT:
        totals_ghg_min_gr, lca_ghg_min_gr, output_ghg_min_gr = opt.optimize_mes(df_data, 0, 1,
                                                            cost_dict, dict_ghg_impacts, dict_limits,
                                                            ghg_init = ghg_init,
                                                            cost_init = cost_init, time_limit = 12*3600,
                                                            heuristics=True,
                                                            export_alias = 'ghg_base'
                                                            )
        
        totals_ghg_min_gr.T.to_excel(r"opt_results\totals_ghg_min_gr.xlsx")
        lca_ghg_min_gr.to_excel(r"opt_results\lca_ghg_min_gr.xlsx")
    else:
        totals_ghg_min_gr = pd.read_excel(r"opt_results\totals_ghg_min_gr.xlsx", index_col=[0]).T
        lca_ghg_min_gr = pd.read_excel(r"opt_results\lca_ghg_min_gr.xlsx", index_col=[0,1,2,3])

    # # 1.4 Cost-Min (Entire MES), Autonomous - from power grid & gas 
    if AUTO_ALL_MIN_OPT:
        totals_auto_all_min_gr, lca_auto_all_min_gr, output_auto_all_min_gr = opt.optimize_mes(df_data, 1, 0, 
                                                            cost_dict, dict_ghg_impacts, dict_limits, autonomous_elect=True, autonomous_gas=True,
                                                            ghg_init = ghg_init,
                                                            cost_init = cost_init
                                                            )
        
        totals_auto_all_min_gr.T.to_excel(r"opt_results\totals_auto_all_min_gr.xlsx")
        lca_auto_all_min_gr.to_excel(r"opt_results\lca_auto_all_min_gr.xlsx")
    else:
        totals_auto_all_min_gr = pd.read_excel(r"opt_results\totals_auto_all_min_gr.xlsx", index_col=[0]).T
        lca_auto_all_min_gr = pd.read_excel(r"opt_results\lca_auto_all_min_gr.xlsx", index_col=[0,1,2,3])

    # # 1.5 Cost-Min (Entire MES), Balanced autonomy, including autonomy from gas and diesel.
    if AUTO_BAL_MIN_OPT:
        totals_auto_bal_min_gr, lca_auto_bal_min_gr, output_auto_bal_min_gr = opt.optimize_mes(df_data, 1, 0, 
                                                            cost_dict, dict_ghg_impacts, dict_limits, autonomous_balanced=True, autonomous_gas=True,
                                                            ghg_init = ghg_init,
                                                            cost_init = cost_init
                                                            )
        
        totals_auto_bal_min_gr.T.to_excel(r"opt_results\totals_auto_bal_min_gr.xlsx")
        lca_auto_bal_min_gr.to_excel(r"opt_results\lca_auto_bal_min_gr.xlsx")
    else:
        totals_auto_bal_min_gr = pd.read_excel(r"opt_results\totals_auto_bal_min_gr.xlsx", index_col=[0]).T
        lca_auto_bal_min_gr = pd.read_excel(r"opt_results\lca_auto_bal_min_gr.xlsx", index_col=[0,1,2,3])
   
    # # 2.1 Cost-Min (bakery industry, unconstrained)
    ghg_init_bi = overview_totals_zero_bi.total_ghg.item()
    cost_init_bi = overview_totals_zero_bi.total_costs.item()

    if COST_OPT_NON_CONSTR_BI: 
        totals_cost_gr_non_constr_bi, lca_cost_gr_non_constr_bi, output_cost_gr_non_constr_bi = opt.optimize_mes(df_data_bi, 1, 0,
                                                            cost_dict, dict_ghg_impacts, dict_limits, road_transport=False,
                                                            ghg_init=ghg_init_bi, cost_init=cost_init_bi,
                                                            )
        
        totals_cost_gr_non_constr_bi.T.to_excel(r"opt_results\totals_cost_gr_non_constr_bi.xlsx")
        lca_cost_gr_non_constr_bi.to_excel(r"opt_results\lca_cost_gr_non_constr_bi.xlsx")
    else:
        totals_cost_gr_non_constr_bi = pd.read_excel(r"opt_results\totals_cost_gr_non_constr_bi.xlsx", index_col=[0]).T
        lca_cost_gr_non_constr_bi = pd.read_excel(r"opt_results\lca_cost_gr_non_constr_bi.xlsx", index_col=[0,1,2,3])

    # # 2.2 Cost-Min (bakery industry, constrained)
    if COST_MIN_OPT_BI:
        totals_cost_min_gr_bi, lca_cost_min_gr_bi, output_cost_min_gr_bi = opt.optimize_mes(df_data_bi, 1, 0,
                                                            cost_dict, dict_ghg_impacts, dict_limits_bi, road_transport=False,
                                                            ghg_init=ghg_init_bi, cost_init=cost_init_bi,
                                                            micro_wt = True
                                                            )
        
        totals_cost_min_gr_bi.T.to_excel(r"opt_results\totals_cost_min_gr_bi.xlsx")
        lca_cost_min_gr_bi.to_excel(r"opt_results\lca_cost_min_gr_bi.xlsx")
    else:
        totals_cost_min_gr_bi = pd.read_excel(r"opt_results\totals_cost_min_gr_bi.xlsx", index_col=[0]).T
        lca_cost_min_gr_bi = pd.read_excel(r"opt_results\lca_cost_min_gr_bi.xlsx", index_col=[0,1,2,3])

    # # 10.9 Scenario -- GHG-80
    diff = totals_cost_min_gr.total_ghg.item() - totals_ghg_min_gr.total_ghg.item()
    reduction_e = 0.8
    constraint_80 = totals_ghg_min_gr.total_ghg.item() + (1-reduction_e) * diff
    totals_ghg_min_gr.total_ghg.item(), totals_cost_min_gr.total_ghg.item(), constraint_80

    if GHG_MIN_80:
        totals_cost_min_ghg_80_gr, lca_cost_min_ghg_80_gr, output_cost_min_ghg_80_gr = opt.optimize_mes(df_data, 1, 0, 
                                                            cost_dict, dict_ghg_impacts, dict_limits, 
                                                            eps_ghg_constraint=constraint_80, time_limit = 20*3600,
                                                            ghg_init = ghg_init,
                                                            cost_init = cost_init,
                                                            heuristics=True)
        
        totals_cost_min_ghg_80_gr.T.to_excel(r"opt_results\totals_cost_min_ghg_80_gr.xlsx")
        lca_cost_min_ghg_80_gr.to_excel(r"opt_results\lca_cost_min_ghg_80_gr.xlsx")
    else:
        totals_cost_min_ghg_80_gr = pd.read_excel(r"opt_results\totals_cost_min_ghg_80_gr.xlsx", index_col=[0]).T
        lca_cost_min_ghg_80_gr = pd.read_excel(r"opt_results\lca_cost_min_ghg_80_gr.xlsx", index_col=[0,1,2,3])

    # Plot the annual and weekly operation for cost and GHG optimization
    file_plot_g = r"weekly_plots\result_{}_0_False_False_0_0_True_ghg_base.xlsx".format(ASSESSMENT_YEAR)
    file_plot_c = r"weekly_plots\result_{}_1_False_False_0_0_True_cost_base.xlsx".format(ASSESSMENT_YEAR)
    times_scen = ["", "summer", "winter"]

    if GEN_FIGS:
        # for cost optimization
        for i in times_scen:
            pf.make_stack_plot(file_plot_c, list_prods_e, list_cons_e, 2 if i == 'summer' else 3, 1, plot_season = i, ax_spec=0, storage_techs=storage_e)
            if i != 'summer':
                pf.make_stack_plot(file_plot_c, list_prods_lt, list_cons_lt, 2 if i == 'summer' else 3, 1, plot_season = i, ax_spec=1, new_fig=False)

            pf.make_stack_plot(file_plot_c, list_prods_ht, list_cons_ht, 2 if i == 'summer' else 3, 1, plot_season = i, ax_spec=1 if i == 'summer' else 2 , new_fig=False)
            plt.savefig("figs/stack_plot_cost_opt_{}.png".format(i), dpi = 300, bbox_inches='tight')
        
        #for GHG optimization
        for i in times_scen:
            pf.make_stack_plot(file_plot_g, list_prods_e, list_cons_e, 3 if i == 'summer' else 4, 1, plot_season = i, ax_spec=0, storage_techs=storage_e)
            if i != 'summer':
                pf.make_stack_plot(file_plot_g, list_prods_lt, list_cons_lt, 3 if i == 'summer' else 4, 1, plot_season = i, ax_spec=1, new_fig=False)

            pf.make_stack_plot(file_plot_g, list_prods_ht, list_cons_ht, 3 if i == 'summer' else 4, 1, plot_season = i, ax_spec=1 if i == 'summer' else 2, new_fig=False)
            pf.make_stack_plot(file_plot_g, list_prods_h2, list_cons_h2, 3 if i == 'summer' else 4, 1, plot_season = i, ax_spec=2 if i == 'summer' else 3, new_fig=False, 
                            storage_techs=storage_h2)
            plt.savefig("figs/stack_plot_GHGs_opt_{}.png".format(i), dpi = 300, bbox_inches='tight')
        
    # # Pareto front
    nr_eps_c = 1 # 0, 0.5, 1
    steps = nr_eps_c + 2

    """Get difference between both optimizations regarding GHGs"""
    # Get total difference between nadir and utopia
    diff = totals_cost_min_gr.total_ghg.item() - totals_ghg_min_gr.total_ghg.item()

    """Make equally sized steps based on this difference regarding GHGs"""
    # Divide in amount of steps over the graph and calculate average delta_d
    n_p = steps
    delta_d = diff/(n_p-1)

    # Now make a list of epsilon constraints based on previous points
    list_epsilon_constraints = []

    for i in range(3,n_p+1):
        eps_constr = totals_ghg_min_gr.total_ghg.item() + delta_d * (i-2)
        list_epsilon_constraints.append(eps_constr)  

    if CALC_PARETO:
        counter = 0
        # Now, calculate tehe results using the epsilon constraints to make it a single optimization problem
        for count, i in enumerate(list_epsilon_constraints):
            print("Status [{}/{}]".format(count+1,len(list_epsilon_constraints)))
            total_result, lca_result, output_result = opt.optimize_mes(df_data, 1, 0, cost_dict, dict_ghg_impacts, 
                                                            dict_limits, eps_ghg_constraint=i, time_limit = 18*3600,
                                                                        ghg_init = ghg_init,
                                                                            cost_init = cost_init)
            # Store
            if counter == 0:
                # add initial cost optimization results to df
                df_init = pd.concat([totals_ghg_min_gr, totals_cost_min_ghg_80_gr, total_result], axis=0)
                df_init_env = pd.concat([lca_ghg_min_gr, lca_cost_min_ghg_80_gr, lca_result], axis=1)    
            else:
                df_add = total_result
                df_init = pd.concat([df_init, df_add], axis=0)

                df_add_env = lca_result
                df_init_env = pd.concat([df_init_env, df_add_env], axis=1)

            counter += 1

        # add initial cost optimization results to df from results of nidir costs
        df_init_out_gr = pd.concat([df_init, totals_cost_min_gr], axis=0)
        df_init_out_env_gr = pd.concat([df_init_env, lca_ghg_min_gr], axis=1)
        
        df_init_out_gr.T.to_excel(r"opt_results\df_init_out_gr.xlsx")
        df_init_out_env_gr.to_excel(r"opt_results\df_init_out_env_gr.xlsx")
    else:
        df_init_out_gr = pd.read_excel(r"opt_results\df_init_out_gr.xlsx", index_col=[0]).T
        df_init_out_env_gr = pd.read_excel(r"opt_results\df_init_out_env_gr.xlsx", index_col=[0,1,2,3])

    # ## 11. Cost comparison
    df_init_cost = pd.concat([overview_totals_zero_bi, #Zero
                            totals_cost_gr_non_constr_bi, # Cost opt non-constrained MB
                            totals_cost_min_gr_bi, # Cost opt constrained MB
                            overview_totals_zero, totals_cost_min_gr, totals_cost_min_gr_constr,
                            totals_ghg_min_gr,
                            totals_auto_all_min_gr, 
                            totals_auto_bal_min_gr
                            ], axis=0).T  

    replacement_list_str = ['BAU (BI)', 'Cost-Min (BI)', 'Cost-Min\n-Constr (BI)', 'BAU', 
                                    'Cost-Min', 'Cost-Min\n-Constr', 'GHG-Min', 'Off-Grid', 
                                    "Balanced\nAutonomy"]
                                    
    df_init_cost.columns.values[:len(replacement_list_str)] = replacement_list_str

    df_select_cost = df_init_cost.loc[
                            df_init_cost.index.str.contains("an_costs")&
                            ~df_init_cost.index.str.contains("total_costs")]/1e3

    df_select_cost = df_select_cost[df_select_cost.sum(axis=1) != 0]

    df_select_cap = df_init_cost.loc[
                            df_init_cost.index.str.contains("cap|p_grid_connection|grid_elect_delta|natural_gas_delta|total_costs|total_ghg|curtailed")&
                            ~df_init_cost.index.str.contains("capex")]
    df_select_cap.to_excel("results\capacity overview.xlsx")

    df_select_env = df_init_cost.loc[df_init_cost.index.str.contains("an_ghg")]/1e3
    df_select_env = df_select_env[df_select_env.sum(axis=1) != 0]

    if CALC_ALL_LCA_IMPACTS:
        # ## 12. Environmental results
        df_init_lca = pd.concat([lca_results_zero_bi, #Zero
                                lca_cost_gr_non_constr_bi, # Cost opt non-constrained MB
                                lca_cost_min_gr_bi, # Cost opt constrained MB
                                lca_results_zero, lca_cost_min_gr, lca_cost_min_gr_constr,
                                lca_ghg_min_gr, 
                                lca_auto_all_min_gr, 
                                lca_auto_bal_min_gr
                                ], axis=1)  

        df_init_lca.columns.values[:len(replacement_list_str)] = replacement_list_str

        # Plot figures for all environmental impact categories
        for i, cat in enumerate(list(my_methods)):
            pf.plot_lca_impact_category(cat, df_init_lca)

        # Plot spider graph based on LCA results, for the entire MES under study
        pf.generate_radar_plot_lca(df_init_lca.iloc[:, 3:], df_init_cost.iloc[:, 3:], nr_rows=3)

        # Results per config
        for col in list(df_init_lca.columns):
            pf.generate_ind_radar_plot_lca(col, df_init_lca)

    # Plot full figure with costs and GHG emissions
    pf.generate_cost_ghg_figure(df_select_cost, df_select_env, df_init_cost)

    # Plot Pareto front for entire MES results.
    pf.plot_pareto_front(df_init_out_gr, [totals_auto_all_min_gr, totals_auto_bal_min_gr], ["Off-Grid", ""] )

    # # 14.1 Grid network cost increase and plot
    df_sens_op_grid_capex_gr, df_sens_op_grid_capex_env_gr = opt.parameter_increase_sensitivity(CALC_SENS_OP_GRID_CAPEX, dict_ghg_impacts, param_change='grid_capex', \
        max_nr=30, step_inc=400, export_results=True, mip_gap=0.0001)

    pf.plot_sensitivity_results_param_inc(df_sens_op_grid_capex_gr, 'grid_capex', ['cap_wind_on', 'cap_pv', 'p_grid_connection', 'cap_bat_en'])

    # # 14. Sensitivity analysis, influence parameters on OFF-GRID systems.
    """
    THIS IS FOR OFF-GRID SYSTEMS
    """
    search_mapping_for_sensitivity = {
                    "dr": "Discount rate",
                    "diesel_price": "Diesel price",
                    "capex": "Capex",
                    "lt": "Lifetime",
                    "om": "O&M costs",
                    "ghi": "Irradiance",
                    "temp_air": "Temperature",
                    "wind_speed": "Windspeed",
                    "km_per_day": "Demand - BEV",
                    "electricity_demand": "Demand - electricity",
                    "bld_heat_demand": "Demand - residential heat",
                    "ind_heat_demand": "Demand - industrial heat",
                    }

    # (1) Local sensitivity analysis, +10%
    if CALC_SENS:
        for counter, sens_parameter in enumerate(list(search_mapping_for_sensitivity.keys())):
            total_result, lca_result = opt.local_sens_analysis(dict_ghg_impacts, \
                sens_parameter, \
                    +0.1, totals_auto_all_min_gr.total_costs.item(), totals_auto_all_min_gr.total_ghg.item())

            total_result.index = [search_mapping_for_sensitivity[sens_parameter]]

            # Store results
            if counter == 0:
                df_init_sens = total_result.copy()
            else:
                df_init_sens = pd.concat([df_init_sens, total_result], axis=0)#, ignore_index=True)

        # Store initial cost optimization results
        df_sens_min_gr = df_init_sens
        df_sens_min_gr.T.to_excel(r"opt_results\df_sens_min_gr.xlsx")
    else:
        df_sens_min_gr = pd.read_excel(r"opt_results\df_sens_min_gr.xlsx", index_col=[0]).T

    # (2) Local sensitivity analysis, -10%
    if CALC_SENS:
        for counter, sens_parameter in enumerate(list(search_mapping_for_sensitivity.keys())):
            total_result, lca_result = opt.local_sens_analysis(dict_ghg_impacts, \
                sens_parameter, \
                    -0.1, totals_auto_all_min_gr.total_costs.item(), totals_auto_all_min_gr.total_ghg.item())

            total_result.index = [search_mapping_for_sensitivity[sens_parameter]]

            # Store results
            if counter == 0:
                df_init_sens = total_result.copy()
            else:
                df_init_sens = pd.concat([df_init_sens, total_result], axis=0)#, ignore_index=True)

        # Store results
        df_sens_out_gr = df_init_sens
        df_sens_out_gr.T.to_excel(r"opt_results\df_sens_out_gr.xlsx")
    else:
        df_sens_out_gr = pd.read_excel(r"opt_results\df_sens_out_gr.xlsx", index_col=[0]).T

    # Plot the results of the local sensitivity analysis, here +10% and -10%
    pf.plot_sensitivity_results(df_sens_min_gr, df_sens_out_gr)

    # # 14. Sensitivity analysis, influence parameters on off-grid systems.
    search_mapping_for_sensitivity_ext = {
                    "bat_en_capex": {"Battery energy capex":-0.5},
                    "chp_mix_capex": {"Advanced CHP capex":-0.5},
                    "electr_capex": {"Electrolyzer capex":-0.5},
                    "pv_capex": {"Solar PV capex":-0.5},
                    "wind_on_capex": {"Onshore wind capex":-0.5},
                    "wind_off_capex": {"Offshore wind capex":-0.5}
                    }

    # Local sensitivity analysis with more extreme value focusing on system design, +-50%
    if CALC_SENS_EXTREMES:
        for counter, (sens_parameter, ch) in enumerate(search_mapping_for_sensitivity_ext.items()):
            total_result, lca_result = opt.local_sens_analysis(dict_ghg_impacts, \
                sens_parameter, \
                    list(ch.values())[0], totals_auto_all_min_gr.total_costs.item(), totals_auto_all_min_gr.total_ghg.item())

            total_result.index = [list(ch.keys())[0]]

            # Store results
            if counter == 0:
                df_init_sens = total_result.copy()
            else:
                df_init_sens = pd.concat([df_init_sens, total_result], axis=0)

        # Store results
        df_sens_out_ext = df_init_sens
        df_sens_out_ext.T.to_excel(r"opt_results\df_sens_out_ext.xlsx")
    else:
        df_sens_out_ext = pd.read_excel(r"opt_results\df_sens_out_ext.xlsx", index_col=[0]).T

    # Include reference off-grid solution for comparison
    df_sens_out_ext_init = pd.concat([totals_auto_all_min_gr, df_sens_out_ext], axis=0)
    
    # Export capacity of system design to Excel
    df_filtered = df_sens_out_ext_init.filter(regex='(?:' + '|'.join(["cap", "p_grid_connection", "grid_elect_delta", "curtailed", \
        "natural_gas_delta", "total_costs"]) + ')').filter(regex='^(?!.*(?:capex))').T
    df_filtered.to_excel("results\capacity_overview_sens.xlsx")

    # Plot heat map for data visualization
    pf.create_heatmap(df_filtered)

    # Export dictionary with capacity limits to include it in the SI of the paper
    energy_tech_dict_full_names = {full_names_mapping_dict[key]: value for key, value in dict_limits.items()}
    pd.DataFrame(list(energy_tech_dict_full_names.items())\
        , columns=['Technology', 'Value [MW(h)]']).to_excel('results\energy_technologies_limit_dict.xlsx')
import pandas as pd
import numpy as np
from private_keys import API_KEY_ENTSOE, KEY_PREMISE

"""Define vars"""
DAYS = 365
MJ_KG_H2 = 120

# KEYS
API_KEY_ENTSOE = API_KEY_ENTSOE # You need a key from ENTSOE to retrieve day-ahead electricity prices
KEY_PREMISE = KEY_PREMISE # You need a key for premise to generate prospective LCA databases
"""YOU WILL ALSO NEED TO HAVE A LOCAL LICENSE KEY FOR GUROBI, see: https://www.gurobi.com/solutions/licensing/ & \
    https://support.gurobi.com/hc/en-us/articles/12872879801105-How-do-I-retrieve-and-set-up-a-Gurobi-license-
"""

# Get cost data.
FILE_NAME_COSTS = r"new_data\technology_costs.xlsx"
COST_DATA = pd.read_excel(FILE_NAME_COSTS, sheet_name = 'costs', index_col = "Parameter", usecols = [0,1,2], nrows=145)
FILE_NAME_ENERGY_IND = r"new_data\CONFIDENTIAL_energy_data_local_industry.xlsx" # one has to specify its own file with industrial power and heat demand data.
FILE_ELECT_DEMAND_RES = r"new_data\elect_demand_res.csv"

# Year of assesment
ASSESSMENT_YEAR = 2022
ASSESS_COUNTRY = "GR" # Greece
LOC_ELECT = ASSESS_COUNTRY  # Example value
LATITUDE = 35.3064
LONGITUDE = 23.5404
FIT = 0.06 # feed-in-tariff, in euro/kWh
SHARE_E_RES = 0.2 # Share of residential heating with electric heating \ 
            #rest here is assumed with oil boilers since they have the biggest share of heating in Greece.
HP_SHARE = 0 # Share of residential heating with air-source heat pumps
RESIDENTS = 1233 # persons in the entire MES

# Decide on maximum capacity of technologies
MAX_CAP_TECHS = 10 #MW
MAX_CAP_TECHS_RES = 5 #MW, for residential purposes
BIOMASS_ENERGY = (8041.61*1.5/2) * 4.5 #in MWh/a, 8,041.61 hectares tonnes olive trees (1.5 tonnes/ha) every 2 years, energy density of 4.5 MWh/t.
HEAT_STORAGE_RES = 0 #MWh, Maximum MWh possible heat storage for residents, we do not consider heat storage for Crete
MAX_H2_EXPORT = 365 #t/a, tonnes, 10 tonnes per day is assumed in other MES paper, however, the case study area is much smaller here.
SHARE_BEV_INIT = 0 # share of initial BEV adoption is zero in Crete.
MAX_GRID_CAP = 5 #MW, Max grid connection capacity
MIN_ELECTROLYZER_CAP = 0.05 #MW, if installed, the minimum capacity installed of the electrolyzer.
MIN_CAP_IND = 0.05  #MW, if installed, the minimum capacity installed of \
                    #high-temperature heat conversion units, needed to ensure \
                    #that not more than two units are installed
AV_COP_ASHP = 3 # It is better to take the power capacity of heat pumps in optimization, \ 
                # not heat capacity. However, LCA activity is in heat capacity... \
                # here we assume for simplicity a COP of 3 for ASHP to convert LCA heat cap to power cap.

# If using a brightway project when calculating LCA impacts yourself
# This can be enabled in case:
# 1. brightway2 and premise are installed 
# 2. One has access to ecoinvent.
# Alternatively, one can set to False, uses exported life cycle GHG emission data (without other impacts)
CALC_ALL_LCA_IMPACTS = True
GENERATE_NEW_LCA_DB = True
PROJECT_NAME = "mes_optimization" #"design_islands_391_GR"
DB_NAME = 'db_energy_system' # Temp DB name for LCA activity
EI_VERSION = "3.9.1"
NAME_REF_DB = "ecoinvent_{}_reference".format(EI_VERSION).replace(".","")
DB_NAME_INIT = 'ecoinvent {} cutoff'.format(EI_VERSION)
LOCATION_DB = r"C:\Users\tterlouw\Ecoinvent\ecoinvent_3_9_1_cutoff\datasets"
DB_NAME_CONS = 'ecoinvent {} consequential'.format(EI_VERSION) 
LOCATION_DB_CONS = r"C:\Users\tterlouw\Ecoinvent\ecoinvent_3_9_1_cutoff_cons\datasets"
CC_METHOD = ('IPCC 2021', 'climate change', 'global warming potential (GWP100)') # standard method to calculate CC impacts.
# Get total number of households in the MES:
TOTAL_HH = np.ceil(RESIDENTS / COST_DATA[NAME_REF_DB].res_home)
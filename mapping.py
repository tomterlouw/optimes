from config import CC_METHOD

my_methods = [
    CC_METHOD,
    ("Own method", "Land transformation", "Land transformation"),
    ('EF v3.1', 'acidification', 'accumulated exceedance (AE)'),
    ('EF v3.1', 'ecotoxicity: freshwater', 'comparative toxic unit for ecosystems (CTUe)'),
    ('EF v3.1', 'ecotoxicity: freshwater, inorganics', 'comparative toxic unit for ecosystems (CTUe)'),
    ('EF v3.1', 'ecotoxicity: freshwater, organics', 'comparative toxic unit for ecosystems (CTUe)'),
    ('EF v3.1', 'energy resources: non-renewable', 'abiotic depletion potential (ADP): fossil fuels'),
    ('EF v3.1', 'eutrophication: freshwater', 'fraction of nutrients reaching freshwater end compartment (P)'),
    ('EF v3.1', 'eutrophication: marine', 'fraction of nutrients reaching marine end compartment (N)'),
    ('EF v3.1', 'eutrophication: terrestrial', 'accumulated exceedance (AE)'),
    ('EF v3.1', 'human toxicity: carcinogenic', 'comparative toxic unit for human (CTUh)'),
    ('EF v3.1', 'human toxicity: carcinogenic, inorganics', 'comparative toxic unit for human (CTUh)'),
    ('EF v3.1', 'human toxicity: carcinogenic, organics', 'comparative toxic unit for human (CTUh)'),
    ('EF v3.1', 'human toxicity: non-carcinogenic', 'comparative toxic unit for human (CTUh)'),
    ('EF v3.1', 'human toxicity: non-carcinogenic, inorganics', 'comparative toxic unit for human (CTUh)'),
    ('EF v3.1', 'human toxicity: non-carcinogenic, organics', 'comparative toxic unit for human (CTUh)'),
    ('EF v3.1', 'ionising radiation: human health', 'human exposure efficiency relative to u235'),
    ('EF v3.1', 'material resources: metals/minerals', 'abiotic depletion potential (ADP): elements (ultimate reserves)'),
    ('EF v3.1', 'ozone depletion', 'ozone depletion potential (ODP)'),
    ('EF v3.1', 'particulate matter formation', 'impact on human health'),
    ('EF v3.1', 'photochemical oxidant formation: human health', 'tropospheric ozone concentration increase'),
    ('EF v3.1', 'water use', 'user deprivation potential (deprivation-weighted water consumption)'),
    ]

name_dict_fig = {
    "an_costs_op_ng": "Operation – natural gas",
    "an_costs_op_ng_res": "Operation – natural gas (residents)",
    "an_costs_capex_boiler": "Investment – gas boiler",
    "an_costs_capex_gb_res": "Investment – res. gas boiler",
    "an_costs_op_grid_abs": "Operation – grid power absorption",
    "an_costs_op_oil_res": "Operation – oil",
    "an_costs_om": "Operation & Maintenance",
    "an_costs_rep": "Replacement costs",
    "an_costs_op_co2": "CO$_{2}$ tax",
    "an_costs_op_grid_inj": "Operation – grid power injection",
    "an_costs_capex_hp": "Investment – heat pump",
    "an_costs_capex_pv": "Investment – PV system",
    "an_costs_capex_bat_en": "Investment – battery energy system",
    "an_costs_capex_bat_p": "Investment – battery power system",
    "an_costs_capex_heat_storage": "Investment – heat storage",
    "an_costs_capex_grid": "Investment – electricity grid",
    "an_costs_capex_p_solar_heat": "Investment – solar thermal",
    "an_costs_op_export_h2": "Export – hydrogen",
    "an_costs_capex_chp_mixer": "Investment – CHP",
    "an_costs_capex_gasif": "Investment – gasifier",
    "an_costs_capex_wind_on": "Investment – onshore wind",
    "an_costs_capex_wind_off": "Investment – offshore wind",
    "an_costs_capex_electrolyzer": "Investment – electrolyzer",
    "an_costs_capex_fc": "Investment – fuel cell",
    "an_costs_capex_h2_ves": "Investment – H$_{2}$ storage",
    "an_costs_op_gasoline_cars": "Operation – gasoline vehicles",
    "Investment – gas boiler_res": "Investment – gas boiler (residents)",
    "Investment – gas boiler_oil_res": "Investment – oil boiler (residents)",
    "an_costs_capex_trans_gas": "Investment – gasoline vehicles",
    "an_costs_op_d_boiler": "Operation – diesel boiler",
    "an_costs_capex_trans_bev": "Investment – BEVs",
    "an_costs_capex_e_boiler": "Investment – electric boiler (HT)",
    "an_costs_op_bev_cars": "Operation – BEVs",
    "an_costs_op_biogas": "Operation – biogas",
    "Operation – natural gas_res": "Operation – natural gas (res.)",
    "an_costs_op_wood": "Operation – wood gasification",
    "an_costs_capex_e_heating": "Investment – electric heating",

    'an_ghg_op_wood': "Operation – wood gasification", 
    'an_ghg_op_d_boiler': "Operation – diesel boiler", 
    'an_ghg_op_ng_res': "Operation – natural gas (residents)",
    'an_ghg_op_oil_res': "Operation – oil",
    'an_ghg_op_grid_abs': "Operation – grid power absorption", 
    'an_ghg_op_grid_inj': "Operation – grid power injection",
    'an_ghg_op_h2': "Export – hydrogen", 
    'an_ghg_op_gasoline_cars': "Operation – gasoline vehicles", 
    'an_ghg_op_chp': "Operation – advanced CHP unit",
    'an_ghg_boiler': "Construction – diesel boiler",
    'an_ghg_hp': "Construction – heat pump",
    'an_ghg_pv': "Construction – PV system",
    'an_ghg_bat_en': "Construction – battery energy system",
    'an_ghg_bat_p': "Construction – battery power system", 
    'an_ghg_heat_storage': "Construction – heat storage",
    'an_ghg_chp_mixer': "Construction – CHP",
    'an_ghg_gasif': "Construction – gasifier", 
    'an_ghg_grid_ins': "Construction – electricity grid",
    'an_ghg_p_solar_heat': "Construction – solar thermal",
    'an_ghg_wind_on': "Construction – onshore wind",
    'an_ghg_wind_off': "Construction – offshore wind",
    'an_ghg_e_heating': "Construction – electric heating",
    'an_ghg_electrolyzer': "Construction – electrolyzer",
    'an_ghg_fc': "Construction – fuel cell", 
    'an_ghg_h2_ves': "Construction – H$_{2}$ storage",
    'an_ghg_e_boiler': "Construction – electric boiler (HT)",
    'an_ghg_gb_res': "Construction – res. gas boiler",
    'an_ghg_bev_cars': "Construction – BEVs",
}

cap_mapping = {
    'cap_wind_off': 'Offshore wind',
    'cap_wind_on': 'Onshore wind',
    'cap_pv': 'Solar PV',
    'cap_bat_en': 'Battery energy storage',
    'cap_bat_p': 'Battery power',
    'cap_h2_ves': 'Hydrogen storage',
    'cap_electrolyzer': 'Electrolyzer',
    'cap_fc': 'Fuel cell',
    'cap_hp': 'Heat pump',
    'cap_d_boiler': 'Diesel boiler',
    'cap_heat_storage': 'Heat storage',
    'cap_solar_heat': 'Solar heat',
    'cap_chp_mix': 'Advanced CHP',
    'cap_gasif': 'Gasifier',
    'cap_gb_res': 'Residential gas boiler',
    'cap_e_heating': 'Electric heating',
    'cap_e_boiler': 'Electric boiler',
    'p_grid_connection': 'Power grid connection',
    'ratio_pv_curtailed': 'Curtailment – solar PV',
    'ratio_wind_curtailed_off': 'Curtailment – offshore wind',
    'ratio_wind_curtailed_on': 'Curtailment – onshore wind',
    'total_costs': 'Annual costs [M€/a]',
}

full_names_mapping_dict = {
    'max_pv': 'max. solar PV',
    'max_wind_on': 'max. wind onshore',
    'max_wind_off': 'max. wind offshore',
    'max_bm_annu': 'max. biomass annual',
    'max_grid_cap': 'max. power grid capacity',
    'st_max': 'max. solar thermal',
    'share_bev': 'share of battery electric vehicles',
    'max_h2_export': 'max. hydrogen export',
    'biogas_potential_MWh': 'biogas potential (MWh)',
    'max_solar_heat_ha': 'max. solar heat in hectare',
    'min_solar_heat_ha': 'min. solar heat in hectare',
    'min_cap_wind_off': 'min. wind offshore capacity',
    'min_cap': 'overall min. capacity',
    'min_cap_ind': 'min. industrial technology capacity',
    'min_cap_res': 'min. residential technology capacity',
    'max_bat': 'max. battery capacity',
    'max_h2_storage': 'max. hydrogen storage capacity',
    'max_electrolyzer': 'max. electrolyzer capacity',
    'min_electrolyzer': 'min. electrolyzer capacity',
    'max_fc': 'max. fuel cell capacity',
    'max_hp': 'max. heat pump capacity',
    'max_boiler': 'max. industrial boiler capacity',
    'max_boiler_res': 'max. residential boiler capacity',
    'max_heat_storage': 'max. heat storage capacity',
    'max_chp_mix': 'max. CHP mix capacity',
    'min_chp_mix': 'min. CHP mix capacity',
    'max_gasif': 'max. gasification capacity',
    'households': 'number of households'
}

column_dict_categories = {
    'Land transformation': "LT",
    'acidification': "AC",
    'climate change': "CC",
    'ecotoxicity: freshwater': "ET$_{F}$",
    'ecotoxicity: freshwater, inorganics': "ET$_{FI}$",
    'ecotoxicity: freshwater, organics': "ET$_{FO}$",
    'energy resources: non-renewable': "ER",
    'eutrophication: freshwater': "EF$_{F}$",
    'eutrophication: marine': "EF$_{M}$",
    'eutrophication: terrestrial': "EF$_{T}$",
    'human toxicity: carcinogenic': "HT$_{C}$",
    'human toxicity: carcinogenic, inorganics': "HT$_{CI}$",
    'human toxicity: carcinogenic, organics': "HT$_{CO}$",
    'human toxicity: non-carcinogenic': "HT$_{NC}$",
    'human toxicity: non-carcinogenic, inorganics': "HT$_{NCI}$",
    'human toxicity: non-carcinogenic, organics': "HT$_{NCO}$",
    'ionising radiation: human health': "IR",
    'material resources: metals/minerals': "MM",
    'ozone depletion': "OD",
    'particulate matter formation': "PM",
    'photochemical oxidant formation: human health': "PF",
    'water use': "WU",
}

color_op_ng = 'chocolate'
color_op_grid_abs = "red"
color_op_grid_inj = "#355E3B"
color_op_chp = "green"
color_an_rep = "yellow"
color_an_om = "turquoise"
color_prod_boiler = 'brown'
color_prod_grid = "black"
color_op_gasoline = 'white'
color_prod_off_wind = 'aqua'
color_prod_on_wind = "aquamarine"
color_prod_pv = "orange"
color_prod_bat = "gold"
color_exp_h2 = "darkblue"
color_prod_heat_stor = 'indigo'
color_prod_chp = "tomato"
color_prod_gasif = 'deeppink'
color_prod_hp = 'red'
color_prod_solar_th = 'chartreuse'
color_prod_h2_stor = 'khaki'
color_prod_elect = 'darkcyan'
color_prod_fc = 'pink'
color_prod_bev = "plum"
color_prod_car_gas = "white"
color_syngas = 'lightgreen'
color_co2 = 'darkgrey'
color_e_heating = 'magenta'
color_biogas = '#7CFC00'
color_oil = 'saddlebrown'
color_diesel = 'brown'

color_dict_cost = {'an_costs_op_ng':color_op_ng,
                   "an_costs_op_ng_res":color_op_ng,
                 'an_costs_op_grid_abs': color_op_grid_abs,
                 "an_costs_op_wood": color_syngas,
                 "an_costs_op_biogas": color_biogas,
                 'an_costs_capex_boiler':color_prod_boiler,
                 'an_costs_capex_gb_res':color_prod_boiler,
                 'an_costs_capex_grid_ins':color_prod_grid,
                 'an_costs_rep':color_an_rep,
                 'an_costs_om': color_an_om,
                 'an_costs_op_grid_inj': color_op_grid_inj,
                 'an_costs_op_export_h2': color_exp_h2,
                 'an_costs_capex_hp': color_prod_hp,
                 'an_costs_capex_pv': color_prod_pv,
                 'an_costs_capex_bat_en': color_prod_bat,
                 'an_costs_capex_bat_p': color_prod_bat,
                 'an_costs_capex_heat_storage': color_prod_heat_stor,
                 'an_costs_capex_chp_mixer': color_prod_chp,
                 'an_costs_capex_gasif': color_prod_gasif,
                 'an_costs_capex_grid':color_prod_grid,
                 'an_costs_capex_p_solar_heat': color_prod_solar_th,
                 'an_costs_capex_wind_on':color_prod_on_wind,
                 'an_costs_capex_wind_off': color_prod_off_wind,
                 'an_costs_capex_electrolyzer': color_prod_elect,
                 'an_costs_capex_fc': color_prod_fc,
                 'an_costs_capex_h2_ves': color_prod_h2_stor,
                   'an_costs_capex_boiler_oil_res':color_oil,
                   'an_costs_capex_e_boiler':color_e_heating,
                 'an_costs_op_co2': color_co2,
                 "an_costs_op_gasoline_cars": color_op_gasoline,
                   "an_costs_cars_gasoline": color_op_gasoline,
                   "an_costs_op_bev_cars": color_prod_bev,
                   "an_costs_capex_trans_bev": color_prod_bev,
                   "an_costs_capex_trans_gas": color_prod_car_gas,
                   "an_costs_capex_e_heating": color_e_heating,
                   "an_costs_op_oil_res":color_oil,
                   "an_costs_op_d_boiler":color_diesel,
                  }

color_dict_env = {
                 'Construction, Natural gas boiler':color_prod_boiler,
                 'Operation, Electricity': color_op_grid_abs,
                 'Operation, heat from natural gas':color_op_ng,
                 "Operation, gas turbine":color_op_chp,
                 'Construction, grid electricity network':color_prod_grid,
                 'Transportation residents, petrol': color_op_gasoline,
                 'Construction, wind offshore': color_prod_off_wind,
                 'Construction, Wind onshore':color_prod_on_wind,
                 'Construction, PV': color_prod_pv,
                 'Construction, battery': color_prod_bat,
                 'Construction, hydrogen storage': color_prod_h2_stor,
                 'Construction, electrolyzer': color_prod_elect,
                 'Construction, fuel cell': color_prod_fc,
                 'Construction, heat pump': color_prod_hp,
                 'Construction, heat storage': color_prod_heat_stor,
                 'Construction, solar heat': color_prod_solar_th,
                 'Construction, gas turbine': color_prod_chp,
                 'Construction, wood gasification': color_prod_gasif,
                 'Operation, syngas production': color_syngas,
                 'Construction, wind onshore':color_prod_on_wind,
                 'Construction, heat pump': color_prod_hp,
                 'Construction, natural gas boiler':color_prod_boiler,
                 'Operation, electricity': color_op_grid_abs,
                 'Credit, hydrogen': color_exp_h2,
                 "Credit, electricity injection": color_op_grid_inj,
                  "Construction, battery electric vehicles": color_prod_bev,
                  "Construction, electric heating": color_e_heating,
                "Operation - heat from oil boiler": color_oil,
                "Operation - heat from diesel boiler":color_diesel,
                'an_ghg_op_wood': color_syngas,
                'an_ghg_op_d_boiler':color_diesel,
                'an_ghg_op_ng_res':color_op_ng,
                'an_ghg_op_oil_res': color_oil,
                'an_ghg_op_grid_abs': color_op_grid_abs,
                'an_ghg_op_grid_inj': color_op_grid_inj,
                'an_ghg_op_h2': color_exp_h2,
                'an_ghg_op_gasoline_cars': color_op_gasoline,
                'an_ghg_op_chp':color_op_chp,
                'an_ghg_boiler':color_prod_boiler,
                'an_ghg_hp': color_prod_hp,
                'an_ghg_pv': color_prod_pv,
                'an_ghg_bat_en': color_prod_bat,
                'an_ghg_bat_p': color_prod_bat,
                'an_ghg_heat_storage': color_prod_heat_stor,
                'an_ghg_chp_mixer': color_prod_chp,
                'an_ghg_gasif': color_prod_gasif,
                'an_ghg_grid_ins':color_prod_grid,
                'an_ghg_p_solar_heat': color_prod_solar_th,
                'an_ghg_wind_on':color_prod_on_wind,
                'an_ghg_wind_off': color_prod_off_wind,
                'an_ghg_e_heating': color_e_heating,
                'an_ghg_electrolyzer': color_prod_elect,
                'an_ghg_fc': color_prod_fc,
                'an_ghg_h2_ves': color_prod_h2_stor,
                'an_ghg_e_boiler':color_e_heating,
                'an_ghg_gb_res':color_prod_boiler,
                'an_ghg_bev_cars': color_prod_bev,
    
             }

hex_op_ng = ''
hex_op_grid_abs = ""
hex_op_grid_inj = ""
hex_op_chp = ""
hex_an_rep = "" 
hex_an_om = ""
hex_prod_boiler = "/"
hex_prod_grid = "/"
hex_op_gasoline =  ""
hex_prod_off_wind = "/"
hex_prod_on_wind = "/"
hex_prod_pv = "/"
hex_prod_bat = "/"
hex_exp_h2 = ""
hex_prod_heat_stor = "/"
hex_prod_chp = "/"
hex_prod_gasif = "/"
hex_prod_hp = "/"
hex_prod_solar_th = "/"
hex_prod_h2_stor = "/"
hex_prod_elect = "/"
hex_prod_fc = "/"
hex_prod_bev = "/"
hex_prod_gas = "/"
hex_prod_e_heating = "/"
hex_syngas = ''
hex_biogas = ''
hex_co2 = ''
hex_oil = ""
hex_prod_oil = "/"
hex_diesel=""

hex_dict_cost = {'an_costs_op_ng':hex_op_ng,
                "an_costs_op_ng_res":hex_op_ng,
                'an_costs_op_grid_abs': hex_op_grid_abs,
                "an_costs_op_wood": hex_syngas,
                "an_costs_op_biogas": hex_biogas,
                'an_costs_capex_boiler':hex_prod_boiler,
                'an_costs_capex_gb_res':hex_prod_boiler,
                'an_costs_capex_grid_ins':hex_prod_grid,
                'an_costs_rep':hex_an_rep,
                'an_costs_om': hex_an_om,
                'an_costs_op_grid_inj': hex_op_grid_inj,
                'an_costs_op_export_h2': hex_exp_h2,
                'an_costs_capex_hp': hex_prod_hp,
                'an_costs_capex_pv': hex_prod_pv,
                'an_costs_capex_bat_en': hex_prod_bat,
                'an_costs_capex_bat_p': hex_prod_bat,
                'an_costs_capex_heat_storage': hex_prod_heat_stor,
                'an_costs_capex_chp_mixer': hex_prod_chp,
                'an_costs_capex_gasif': hex_prod_gasif,
                'an_costs_capex_grid':hex_prod_grid,
                'an_costs_capex_p_solar_heat': hex_prod_solar_th,
                'an_costs_capex_wind_on':hex_prod_on_wind,
                'an_costs_capex_wind_off': hex_prod_off_wind,
                'an_costs_capex_electrolyzer': hex_prod_elect,
                'an_costs_capex_fc': hex_prod_fc,
                'an_costs_capex_h2_ves': hex_prod_h2_stor,
                'an_costs_op_co2': hex_co2,
                "an_costs_op_gasoline_cars": hex_op_gasoline,
                "an_costs_cars_gasoline": hex_op_gasoline,
                "an_costs_capex_trans_gas": hex_prod_gas,
                "an_costs_op_bev_cars": "",
                "an_costs_capex_trans_bev": hex_prod_bev,
                "an_costs_capex_e_heating": hex_prod_e_heating,
                "an_costs_op_oil_res":hex_oil,
                'an_costs_capex_boiler_oil_res':hex_prod_oil,
                'an_costs_capex_e_boiler':hex_prod_e_heating,
                "an_costs_op_d_boiler":hex_diesel,
                  }

hex_dict_env = {
                'Construction, Natural gas boiler':hex_prod_boiler,
                'Operation, Electricity': hex_op_grid_abs,
                'Operation, heat from natural gas':hex_op_ng,
                "Operation, gas turbine":hex_op_chp,
                'Construction, grid electricity network':hex_prod_grid,
                'Transportation residents, petrol': hex_op_gasoline,
                'Construction, wind offshore': hex_prod_off_wind,
                'Construction, Wind onshore':hex_prod_on_wind,
                'Construction, PV': hex_prod_pv,
                'Construction, battery': hex_prod_bat,
                'Construction, hydrogen storage': hex_prod_h2_stor,
                'Construction, electrolyzer': hex_prod_elect,
                'Construction, fuel cell': hex_prod_fc,
                'Construction, heat pump': hex_prod_hp,
                'Construction, heat storage': hex_prod_heat_stor,
                'Construction, solar heat': hex_prod_solar_th,
                'Construction, gas turbine': hex_prod_chp,
                'Construction, wood gasification': hex_prod_gasif,
                'Operation, syngas production': hex_syngas,
                'Construction, wind onshore':hex_prod_on_wind,
                'Construction, heat pump': hex_prod_hp,
                'Construction, natural gas boiler':hex_prod_boiler,
                'Operation, electricity': hex_op_grid_abs,
                'Credit, hydrogen': hex_exp_h2,
                "Credit, electricity injection": hex_op_grid_inj,
                "Construction, battery electric vehicles": hex_prod_bev,
                "Construction, electric heating": hex_prod_e_heating,
                "Operation - heat from oil boiler": hex_oil,
                "Operation - heat from diesel boiler":hex_diesel,
                'an_ghg_op_wood': hex_syngas,
                'an_ghg_op_d_boiler':hex_diesel,
                'an_ghg_op_ng_res':hex_op_ng,
                'an_ghg_op_oil_res': hex_oil,
                'an_ghg_op_grid_abs': hex_op_grid_abs,
                'an_ghg_op_grid_inj': hex_op_grid_inj,
                'an_ghg_op_h2': hex_exp_h2,
                'an_ghg_op_gasoline_cars': hex_op_gasoline,
                'an_ghg_op_chp':hex_op_chp,
                'an_ghg_boiler':hex_prod_boiler,
                'an_ghg_hp': hex_prod_hp,
                'an_ghg_pv': hex_prod_pv,
                'an_ghg_bat_en': hex_prod_bat,
                'an_ghg_bat_p': hex_prod_bat,
                'an_ghg_heat_storage': hex_prod_heat_stor,
                'an_ghg_chp_mixer': hex_prod_chp,
                'an_ghg_gasif': hex_prod_gasif,
                'an_ghg_grid_ins':hex_prod_grid,
                'an_ghg_p_solar_heat': hex_prod_solar_th,
                'an_ghg_wind_on':hex_prod_on_wind,
                'an_ghg_wind_off': hex_prod_off_wind,
                'an_ghg_e_heating': hex_prod_e_heating,
                'an_ghg_electrolyzer': hex_prod_elect,
                'an_ghg_fc': hex_prod_fc,
                'an_ghg_h2_ves': hex_prod_h2_stor,
                'an_ghg_e_boiler':hex_prod_e_heating,
                'an_ghg_gb_res':hex_prod_boiler,
                'an_ghg_bev_cars': hex_prod_bev,
             }

contribution_mapping_system = {
    "photovoltaic open ground installation, 570 kWp, multi-Si, on open ground":"Construction, PV",
    "li-ion (NMC)":"Construction, battery",    
    "battery management system, kWh":"Construction, battery",   
    "energy management system, kWh":"Construction, battery",        
    "power conditioning system, container system":"Construction, battery",       
    'market for wind power plant, 2MW, offshore, fixed parts':"Construction, wind offshore",     
    'market for wind power plant, 2MW, offshore, moving parts':"Construction, wind offshore",
    "market for wind turbine, 2MW, onshore":"Construction, wind onshore",
    "market for wind turbine network connection, 2MW, onshore":"Construction, wind onshore",
    "wind turbine network connection construction, 4.5MW, onshore":"Construction, grid electricity network",
    'electrolyzer production, 1MWe, PEM, Stack':"Construction, electrolyzer", 
    'electrolyzer production, 1MWe, PEM, Balance of Plant':"Construction, electrolyzer", 
    'market for electricity, low voltage':"Operation, electricity",    
    'market group for electricity, low voltage':"Operation, electricity",  
    'market for electricity, low voltage, Crete':"Operation, electricity",    
    "market for electricity, low voltage, grid injection":"Credit, electricity injection", 
    "high pressure hydrogen storage tank": "Construction, hydrogen storage",
    "heat pump, air-water, 7kW": "Construction, heat pump",
    "micro gas turbine production, 100kW electrical": "Construction, gas turbine",
    "advanced gas turbine, 400kWe": "Construction, gas turbine",
    "mixed CHP unit, operation": "Operation, gas turbine",
    "gas boiler production": "Construction, natural gas boiler",
    "heat production, natural gas, w\o infrastructure": "Operation, heat from natural gas", 
    "synthetic gas production, from wood, at fixed bed gasifier, w\o infrastructure": "Operation, syngas production",
    "synthetic gas factory construction": "Construction, wood gasification",
    "heat storage production, 2000l": "Construction, heat storage",
    "solar collector system installation, Cu flat plate collector, one-family house, combined system, w\o storage": "Construction, solar heat",
    "fuel cell system assembly, 1 kWe, proton exchange membrane (PEM)": "Construction, fuel cell",
    "market for transport, passenger car, medium size, petrol, EURO 5": "Transportation residents, petrol",
    'transport, passenger car, electric, w\o fuel': "Construction, battery electric vehicles",
    "hydrogen production, steam reforming": "Credit, hydrogen",
    "electric storage heater, 5kW": "Construction, electric heating",
    "heat production, light fuel oil, at boiler 10kW condensing, non-modulating": "Operation - heat from oil boiler",
    "heat production, light fuel oil, at industrial furnace 1MW": "Operation - heat from diesel boiler",
    }

dict_units = {
    ('IPCC 2021', 'climate change', 'global warming potential (GWP100)'): "kg CO$_{2}$-eq.",  
    ("Own method", "Land transformation", "Land transformation"): "m$^{2}$ land", 
    ('EF v3.1', 'acidification', 'accumulated exceedance (AE)'): "mol H+ -eq.", 
    ('EF v3.1', 'ecotoxicity: freshwater', 'comparative toxic unit for ecosystems (CTUe)'): "CTUe",  
    ('EF v3.1', 'ecotoxicity: freshwater, inorganics', 'comparative toxic unit for ecosystems (CTUe)'): "CTUe",  
    ('EF v3.1', 'ecotoxicity: freshwater, organics', 'comparative toxic unit for ecosystems (CTUe)'): "CTUe",  
    ('EF v3.1', 'energy resources: non-renewable', 'abiotic depletion potential (ADP): fossil fuels'): "MJ, net calorific value",  
    ('EF v3.1', 'eutrophication: freshwater', 'fraction of nutrients reaching freshwater end compartment (P)'): "kg P-eq.",
    ('EF v3.1', 'eutrophication: marine', 'fraction of nutrients reaching marine end compartment (N)'): "kg N-eq.",  
    ('EF v3.1', 'eutrophication: terrestrial', 'accumulated exceedance (AE)'): "mol N-eq." ,
    ('EF v3.1', 'human toxicity: carcinogenic', 'comparative toxic unit for human (CTUh)'): "CTUh",  
    ('EF v3.1', 'human toxicity: carcinogenic, inorganics', 'comparative toxic unit for human (CTUh)'): "CTUh",  
    ('EF v3.1', 'human toxicity: carcinogenic, organics', 'comparative toxic unit for human (CTUh)'): "CTUh",  
    ('EF v3.1', 'human toxicity: non-carcinogenic', 'comparative toxic unit for human (CTUh)'): "CTUh",  
    ('EF v3.1', 'human toxicity: non-carcinogenic, inorganics', 'comparative toxic unit for human (CTUh)'): "CTUh",  
    ('EF v3.1', 'human toxicity: non-carcinogenic, organics', 'comparative toxic unit for human (CTUh)'): "CTUh",  
    ('EF v3.1', 'ionising radiation: human health', 'human exposure efficiency relative to u235'): "kBq U$_{235}$-eq.",  
    ('EF v3.1', 'material resources: metals/minerals', 'abiotic depletion potential (ADP): elements (ultimate reserves)'): "kg Sb-eq.",  
    ('EF v3.1', 'ozone depletion', 'ozone depletion potential (ODP)'): "kg CFC-11-eq.",  
    ('EF v3.1', 'particulate matter formation', 'impact on human health'): "disease incidence",  
    ('EF v3.1', 'photochemical oxidant formation: human health', 'tropospheric ozone concentration increase'): "kg NMVOC-eq.",  
    ('EF v3.1', 'water use', 'user deprivation potential (deprivation-weighted water consumption)'): "m$^{3}$ world eq. deprived",  
             }

# This is some mapping for the stack plots
list_prods_e = ['p_grid_abs', 'p_battdis', 'p_chp_mix','pv_total', "wind_total_off", "wind_total_on", 'p_fc']
list_cons_e = ['p_battch', 'p_grid_inj',  "f_bev_1", "f_bev_2", "f_bev_3","f_elect","f_hp","electricity_demand", "f_e_heating", "f_e_boiler"]
list_prods_ht = ['p_d_boiler', 'p_th_chp_mix', 'p_e_boiler']
list_cons_ht = ["ind_heat_demand"]
list_prods_lt = ['p_fc_heat', 'p_gb_res', 'p_solar_heat', 'p_th_hp', "heat_stor_dis", "p_e_heating"]
list_cons_lt = ['heat_stor_ch', 'bld_heat_demand']
list_prods_h2 = ['p_elect', 'p_h2_ves']
list_cons_h2 = ['f_chp_h2', 'p_h2_ves', 'f_fc', 'export_h2']
storage_e = ["E_batt"]
storage_h2 = ["E_h2_ves"]
storage_th_lt = ['E_heat_stor']

dict_color_n = {'p_grid_abs':'brown',
                 'p_battdis': "crimson",
                 'p_chp_mix': 'darkgreen',
                 'pv_total': 'yellow',
                 'wind_total_off': 'navy',
                 'wind_total_on': 'lightblue',
                 'p_fc': 'coral',
                 'p_battch': "crimson",
                 'p_grid_inj': '#355E3B',
                 'f_bev_1': 'aqua',
                 'f_bev_2': '#3af2f2',
                 'f_bev_3': '#11c1c1',
                 'f_elect': 'darkorange',
                 'f_hp': 'purple',
                 'f_e_heating': 'magenta',
                 'f_e_boiler': "orchid",
                 'electricity_demand': '#838383',
                
                 'p_d_boiler': 'brown',
                 'p_th_chp_mix': 'darkgreen',
                 'p_e_boiler': "orchid",
                 'ind_heat_demand': '#838383',  
                
                 'heat_stor_dis': "red",
                 'heat_stor_ch': "darkred",
                 'p_th_hp': 'green',
                 'p_solar_heat':'orange',
                 'p_fc_heat':'coral',
                 'p_gb_res':'darkblue',
                 'p_e_heating': 'magenta',
                 'bld_heat_demand': '#838383',  
                
                 'p_elect': 'blue',
                 'f_chp_h2':'purple',
                 'p_h2_ves_dis':'red',
                 'p_h2_ves':'darkred',
                 'f_fc':'yellow',
                 'export_h2': 'darkblue',  
               }

# Define a replacement dictionary
replacement_dict = {
    "p_grid_abs": "Grid - import",
    "p_grid_inj": "Grid - export",
    "p_battdis": "Battery - discharge",
    "p_battch": "Battery - charge",
    "pv_total": "Generation - solar PV",
    "wind_total_off": "Generation - Off. wind",
    "wind_total_on": "Generation - On. wind",
    "p_chp_mix": "Conversion - adv. CHP",
    "f_bev_1": "Consumption - BEV (1)",
    "f_bev_2": "Consumption - BEV (2)",
    "f_bev_3": "Consumption - BEV (3)",
    "f_elect": "Consumption - electrolyzer",
    "f_hp": "Consumption - HP",
    "electricity_demand": "Power demand",
    "p_fc": "Conversion - fuel cell",
    "p_th_chp_mix": "Conversion - adv. CHP",
    "p_d_boiler": "Conversion - diesel boiler",
    "ind_heat_demand": "Industrial heat demand",
    "bld_heat_demand": "Residential heat demand",
    "heat_stor_ch": "Heat stor. - charge",
    "heat_stor_dis": "Heat stor. - discharge",
    "p_th_hp": "Generation - HP",
    "p_fc_heat": "Conversion - fuel cell",
    "p_gb_res": "Conversion - gas boiler",
    "p_solar_heat": "Generation - solar heat",
    "Conversion - fuel cell_heat": "Conversion - fuel cell",
    "Conversion - gas boiler_res": "Conversion - gas boiler",
    "p_elect": "Conversion - electrolyzer",
    "f_fc": "Conversion - fuel cell",
    "f_chp_h2": "H$_2$ to advanced CHP",
    "export_h2": "H$_2$ export",
    "p_h2_ves_dis": "H$_2$ storage - disch.",
    "p_h2_ves": "H$_2$ storage - ch.",
    "f_e_boiler": "Consumption - electric boiler",
    "p_e_boiler": "Generation - electric boiler",
    "f_e_heating": "Consumption - res. electrical heat",
    "p_e_heating": "Generation - res. electrical heat"}
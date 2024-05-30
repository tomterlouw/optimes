import numpy as np
import pandas as pd
import pvlib
from pvlib import pvsystem
import datetime as dt
import windpowerlib as wp
from windpowerlib import ModelChain, WindTurbine

from config import (
    API_KEY_ENTSOE, NAME_REF_DB, ASSESSMENT_YEAR, MAX_CAP_TECHS_RES, MAX_H2_EXPORT,
    MAX_CAP_TECHS, BIOMASS_ENERGY, HEAT_STORAGE_RES, MAX_H2_EXPORT, SHARE_BEV_INIT,
    MAX_GRID_CAP, MIN_ELECTROLYZER_CAP, MIN_CAP_IND, RESIDENTS, COST_DATA, TOTAL_HH,
    AV_COP_ASHP, FILE_NAME_ENERGY_IND, FILE_ELECT_DEMAND_RES, FIT
    )

# API ENTSO-E
from entsoe import EntsoePandasClient
CLIENT = EntsoePandasClient(api_key=API_KEY_ENTSOE)

class EnergyDataProcessor:
    """
    Class to handle data processing for energy optimization.

    Attributes:
        assess_country (str): The country being assessed.
        get_year (int): The year of assessment.
        lat (float): self.lat of the location.
        lon (float): self.lon of the location.
        fit (float): Feed-in-tariff.
        constrained (boolean): if True, the MES is constrained by limitations in terms of capacities that can be installed.
        residential sector (boolean): if True, the residential sector is included.
    """
    
    def __init__(self, assess_country, get_year, lat, lon, fit, constrained, residential_sector):
        """
        Initializes the EnergyDataProcessor class.

        Parameters:
            assess_country (str): The country being assessed.
            get_year (int): The year of assessment.
            lat (float): self.lat of the location.
            lon (float): self.lon of the location.
            fit (float): Feed-in-tariff.
            constrained (boolean): if True, the MES is constrained by limitations in terms of capacities that can be installed.
            residential sector (boolean): if True, the residential sector is included.
        """
        self.assess_country = assess_country
        self.get_year = get_year
        self.lat = lat
        self.lon = lon
        self.fit = fit
        self.constrained = constrained
        self.residential_sector = residential_sector

    def get_weather_data(self, start_year=2007, end_year=2016):
        """
        Gets location-specific weather data from PVGIS.

        Args:
            self.lat (float): self.lat [degrees].
            self.lon (float): self.lon [degrees].
            start_year (int): start of typical meteorogical year (TMY) data, std is 2007 [int].
            end_year (int): start of typical meteorogical year (TMY) data, std is 2016 [int].
        Returns:
            data: dataframe with annual data based on typical meteorogical year data of PVGIS, with:
                    Date & time (UTC for normal CSV, local timezone time for the EPW format)
                    temp_air - Dry bulb (air) temperature.
                    relative_humidity [%] - Relative Humidity.
                    ghi [W/m2] - Global horizontal irradiance.
                    dni [W/m2] - Direct (beam) irradiance.
                    dhi [W/m2] - Diffuse horizontal irradiance.
                    IR(h) [W/m2] - Infrared radiation downwards.
                    wind_speed [m/s] - Windspeed.
                    wind_direction [°] - Wind direction.
                    pressure [Pa] - Surface (air) pressure.
            elevation: elevation of location [meters]
        """

        all_data = pvlib.iotools.get_pvgis_tmy(float(self.lat), float(self.lon), startyear=start_year, endyear=end_year, 
                                         outputformat='csv')
        elevation = all_data[2]['elevation'] 
        data = all_data[0] 
        data = data.reset_index()
        data['time(UTC)'] = data['time(UTC)'].astype(str)

        data['time(UTC)'] = data['time(UTC)'].str[:-9]#.astype(float)
        data['time(UTC)'] = data['time(UTC)'].apply(lambda x: str(ASSESSMENT_YEAR) + x[4:])#.astype(float)
        data['time(UTC)'] = data['time(UTC)'].apply(lambda x: 
                                                    dt.datetime.strptime(str(x),'%Y-%m-%d %H:%M'))  

        data = data.set_index('time(UTC)')

        data['ghi'] = data['ghi'].fillna(0)
        data['dni'] = data['dni'].fillna(0)
        data['dhi'] = data['dhi'].fillna(0)

        data['wind_speed'][data['wind_speed']<0] = 0

        return data, elevation

    def get_renewable_profiles(self, weather_data):
        """
        Calculates renewable electricity generation profiles for solar PV, onshore wind, and offshore wind for TMY per kWp.

        Args:
            weather_data (dataframe): dataframe with weather information for one year [pd.Dataframe].
            self.lat (float): self.lat, in decimal degrees, between -90 and 90, north is positive (ISO 19115) [degrees].
            self.lon (float): self.lon, in decimal degrees, between -180 and 180, east is positive (ISO 19115) [degrees].
        Returns:
            pv_MW (np.array): the pv profile per kWp module modelled [kWp].
            wind_MW_on (np.array): the wind profile for on-shore per kWp module modelled [kWp].
            wind_MW_off (np.array): the wind profile for off-shore per kWp module modelled [kWp] [kWp].
        """

        #PV
        pv_MW_start = pv_profile_generator_tmy(weather_data, self.lat, self.lon)
        pv_MW = np.array(pv_MW_start)

        # WIND
        wind_MW_start_on = wind_profile_generator_tmy(weather_data, self.lat, self.lon, turbine_spec="V90/2000", hub_height = 80)
        wind_MW_on = np.array(wind_MW_start_on)

        wind_MW_start_off = wind_profile_generator_tmy(weather_data, self.lat, self.lon, offshore=True)
        wind_MW_off = np.array(wind_MW_start_off)

        print("pv cf '{}', wind cf off '{}', wind cf on '{}'".format(round(pv_MW.mean(),3), 
                                                                                                          round(wind_MW_off.mean(),3),
                                                                                                          round(wind_MW_on.mean(),3)))
        return pv_MW, wind_MW_off, wind_MW_on

    def get_elect_prices(self, loc, year, fit=False, share_inj = 0.25):
        """
        Gets the day-ahead price electricity price for a specific year from ENTSO-E.

        Args:
            loc (str): location [str].
            year (int): year of assessment [-].
            fit (boolean, or float): if a float, then this will set as the feed-in-tariff. Otherwise, x% of electricity cost will be taken, specified with 'share_inj' [euro/kW].
            share_inj (float):  x% of electricity cost will be taken as feed-in-tariff [-].
        Returns:
            da_prices (np.array): hourly electricity price, based on day-ahead price [kWp].
            fit_injection (np.array): hourly feed-in-tariff used.
        """
        ei_loc = loc.split('_')[-1]

        # We assess three countries here, for Norway, there are different electricity markets, here we choose 'NO_2' as this fits to the location.
        dict_entsoe_coupling = {"GR":"GR", "NO":"NO_2", "GB":"GB"}

        # Electricity prices timestamps for one year.
        start = pd.Timestamp('{}0101'.format(year), tz='utc')
        end = pd.Timestamp('{}0101'.format(year+1), tz='utc')

        # methods that return Pandas Series, convert to UTC
        df_prices = CLIENT.query_day_ahead_prices(dict_entsoe_coupling[ei_loc], start=start,end=end)
        df_prices.index = pd.to_datetime(df_prices.index, utc=True)

        # Give a warning with unexpected number of timeslots.
        if int(len(df_prices)) < 8760:
            print("WARNING: '{}' missing timestamps, they are forward filled".format(8760-len(df_prices)))
            df_prices = df_prices.resample('H').ffill().reset_index().set_index('index')
            df_prices = df_prices.iloc[:,0]

        # In great britian, GBP is more worth than euro, account fir this.
        if loc=="GB":
            # compensate for GBP currency
            da_prices = np.array((df_prices/1e3)[:8760]) * 1.15  
        else:
            da_prices = np.array((df_prices/1e3)[:8760])

        # Put into arrays
        if fit:
            fit_injection = fit      
        else:
            fit_injection = np.array(da_prices*share_inj)

        print("Average electricity price '{}' euro/kWh, export '{}' €/kWh".format(round(da_prices.mean(),4),round(np.mean(fit_injection),4)))

        return da_prices, fit_injection

    def get_energy_demands(self, weather_data):
        """
        Gets energy demands for a specific region.

        Parameters:
            assess_country (str): The country being assessed.
            weather_data (pd.DataFrame): weather data.

        Returns:
            (cop_array, elect_demand_ind, ind_heat_demand, bld_heat_demand, elect_demand_res): Energy demand data.
        """
        cop_array, elect_demand_ind, ind_heat_demand, bld_heat_demand, elect_demand_res = get_energy_demands_region(weather_data.temp_air)
        elect_demand_res.loc[:] = 0 if not self.residential_sector else elect_demand_res #electricity demand for residential households is zero if not residential sector
        bld_heat_demand.loc[:] = 0 if not self.residential_sector else bld_heat_demand # building heat demand for residential households is zero if not residential sector
        
        return cop_array, elect_demand_ind, ind_heat_demand, bld_heat_demand, elect_demand_res

    def process_data(self):
        """
        Processes data and returns a DataFrame as input for an optimization problem.

        Returns:
            df_data (pd.DataFrame): Processed data.
        """
        # Get weather data and elevation
        weather_data, elevation = self.get_weather_data()

        # Get renewable profiles
        pv_MW_array, wind_MW_array_off, wind_MW_array_on = self.get_renewable_profiles(weather_data)

        # Get electricity prices and revenue from injections
        da_prices, rev_inj = self.get_elect_prices(self.assess_country, self.get_year, fit=FIT)

        # Get energy demands
        cop_array, elect_demand_ind, ind_heat_demand, bld_heat_demand, elect_demand_res = self.get_energy_demands(weather_data)

        # Create DataFrame with data
        TPs = pd.read_excel(r"new_data\TP_BEVs.xlsx", index_col=[0])
        TP_1 = list(TPs.TP1) * 365
        TP_2 = list(TPs.TP2) * 365
        TP_3 = list(TPs.TP3) * 365
        
        wind_MW_start_on_micro = wind_profile_generator_tmy(weather_data, self.lat, self.lon, turbine_spec="E48/800", hub_height = 40)
        
        df_data = pd.DataFrame(data={'pv_MW_array': pv_MW_array, 'wind_MW_array_off': wind_MW_array_off,
                                    'wind_MW_array_on': wind_MW_array_on, 'wind_MW_array_on_micro': wind_MW_start_on_micro,
                                    "temperature": weather_data.temp_air, "irradiance": weather_data['ghi'] / 1000,
                                    "ind_heat_demand": np.array(ind_heat_demand), "bld_heat_demand": np.array(bld_heat_demand),
                                    "elect_demand_ind": np.array(elect_demand_ind),
                                    "cop_array": cop_array, "electricity_demand": np.array(
                                        elect_demand_res.values + elect_demand_ind.values),
                                    "elect_demand_res": np.array(elect_demand_res), "grid_abs_price": da_prices,
                                    "rev_inj": rev_inj, "TP_1": TP_1, "TP_2": TP_2, "TP_3": TP_3,
                                    }, index=pd.date_range('1/1/{} 00:00'.format(ASSESSMENT_YEAR), periods=8760, freq='H'))

        dict_limits = get_max_caps_regions(constrained= self.constrained, residential_sector = self.residential_sector)

        return df_data, dict_limits

class ReferenceEnergySystemCapacities:
    """
    Class for simulating a reference energy system.
    """

    def __init__(self, df_data, parm, dict_limits, total_hh,
                 share_res_oil_boiler, share_res_hp, share_res_e_heating, 
                 share_res_ng_boiler, share_ind_boiler):
        """
        Initialize a new instance of ReferenceEnergySystemCapacities.

        Parameters:
            df_data (DataFrame): DataFrame containing energy demand data.
            parm (dict): Dictionary of parameters for the energy system.
            dict_limits (dict): Dictionary of limits for the energy system.
            total_hh (int): Total number of households in the system.
            share_res_oil_boiler (float): Share of residential heating from oil boilers.
            share_res_hp (float): Share of residential heating from heat pumps.
            share_res_e_heating (float): Share of residential heating from electric heating.
            share_res_ng_boiler (float): Share of residential heating from natural gas boilers.
            share_ind_boiler (float): Share of industrial heating from boilers.
        """
        # Assigning input values to instance variables
        self.df_data = df_data
        self.parm = parm
        self.dict_limits = dict_limits
        self.total_hh = total_hh
        
        # Share of energy sources
        self.share_res_oil_boiler = share_res_oil_boiler
        self.share_res_hp = share_res_hp
        self.share_res_e_heating = share_res_e_heating
        self.share_res_ng_boiler = share_res_ng_boiler
        self.share_ind_boiler = share_ind_boiler
        
        if (share_res_oil_boiler+share_res_hp+share_res_e_heating+share_res_ng_boiler) % 1 > 0:
            print("WARNING, sum is not equal to one or zero: '{}'".format((share_res_oil_boiler+share_res_hp+share_res_e_heating+share_res_ng_boiler)))

    def determine_capacities(self):
        """
        Calculate various capacities related to the energy system.
        """
        # Calculate grid connection capacity
        self.p_grid_connection = (
                self.df_data.electricity_demand.max() +
                (self.df_data.bld_heat_demand * self.share_res_hp / self.df_data['cop_array']).max() +
                (self.share_res_e_heating * self.df_data.bld_heat_demand / self.parm['e_heating_eff']).max()
            ) if self.share_res_e_heating > 0 else self.df_data.electricity_demand.max() + (self.df_data.bld_heat_demand * self.share_res_hp / self.df_data['cop_array']).max()

        # Calculate industrial boiler capacity, here assumed that the entire industrial heat demand is met by this boiler
        self.cap_d_boiler_ind = self.df_data.ind_heat_demand.max() / self.parm['d_boiler_eff']

    def calculate_residential_heating_capacities(self):
        """
        Calculate capacities related to residential heating sources.
        """
        # Calculate capacities of different residential heating sources
        self.cap_boiler_oil_res = self.share_res_oil_boiler * self.df_data.bld_heat_demand.max() / self.parm['o_boiler_eff']
        self.cap_hp_res = self.share_res_hp * self.df_data.bld_heat_demand.max() / AV_COP_ASHP # power capacity of HP
        self.cap_e_heating_res = self.share_res_e_heating * self.df_data.bld_heat_demand.max() / self.parm['e_heating_eff']
        self.cap_boiler_ng_res = self.share_res_ng_boiler * self.df_data.bld_heat_demand.max() / self.parm['gb_eff_res']

        # Calculate energy and fuel consumption for oil boiler in residential heating
        self.p_boiler_oil_res = self.share_res_oil_boiler * self.df_data.bld_heat_demand.sum()
        self.f_boiler_oil_res = self.share_res_oil_boiler * self.df_data.bld_heat_demand.sum() / self.parm['o_boiler_eff']

        # Calculate energy and fuel consumption for natural gas boiler in residential heating
        self.p_boiler_gb_res = self.share_res_ng_boiler * self.df_data.bld_heat_demand.sum()
        self.f_boiler_gb_res = self.share_res_ng_boiler * self.df_data.bld_heat_demand.sum() / self.parm['gb_eff_res']

    def calculate_grid_absorption(self):
        """
        Calculate the absorption of electricity from the grid.
        """
        # Calculate overall grid absorption
        self.p_grid_abs_overall = self.df_data.electricity_demand.sum()

        # Calculate grid absorption from heat pumps
        self.p_grid_abs_hp = (self.df_data.bld_heat_demand * self.share_res_hp / self.df_data['cop_array']).sum()

        # Calculate grid absorption from electric heating
        self.p_grid_abs_eh = (self.df_data.bld_heat_demand * self.share_res_e_heating / self.parm['e_heating_eff']).sum()

        # Calculate total grid absorption
        self.p_grid_abs = self.p_grid_abs_overall + self.p_grid_abs_hp + self.p_grid_abs_eh

    def calculate_industrial_boiler(self):
        """
        Calculate industrial boiler parameters.
        """
        # Calculate industrial boiler power and fuel consumption
        self.p_d_boiler = self.df_data.ind_heat_demand.sum()
        self.f_d_boiler = self.df_data.ind_heat_demand.sum() / self.parm['d_boiler_eff']

    def calculate_transportation(self):
        """
        Calculate parameters related to transportation.
        """
        # Calculate share of gasoline vehicles
        self.share_gasoline = 1 - self.dict_limits['share_bev']

        # Calculate kilometers traveled by petrol and electric vehicles
        self.petrol_car_kms = self.share_gasoline * 365 * self.parm['km_per_day'] * self.total_hh
        self.bev_car_kms = self.dict_limits['share_bev'] * 365 * self.parm['km_per_day'] * self.total_hh

        # Calculate additional grid absorption due to electric vehicle usage
        self.additional_grid_abs = self.bev_car_kms * self.parm['bev_MWh_km']

        # Update total grid absorption with additional electric vehicle demand
        self.p_grid_abs += self.additional_grid_abs

def calc_cop(Tamb: float, Tsupp: float) -> float:
    """
    Determines the coefficient of performance from an example HP from data in Fischer et al. (2017)
    https://www.sciencedirect.com/science/article/pii/S0360544216315572 

    Args:
        Tamb (float): ambient temperature [degrees C].
        Tsupp (float): heat supply temperature [degrees C].
    Returns:
        float: coefficient of performance.
    """
    T_lift = Tsupp - Tamb
    z_1 = 5.06 - 0.05 * T_lift
    z_2 = 0.00006 * T_lift**2
    return z_1 + z_2

def get_roughness_length(latitude, longitude, overseas=False):
    """
    Determines roughness length, now simply assumed to be 0.05 for overseases and 0.15 onshore. A tif file can be used instead...
    Args:
        latitude (float): latitude [degrees].
        longitude (float): longitude [degrees].
        overseas (bool): whether location is offshore.
    Returns:
        float: rouhgness length.
    """
    if overseas:
        print("WARNING: roughness length set to 0.05")
        rix = 0.05
    else:
        print("WARNING: roughness length set to 0.15")
        rix = 0.15
    return rix

def pv_profile_generator_tmy(weather_data, latitude: float, longitude: float,
                             module_spec="Sharp_NDQ235F4__2013_", type_pv = "open_rack"):
    """
    Calculates PV generation per kWp.

    Args:
        weather_data (dataframe): dataframe with weather information for one year [pd.Dataframe].
        latitude (float): latitude, in decimal degrees, between -90 and 90, north is positive (ISO 19115) [degrees].
        longitude (float): longitude, in decimal degrees, between -180 and 180, east is positive (ISO 19115) [degrees].
        module_spec (str): type of solar module to be modelled [str] default is Sharp_NDQ235F4__2013_. The PV module modelled, here standard the Sharp, since parameters are available and multi-Si
        type_pv (str): str, default is roof_mounted. How the module will be installed, here assumed to be roof mounted. However, can be modified in the future
    Returns:
        total_profile_tmy (dataframe): the PV profile per kWp panel modelled, for a TMY.
    """
    
    # Assuming that panels are oriented to the south with angle of 35 deg to surface
    sa = 180 #Surface azimuth angle, 180 (south), 270(west), 90(east)
    st = 35 #angle the roof makes with the surface, assumed to be 35 degrees
    
    # If TMY, delete first substrings from rows to avoid inconsitencies
    temperature = weather_data['temp_air'].astype(float)
    windspeed = weather_data['wind_speed'].astype(float)
    ghi = weather_data['ghi']
    dni = weather_data['dni']
    dhi = weather_data['dhi']
    pressure = weather_data['pressure'].astype(float)
    
    #solar details for location
    sp = pvlib.solarposition.get_solarposition(weather_data.index, float(latitude), float(longitude), pressure=pressure)
    irradiance = pvlib.irradiance.get_total_irradiance(st, sa, sp.apparent_zenith, sp.azimuth, dni, ghi, dhi)
    
    # Get db with modules, choose the one in the function and get specs
    modules = pvsystem.retrieve_sam('SandiaMod').T
    modules_multi_si = modules[modules['Material']=="mc-Si"]
    multi_si = modules_multi_si[modules_multi_si.index.str.contains(module_spec)].T
    multi_si = multi_si[module_spec]
    
    #relative airmass and aoi
    airmass = pvlib.atmosphere.get_relative_airmass(sp.apparent_zenith)
    airmass_absolute = pvlib.atmosphere.get_absolute_airmass(airmass, pressure=pressure)

    #aoi and celltemp
    aoi = pvlib.irradiance.aoi(st, sa, sp.zenith, sp.azimuth)

    if type_pv == "roof_mounted":
        temperature_model_parameters = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS['sapm']['close_mount_glass_glass']
        celltemp = pvlib.temperature.sapm_cell(irradiance.poa_global, temperature,windspeed, **temperature_model_parameters)
    elif type_pv == "open_rack":
        temperature_model_parameters = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS['sapm']['open_rack_glass_polymer']
        celltemp = pvlib.temperature.sapm_cell(irradiance.poa_global, temperature,windspeed, **temperature_model_parameters)
    else:
        raise ValueError("Invalid type of PV panel installation")
        
    effective_irradiance = pvlib.pvsystem.sapm_effective_irradiance(irradiance.poa_direct, irradiance.poa_diffuse, 
                                                                    airmass_absolute, aoi, multi_si)
    
    #Multi-Si
    dc = pvlib.pvsystem.sapm(effective_irradiance, celltemp, multi_si)

    #Multi-Si specifications
    Imp_multiSi = multi_si.Impo #8.1243
    Vmp_multiSi = multi_si.Vmpo #29.1988
    kWp = (Imp_multiSi * Vmp_multiSi) / 1000
    
    #use pv Watts model to convert dc to ac power (inverter model)
    sapm_inverters = pvlib.pvsystem.retrieve_sam('cecinverter')
    inverter = sapm_inverters['ABB__MICRO_0_25_I_OUTD_US_208__208V_']
    ac = pvlib.inverter.sandia(dc['v_mp'], dc['p_mp'], inverter)

    #convert to kW and fill NaN values
    ac_mod = ac.fillna(0) / 1000
    ac_mod[ac_mod < 0] = 0

    # Generate profile per 1 kWp
    total_profile_tmy = (1/kWp) * ac_mod
    
    print("Calculated PV profile with a capacity factor of '{}'".format(round(total_profile_tmy.mean(),4)))
    
    return total_profile_tmy

def wind_profile_generator_tmy(weather_data, latitude: float, longitude: float, turbine_spec="E-126/4200", hub_height = 135, 
                               roughness_length = 0.15, offshore=False):
    """
    Calculates wind generation per kWp.

    Args:
        weather_data (dataframe): dataframe with weather information for one year [pd.Dataframe].
        latitude (float): latitude, in decimal degrees, between -90 and 90, north is positive (ISO 19115) [degrees].
        longitude (float): longitude, in decimal degrees, between -180 and 180, east is positive (ISO 19115) [degrees].
        turbine_spec (str): type of wind turbine to be modelled [str] default is E-126/4200.
        hub_height (float): float, default is 135 meter.
        roughness_length (float): standard roughtness length applied.
        offshore (bool): if True, the location is offshore.
    Returns:
        total_profile (dataframe): the wind profile per kWp turbine modelled, for a TMY.
    """
    
    df = pd.DataFrame(data=[], 
                  index = pd.date_range('1/1/{} 00:00'.format(ASSESSMENT_YEAR), periods=8760, freq='H'), 
                 columns=np.arange(0,5))
    df.columns = [
        ['pressure','temperature','wind_speed','roughness_length','temperature'],
        [0, 2, 10, 0, 10]]
    df.columns.names = ['variable_name', 'height']

    """
        DataFrame with time series for wind speed `wind_speed` in m/s,
        temperature `temperature` in K, roughness length `roughness_length`
        in m, and pressure `pressure` in Pa.
        The columns of the DataFrame are a MultiIndex where the first level
        contains the variable name (e.g. wind_speed) and the second level
        contains the height at which it applies (e.g. 10, if it was
        measured at a height of 10 m).
    """
    # Replace new data in df
    df['temperature',2] = weather_data.temp_air + 273.15
    df['temperature',10] = wp.temperature.linear_gradient(weather_data['temp_air']+ 273.15, 2, 10)
    df['wind_speed',10] = weather_data.wind_speed
    df['roughness_length', 0] = get_roughness_length(latitude, longitude, overseas=offshore)
    df['pressure', 0] = weather_data.pressure
    
    # Forward fill in case empty or NaN values
    df.ffill(inplace=True)

    # specification of wind turbine where power curve is provided in the
    # oedb turbine library
    enercon_e126 = {
            'turbine_type': turbine_spec,  # turbine type as in oedb turbine library
            'hub_height': hub_height  # in m
        }
    # initialize WindTurbine object
    e126 = WindTurbine(**enercon_e126)
    
    # own specifications for ModelChain setup
    modelchain_data = {
        'wind_speed_model': 'logarithmic',      # 'logarithmic' (default),
                                                # 'hellman' or
                                                # 'interpolation_extrapolation'
        'density_model': 'ideal_gas',           # 'barometric' (default), 'ideal_gas'
                                                #  or 'interpolation_extrapolation'
        'temperature_model': 'linear_gradient', # 'linear_gradient' (def.) or
                                                # 'interpolation_extrapolation'
        'power_output_model':
            'power_curve',                      # 'power_curve' (default) or
                                                # 'power_coefficient_curve'
        'density_correction': True,             # False (default) or True
        'obstacle_height': 0,                   # default: 0
        'hellman_exp': None}                    # None (default) or None

    # initialize ModelChain with own specifications and use run_model method to
    # calculate power output
    mc_e126 = ModelChain(e126, **modelchain_data).run_model(df)
    
    # write power output time series to WindTurbine object
    e126.power_output = mc_e126.power_output
    
    # Calculcate power output per kW, and modify export file
    total_profile = e126.power_output/e126.nominal_power
    
    # Sometimes it can happen that the pwoer output is slightly higher than the nomial power, i.e. make sure that this is not allowed;
    total_profile[total_profile > 1] = 1 

    print("Calculated wind profile with a capacity factor of '{}'".format(round(total_profile.mean(),3)))
    return total_profile

def get_max_caps_regions(st_max = MAX_CAP_TECHS, st_max_res = MAX_CAP_TECHS_RES,
                         min_cap=0, min_cap_res=0, av_p_waste_capita = 5,
                         residential_sector = True,
                         constrained=False):
    
    """
    Get maximum capacity constraints for a specified region.

    Parameters:
        st_max (int): Maximum capacity for technologies (in MW).
        st_max_res (int): Maximum capacity for residential technologies (in MW).
        min_cap (int): Minimum capacity for certain technologies (in MW).
        min_cap_res (int): Minimum capacity for residential technologies (in MW).
        av_p_waste_capita (int): Average waste production per capita (in W).
        residential_sector (boolean): Include residential sector. Default is True.
        mb (bool): If True, applies strict constraints for Manna Bakery.

    Returns:
        dict: A dictionary containing the maximum capacity constraints for various technologies.
    """
    
    if constrained:
        dict_reg_limits = {
            "max_pv": 0.5,
            "max_wind_on": 0.06,
            "max_wind_off": 0,
            "max_bm_annu": BIOMASS_ENERGY,
            "max_grid_cap": 2,
            "st_max": st_max,
            "share_bev": SHARE_BEV_INIT,
            "max_h2_export": MAX_H2_EXPORT,
            "biogas_potential_MWh":RESIDENTS * 8760 * av_p_waste_capita / 1e6 if residential_sector else 0, # to MW
            "max_solar_heat_ha": st_max,#m2
            "min_solar_heat_ha": (1/10000), #1 m2
            "min_cap_wind_off": min_cap, 
            'min_cap': min_cap,  
            'min_cap_ind': MIN_CAP_IND,
            'min_cap_res': min_cap_res,
            'max_bat': st_max*10,
            'max_h2_storage': st_max*10,
            'max_electrolyzer': st_max,
            'min_electrolyzer': MIN_ELECTROLYZER_CAP,
            'max_fc': st_max,
            'max_hp': st_max_res,
            'max_boiler': st_max,
            'max_boiler_res': st_max_res,
            'max_heat_storage': HEAT_STORAGE_RES, 
            'max_chp_mix': 5,
            'min_chp_mix': min_cap,
            'max_gasif': st_max, 
            'households': np.ceil(RESIDENTS / COST_DATA[NAME_REF_DB].res_home) if residential_sector else 0}
    else:
        dict_reg_limits = {
            "max_pv": st_max,
            "max_wind_on": st_max,
            "max_wind_off": st_max,
            "max_bm_annu": BIOMASS_ENERGY, 
            "max_grid_cap": MAX_GRID_CAP,
            "st_max": st_max,
            "share_bev": SHARE_BEV_INIT,
            "max_h2_export": MAX_H2_EXPORT,
            "biogas_potential_MWh":RESIDENTS * 8760 * av_p_waste_capita / 1e6 if residential_sector else 0,
            "max_solar_heat_ha": st_max,#m2
            "min_solar_heat_ha": (1/10000), #1 m2
            "min_cap_wind_off": min_cap, 
            'min_cap': min_cap,    
            'min_cap_ind': MIN_CAP_IND,
            'min_cap_res': min_cap_res,
            'max_bat': st_max*10,
            'max_h2_storage': st_max*10,
            'max_electrolyzer': st_max,
            'min_electrolyzer': MIN_ELECTROLYZER_CAP,
            'max_fc': st_max,
            'max_hp': st_max_res,
            'max_boiler': st_max,
            'max_boiler_res': st_max_res,
            'max_heat_storage': HEAT_STORAGE_RES,
            'max_chp_mix': st_max,
            'min_chp_mix': min_cap,
            'max_gasif': st_max, 
            'households': np.ceil(RESIDENTS / COST_DATA[NAME_REF_DB].res_home) if residential_sector else 0
        }
    
    return dict_reg_limits

def get_energy_demands_region(temperature, sec_db = NAME_REF_DB,
                              file_name_ind = FILE_NAME_ENERGY_IND, 
                              file_elect_demand_res = FILE_ELECT_DEMAND_RES):
    """
    Calculates energy demands for a specific location and temperature.

    Args:
        temperature (np.array): The temperature at the location.
        sec_db (str, optional): Name of the reference database. Default is NAME_REF_DB.
        file_name_elect (str, optional): Name of file with power demand of industry. Default is FILE_NAME_ELECT.        
        file_elect_demand_res (str, optional): Name of file with power demand of residential households. Default is FILE_ELECT_DEMAND_RES.

    Returns:
        cop_array (numpy.array): Coefficient of Performance array based on temperature.
        elect_demand_ind (numpy.array): Individual electricity demand array.
        ind_heat_demand (numpy.array): Individual heat demand array.
        bld_heat_demand (numpy.array): Building heat demand array.
        elect_demand_res (numpy.array): Residential electricity demand array.
    """
    
    cop_array = np.array(calc_cop(temperature, COST_DATA[sec_db].hp_supp_t))
    
    all_bld_heat_loss = TOTAL_HH * COST_DATA[sec_db].heat_loss_buildings
    temperature_eff = COST_DATA[sec_db].heat_threshold - temperature
    temperature_eff[temperature_eff<0] = 0
    bld_heat_demand = temperature_eff * all_bld_heat_loss * COST_DATA[sec_db].transmission_loss / 1e3 #from kW to MW

    # Data from industry
    data = pd.read_excel(file_name_ind, sheet_name = "Hourly Data", index_col = "DATETIME") 
    
    # electricity demand industry
    elect_demand_ind = data['Electricity Consumption (kWh)'][:8760] / 1e3
    elect_demand_ind[elect_demand_ind<0] = 0
    elect_demand_ind[elect_demand_ind<0.0001] = np.nan
    elect_demand_ind.fillna(method='ffill',inplace=True)
    
    # heat demand industry   
    ind_heat_demand = data['Heat Consumption (kWh)'][:8760] / 1e3
    ind_heat_demand[ind_heat_demand<0] = 0
    ind_heat_demand.fillna(method='ffill', inplace=True)

    # electricity demand residents
    elect_demand_res = pd.read_csv(file_elect_demand_res, index_col = [0])
    elect_demand_res = elect_demand_res.iloc[:,0]
            
    return cop_array, elect_demand_ind, ind_heat_demand, bld_heat_demand, elect_demand_res
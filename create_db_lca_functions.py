import brightway2 as bw
from brightway2 import *
import bw2io
from premise import *
from platform import python_version
import premise
from bw_recipe_2016 import WaterConsumption, get_biosphere_database
from bw2data import Database, Method
from premise_gwp import add_premise_gwp

from config import (NAME_REF_DB, ASSESS_COUNTRY, DB_NAME_INIT, EI_VERSION,
                    PROJECT_NAME, NAME_REF_DB, KEY_PREMISE, ASSESSMENT_YEAR)

print("Using Python v.({}) and premise v.{}.".format(python_version(), str(premise.__version__).replace(", ", ".")))

def import_additional_lcias():
    """
    Import additional LCIA methods related to water consumption and land transformation.

    This function imports and applies LCIA methods for water consumption and creates a new environmental impact
    category for land transformation. It utilizes functions from the 'premise_gwp' and 'bw_recipe_2016' packages.

    Steps:
    1. Adds premise global warming potential (GWP) methods using 'add_premise_gwp'.
    2. Retrieves the biosphere database using 'get_biosphere_database'.
    3. Creates and applies LCIA methods for water consumption.
    4. Defines a new `LCIA' method for land transformation and writes impact factors.

    Returns:
    None

    Example:
    >>> import_additional_lcias()
    """
    # define project
    bw.projects.set_current(PROJECT_NAME) #Creating/accessing the project
    bw.bw2setup() #Importing elementary flows, LCIA methods and some other data

    # Step 1: Add premise global warming potential (GWP) methods
    add_premise_gwp()

    # Step 2: Retrieve the biosphere database
    biosphere = get_biosphere_database()

    # Step 3: Create and apply LCIA methods for water consumption
    gw = WaterConsumption(None, biosphere)
    gw.apply_strategies()
    gw.write_methods(overwrite=True)
    gw.data[0]

    # Step 4: Create a new LCIA method for land transformation
    my_cfs_land = []

    for bio in Database("biosphere3"):
        if "Transformation, from" in bio['name'] and "square meter" == bio['unit']:   
            line = (bio.key, 1)
            my_cfs_land.append(line)

    my_method = Method(("Own method", "Land transformation", "Land transformation"))
    my_metadata = {"unit": "m2", "meaning": "to represent land transformation"}
    my_method.register(**my_metadata)
    my_method.write(my_cfs_land)

def import_ecoinvent_database(db_name, location_path, overwrite=True):
    """
    Import an Ecoinvent database and create default LCIA methods if not already imported.

    This function imports an Ecoinvent database using the specified database name and location path.
    It checks whether the database is already imported and, if not, imports it, applies strategies,
    provides statistics, and writes the database. Additionally, it creates default LCIA methods and core migrations.

    Parameters:
    - db_name (str): Database name for the Ecoinvent dataset.
    - location_path (str): Location path to the datasets subfolder of the unzipped Ecoinvent file.
    - overwrite (bool, optional): If True, overwrites existing LCIA methods. Default is False.

    Returns:
    None

    Example:
    >>> import_ecoinvent_database("ecoinvent_init", "/path/to/ecoinvent/datasets", overwrite=True)
    """
    # Check if the database is already imported
    if db_name in bw.databases:
        print(f"{db_name} has already been imported.")
    else:
        # Import and process the Ecoinvent database
        ei_importer = bw.SingleOutputEcospold2Importer(location_path, db_name)
        ei_importer.apply_strategies()
        ei_importer.statistics()
        ei_importer.write_database()
        bw2io.create_default_lcia_methods(overwrite=overwrite)
        bw2io.create_core_migrations()

def generate_future_ei_dbs(scenarios = ["SSP2-Base", "SSP2-PkBudg1150","SSP2-PkBudg500"], iam = 'remind',
                           start_yr=2020, end_yr = 2050, step = 15, endstring="base"):
    """
    Generate Ecoinvent scenario models with specified parameters.

    This function generates Ecoinvent scenario models based on the specified year, scenario, and additional settings.
    It avoids adding duplicated databases by checking the existing databases in Brightway2.

    Parameters:
    - scenarios (list): The scenarios for which the models are generated. Default is: 
                ["SSP2-Base",
                "SSP2-PkBudg1150",
                "SSP2-PkBudg500"] corresponding to baseline, 2 degrees C, and 1.5 degrees C.
    - iam (str): IAM chosen, can be 'remind' or 'image'.
    - start_yr (int): The starting year for the scenarios.
    - end_yr (int): The end year for the scenarios.
    - step (int): step between scenario years.
    - endstring (str, optional): A suffix to differentiate the generated databases. Default is "base".

    Returns:
    tuple: A tuple containing two lists -
        1. List of dictionaries specifying the models for the scenarios.
        2. List of database names generated based on the specified parameters.

    Example:
    >>> generate_future_ei_dbs("SSP2-Base")
    ([{'model': 'remind', 'pathway': 'SSP2-Base', 'year': 2030},
      {'model': 'remind', 'pathway': 'SSP2-Base', 'year': 2050}],
     ['ecoinvent_remind_SSP2-Base_2030_custom', 'ecoinvent_remind_SSP2-Base_2050_base'])
    """

    list_years = [start_yr + i * step for i in range(1, int((end_yr - start_yr) / step) + 1)]

    list_spec_scenarios = []
    list_names = []

    for pt in scenarios:
        for yr in list_years:
            string_db = "ecoinvent_{}_{}_{}_{}".format(iam, pt, yr, endstring)

            if yr == ASSESSMENT_YEAR and pt == "SSP2-Base":
                dict_spec = {"model": iam, "pathway": pt, "year": yr,
                                "exclude": ["update_electricity", "update_cement", "update_steel", "update_dac",
                                            "update_fuels", "update_emissions", "update_two_wheelers"
                                            "update_cars", "update_trucks", "update_buses"]}

                if string_db not in bw.databases:
                    list_spec_scenarios.append(dict_spec)
                    list_names.append(string_db)
                else:
                    print("Avoid duplicated db and therefore following db not added: '{}'".format(string_db))
            else:
                dict_spec = {"model": iam, "pathway": pt, "year": yr}

                if string_db not in bw.databases:
                    list_spec_scenarios.append(dict_spec)
                    list_names.append(string_db)
                else:
                    print("Avoid duplicated db and therefore following db not added: '{}'".format(string_db))

    return list_spec_scenarios, list_names

# ### generate the database which we are going to use, as premise include many novel datasets. Add some datasets that we generated ourselves.
def generate_reference_database():
    """
    Generate a reference database based on specified parameters.

    This function generates a reference database based on the specified parameters.
    It deletes the existing reference database with the same name if it exists and then creates a new one.

    Returns:
    None

    Example:
    >>> generate_reference_database()
    """
    # Delete old reference database with the same name
    for db_name in list(bw.databases):
        if NAME_REF_DB in db_name:
            print(db_name)
            del bw.databases[db_name]

    # Create a new reference database using NewDatabase
    ndb = NewDatabase(
        scenarios=[{"model": "remind", "pathway": 'SSP2-Base', "year": "2022",
                    "exclude": ["update_electricity", "update_cement", "update_steel", "update_dac",
                                "update_fuels", "update_emissions"]}],
        source_db=DB_NAME_INIT,
        source_version=EI_VERSION,
        key=KEY_PREMISE,
        additional_inventories=[
            {"filepath": r"new_data\lci-NMC-battery-import.xlsx", "ecoinvent version": EI_VERSION},
            {"filepath": r"new_data\CONFIDENTIAL_AURELIA.xlsx", "ecoinvent version": EI_VERSION},
            {"filepath": r"new_data\lci-add.xlsx", "ecoinvent version": EI_VERSION},
            {"filepath": r"new_data\lci-BOS-import.xlsx", "ecoinvent version": EI_VERSION}]
    )

    # Write the new reference database to Brightway2
    ndb.write_db_to_brightway(name=NAME_REF_DB)

def generate_specialized_datasets():
    """
    Generate specialized datasets for specific activities and locations based on specified parameters.

    This function generates specialized datasets for heat production from natural gas, solar thermal, wood production,
    and electric vehicles based on the specified parameters. It includes logic to avoid duplicates and modify
    existing datasets to exclude specific infrastructure.

    Returns:
    None

    Example:
    >>> generate_specialized_datasets()
    """

    for sec_db in list(bw.databases):
        if "ecoinvent_" in str(sec_db):
            print(sec_db)
            check_act = [act for act in bw.Database(sec_db) if 'heat production, natural gas, w\o infrastructure' in act['name']
                       and "Europe without Switzerland" in act['location']]#[0]

            if len(check_act) > 0:
                #check_act.delete()
                print("Activity 'heat production, natural gas, w\o infrastructure' already generated")
            else:
                # Get activity which we want to modify
                spec_act = [act for act in bw.Database(sec_db) if 
                            'heat production, natural gas, at boiler condensing modulating <100kW' == act['name']
                       and "Europe without Switzerland" == act['location']][0]
                # Copy the process
                spec_act = spec_act.copy()
                # Set a new name and reference product to the copied activity
                spec_act['name'] = "heat production, natural gas, w\o infrastructure"
                spec_act['reference product'] = "heat production, natural gas, w\o infrastructure"

                # Delete exchanges which have to do with infrastructure
                for exc in spec_act.exchanges():
                    #print(exc)
                    if ('market for gas boiler' == exc['name']) or ('gas boiler production' == exc['name']):
                        exc.delete()
                        print("'gas boiler' deleted")

                spec_act.save()

    # ### Make a new dataset for solar thermal and exclude the infrastructure that should be not accounted for.
    for sec_db in list(bw.databases):
        if "ecoinvent_" in str(sec_db):
            print(sec_db)

            check_act = [act for act in bw.Database(sec_db) if 'solar collector system installation, Cu flat plate collector, one-family house, combined system, w\o storage' == act['name']
                       and "CH" == act['location']]#[0]

            if len(check_act) > 0:
                #check_act.delete()
                print("Activity 'solar collector system installation, Cu flat plate collector, one-family house, combined system, w\o storage' already generated")
            else:
                # Get activity which we want to modify
                spec_act = [act for act in bw.Database(sec_db) if 
                            'solar collector system installation, Cu flat plate collector, one-family house, combined system' == act['name']
                       and "CH" == act['location']][0]
                # Copy the process
                spec_act = spec_act.copy()
                # Set a new name and reference product to the copied activity
                spec_act['name'] = "solar collector system installation, Cu flat plate collector, one-family house, combined system, w\o storage"

                # Delete exchanges which have to do with infrastructure
                for exc in spec_act.exchanges():
                    if 'heat storage production, 2000l' == str(exc['name']):
                        exc.delete()
                        print("'heat storage production, 2000l' deleted")

                spec_act.save()

    # ### Make a location-specific wood dataset, for the specific location
    for sec_db in list(bw.databases):
        if "ecoinvent_" in str(sec_db):
            print(sec_db)
            check_act = [act for act in bw.Database(sec_db) if 'synthetic gas production, from wood, at fixed bed gasifier, w\o infrastructure' in act['name']
                       and ASSESS_COUNTRY == act['location'] and act['reference product'] == 'synthetic gas']#[0]

            if len(check_act) > 0:
                #check_act.delete()
                print("Activity 'synthetic gas production, from wood, at fixed bed gasifier, w\o infrastructure '{}' already generated".format(ASSESS_COUNTRY))
            else:
                # Get activity which we want to modify
                spec_act = [act for act in bw.Database(sec_db) if 
                            'synthetic gas production, from wood, at fixed bed gasifier' == act['name']
                       and "CH" == act['location']][0]
                # Copy the process
                spec_act = spec_act.copy()
                # Set a new name and reference product to the copied activity
                spec_act['name'] = "synthetic gas production, from wood, at fixed bed gasifier, w\o infrastructure"
                spec_act['reference product'] = "synthetic gas"
                spec_act['location'] = ASSESS_COUNTRY
                spec_act.save()

                # Delete exchanges which have to do with infrastructure
                for exc in spec_act.exchanges():
                    print(exc)
                    if 'synthetic gas factory construction' == exc.input['name']:
                        print(exc['name'])
                        exc.delete()
                        print("'synthetic gas factory construction' deleted")
                        
                    elif 'industrial furnace production, natural gas' == exc.input['name']:
                        print(exc['name'])
                        exc.delete()
                        print("'industrial furnace production, natural gas' deleted")

                    elif 'market for wood chips, dry, measured as dry mass' == exc.input['name'] and ("RER" != exc.input['location']):
                        print(exc['name'])

                        amount = exc['amount']
                        exc.delete() 
                        select_act = [act for act in bw.Database(sec_db) if 'market for wood chips, dry, measured as dry mass' == act['name']
                                      and "RER" == act['location']][0]

                        spec_act.new_exchange(input=select_act.key, amount = amount, type='technosphere').save()     
                        print("'market for wood chips, dry, measured as dry mass' changed")

                    elif ('market for electricity, medium voltage' == exc.input['name']) and (ASSESS_COUNTRY != exc.input['location']):
                        print(exc['name'])

                        amount = exc['amount']
                        select_act = [act for act in bw.Database(sec_db) if 'market for electricity, medium voltage' == act['name']
                                      and ASSESS_COUNTRY == act['location']][0]
                        spec_act.new_exchange(input=select_act.key, amount = amount, type='technosphere').save() 
                        exc.delete()     
                        print("'market for electricity, medium voltage' changed")

                    elif 'market for tap water' == exc.input['name'] and ("Europe without Switzerland" != exc.input['location']):
                        amount = exc['amount']
                        select_act = [act for act in bw.Database(sec_db) if 'market for tap water' == act['name']
                                      and 'Europe without Switzerland' == act['location']][0]
                        spec_act.new_exchange(input=select_act.key, amount = amount, type='technosphere').save() 
                        exc.delete()     
                        print("'market for tap water' changed")

                    elif 'market for wastewater, average' == exc.input['name'] and ("Europe without Switzerland" != exc.input['location']):
                        print(exc['name'])
                        amount = exc['amount']
                        select_act = [act for act in bw.Database(sec_db) if 'market for wastewater, average' == act['name']
                                      and 'Europe without Switzerland' == act['location']][0]
                        spec_act.new_exchange(input=select_act.key, amount = amount, type='technosphere').save() 
                        exc.delete()    
                        print("'market for wastewater, average' changed")

                    elif 'market for wood ash mixture, pure' == exc.input['name'] and ("Europe without Switzerland" != exc.input['location']):
                        print(exc['name'])
                        amount = exc['amount']
                        select_act = [act for act in bw.Database(sec_db) if 'market for wood ash mixture, pure' == act['name']
                                      and 'Europe without Switzerland' == act['location']][0]
                        spec_act.new_exchange(input=select_act.key, amount = amount, type='technosphere').save() 
                        exc.delete()    
                        print("'market for wood ash mixture, pure' changed")
                        
                spec_act.save()


    # ### Make a new dataset that does not include for electricity demand in BEVs, as this is seperately determined within the optimization problem
    for sec_db in list(bw.databases):
        if "ecoinvent_" in str(sec_db):
            print(sec_db)

            check_act = [act for act in bw.Database(sec_db) if 'transport, passenger car, electric, w\o fuel' in act['name']
                       and "GLO" in act['location']]

            if len(check_act) > 0:
                #check_act.delete()
                print("Activity 'transport, passenger car, electric, w\o fuel' already generated")
            else:
                # Get activity which we want to modify
                spec_act = [act for act in bw.Database(sec_db) if 
                            'transport, passenger car, electric' == act['name']
                       and "GLO" == act['location']][0]
                # Copy the process
                spec_act = spec_act.copy()
                # Set a new name and reference product to the copied activity
                spec_act['name'] = "transport, passenger car, electric, w\o fuel"
                spec_act['reference product'] = "transport, passenger car, electric, w\o fuel"

                # Delete exchanges which have to do with infrastructure
                for exc in spec_act.exchanges():
                    if 'market group for electricity, low voltage' == str(exc['name']):
                        print(exc['name'])
                        exc.delete()
                        print("'market group for electricity, low voltage' deleted")

                spec_act.save()

def generate_prospective_lca_dbs(list_spec_scenarios, list_names):
    """
    Generate and update future LCA databases for prospective LCA.

    This function generates and updates future LCA databases based on specified scenarios if needed, and writes them to Brightway2.

    Parameters:
    - list_spec_scenarios (list): List of dictionaries specifying scenarios for new databases.
    - list_names (list): List of names specifying scenario names for new databases.

    Returns:
    None

    Example:
    >>> generate_and_update_lca_databases([{"model": "remind", "pathway": 'SSP2-Base', "year": "2035"}], ['ecoinvent_remind_SSP2-Base_2035_base'])
    """

    if len(list_spec_scenarios) > 0:
        ndb = NewDatabase(
            scenarios=list_spec_scenarios,
            source_db=DB_NAME_INIT,
            source_version=EI_VERSION,
            key=KEY_PREMISE,
            additional_inventories=[
                {"filepath": r"new_data\lci-NMC-battery-import.xlsx", "ecoinvent version": EI_VERSION},
                {"filepath": r"new_data\CONFIDENTIAL_AURELIA.xlsx", "ecoinvent version": EI_VERSION},
                {"filepath": r"new_data\lci-add.xlsx", "ecoinvent version": EI_VERSION},
                {"filepath": r"new_data\lci-BOS-import.xlsx", "ecoinvent version": EI_VERSION}]
        )

        print("START UPDATING")
        ndb.update()

        print("START WRITING")
        ndb.write_db_to_brightway(name = list_names)

# For example, to generate prospective db for the selected scenarios for 2035 and 2050
#list_spec_scenarios, list_names = generate_future_ei_dbs(scenarios = ["SSP2-Base", "SSP2-PkBudg1150","SSP2-PkBudg500"], iam = 'remind',
                           #start_yr=2020, end_yr = 2050, step = 15, endstring="base")
# generate_prospective_lca_dbs(list_spec_scenarios, list_names)
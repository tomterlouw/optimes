# ``OPTIMES: OPTImal design of Multi-Energy Systems.``
``optimes`` is a  repository that allows for optimizing multi-energy systems, code which has been developed during H2020 project ROBINSON (see Acknowledgements).
 
 This code can be used to optimize a (decentralized) multi-energy system (MES) considering different energy sectors: residential households, local industry, and personal mobility. To achieve this, a mixed integer linear program is formulated that allows to optimally designs the MES on costs and/or a life cycle environmental impact category. The full formulation and explanation are provided in the scientific references below.

 Here, the problem is applied to the Kissamos region in Crete (Greece), but an (almost) identical problem has been applied to Eigerøy (Norway) during H2020 ROBINSON (https://www.robinson-h2020.eu/). It is worth noting that confidential data for the local industry (bakery industry) and the life cycle inventory for the advanced CHP unit have been removed from this repository; they should be replaced with own (or dummy) data.

<div style="text-align:center">
<img src="https://ars.els-cdn.com/content/image/1-s2.0-S0306261923007389-ga1_lrg.jpg" height="450"/>
</div>
Figure 1. Example of a MES system, in this figure, applied to Eigerøy (Norway). Figure is obtained from Terlouw et al. (2023), i.e. the first reference in `Citing and scientific references`.

The ``main.py`` file has been used to generate all results for the case study in Crete. The set-up of the problem and the optimization problem can be found in ``energy_data_processor.py`` and ``opt_lca_functions.py``, respectively. The ``config.py`` file can be used to set main parameters and to connect to local files needed for the optimization.

Dependencies
====================

To run this repository or to use the code one needs:
- ``API_KEY_ENTSOE``: A key from ENTSOE to retrieve day-ahead electricity prices: https://transparency.entsoe.eu/content/static_content/Static%20content/web%20api/Guide.html. 
- ``KEY_PREMISE``: A key for ``premise`` to generate (prospective) LCA databases, including additional life cycle inventories, see https://github.com/polca/premise. Alternatively, one could focus on costs and GHG emissions and set ``CALC_ALL_LCA_IMPACTS`` to ``False`` to still allow for an optimization without calculating other LCA impact categories.
- Thus, if considering other environmental burdens in the optimization problem, one also needs access and a license for the ecoinvent database to generate the background LCA database and retrieve LCA impact factors:
https://ecoinvent.org/database/.
- A local license key for Gurobi to optimize the MES, see: https://support.gurobi.com/hc/en-us/articles/12872879801105-How-do-I-retrieve-and-set-up-a-Gurobi-license 
- Additional Python package and dependencies are listed in the file ``env_optimes.yml``.

Citing and scientific references
====================

Please cite (one of) the following papers if using this repository, data, and/or (some of) the code:
- Terlouw, T., Savvakis, N. Bauer, C., McKenna, R., Arampatzis, G. (2024). Designing multi-energy systems in Mediterranean regions towards energy autonomy. Applied Energy, 377, Volume B, 124458, https://doi.org/10.1016/j.apenergy.2024.124458. 
- Terlouw, T., Gabrielli, P., AlSkaif, T., Bauer, C., McKenna, R., & Mazzotti, M. (2023). Optimal economic and environmental design of multi-energy systems. Applied Energy, 347, 121374. https://doi.org/10.1016/j.apenergy.2023.121374. 

For the license, please refer to the ``LICENSE`` file.

Contributing
====================

Your contribution is welcome, of course! Please discuss collaborations and larger changes by contacting Tom Terlouw via tom.terlouw@psi.ch.

Acknowledgements
====================

The code is based on multiple scientific works (see ``Citing and scientific references``) performed during H2020 ROBINSON project:
https://www.robinson-h2020.eu/.
This work has received funding from the European Union's Horizon 2020 Research and Innovation Programme under Grant Agreement N. 957752 (ROBINSON).
# ``OPTIMES: OPTImal design of Multi-Energy Systems.``
Code and repository for the optimization of multi-energy systems.
 
 This code can be used to optimize a (decentralized) multi-energy system (MES) considering different energy sectors: residential households, local industry, and personal mobility. To achieve this, a mixed integer linear program is formulated that optimally designs the MES on costs and/or a life cycle environmental impact category. The full formulation and explanation are provided in the scientific references below.

 Here, the problem is applied to the Kissamos region in Crete (Greece), but an (almost) identical problem has been applied to Eigerøy (Norway) during H2020 ROBINSON (https://www.robinson-h2020.eu/). It is worth noting that confidential data for the local industry (bakery industry) and the life cycle inventory for the advanced CHP unit are removed from this repository; they should be exluded or replaced with own or dummy data.

<div style="text-align:center">
<img src="https://ars.els-cdn.com/content/image/1-s2.0-S0306261923007389-ga1_lrg.jpg" height="450"/>
</div>
Figure 1. Example of a MES system, in this figure, applied to Eigerøy (Norway). Figure is obtained from Terlouw et al. (2023), i.e. the first reference in `Citing and scientific references'.

The ``main.py`` file has been used to generate all results for the case study in Crete. The set-up of the problem and the optimization problem can be found in ``energy_data_processor.py`` and ``opt_lca_functions.py``, respectively. The ``config.py`` file can be used to set main parameters and to connect to local files needed for the optimization.

Dependencies
====================

To run this repository or to use the code one needs:
- ``API_KEY_ENTSOE``: You need a key from ENTSOE to retrieve day-ahead electricity prices: https://transparency.entsoe.eu/content/static_content/Static%20content/web%20api/Guide.html. 
- ``KEY_PREMISE``: You need a key for ``premise`` to generate (prospective) LCA databases, including additional life cycle inventories, see https://github.com/polca/premise. 
- You will also need a local license key for Gurobi, see: https://support.gurobi.com/hc/en-us/articles/12872879801105-How-do-I-retrieve-and-set-up-a-Gurobi-license 
- If considering environmental burdens in the optimization problem, one also needs access and a license for the ecoinvent database to access the LCA database:
https://ecoinvent.org/database/.
- Python package and dependencies are listed in the file ``env_optimess.yml``.

Citing and scientific references
====================

Please cite (one of) the following papers if using this repository, data, and/or (some of) the code:
- Terlouw, T., Gabrielli, P., AlSkaif, T., Bauer, C., McKenna, R., & Mazzotti, M. (2023). Optimal economic and environmental design of multi-energy systems. Applied Energy, 347, 121374. https://doi.org/10.1016/j.apenergy.2023.121374  
- Terlouw, T., Savvakis, N. Bauer, C., McKenna, R., Arampatzis, G. (2024). Designing multi-energy systems in Mediterranean regions towards energy autonomy. Applied Energy [in review].

For the license, please see the ``LICENSE`` file.

Contributing/improvements/collaborations
====================

Your contribution is welcome, of course! Please follow the pull request workflow, even for minor changes. Please discuss larger changes first by contacting Tom Terlouw via tom.terlouw@psi.ch.

Acknowledgements
====================

The code is based on multiple scientific works performed during H2020 ROBINSON project:
https://www.robinson-h2020.eu/.
This work has received funding from the European Union's Horizon 2020 Research and Innovation Programme under Grant Agreement N. 957752 (ROBINSON).

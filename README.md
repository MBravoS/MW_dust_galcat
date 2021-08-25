Code for the Milky Way dust and Large Scale Structure reconstruction from galaxy surveys described in [Bravo et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021arXiv210608818B/abstract).

The four `run_surveyname.sh` are the (PBS) submission scripts used to run the code in the Geryon HPC located in Universidad Católica.
`run_lsst_deep.sh` was used to generate the results in the main text, while `run_lsst.sh` and `run_desi.sh` for the results in the appendix.
The code is designed with some degree of modularity, using this code for different surveys would entail just adding the code to generate the appropiate input to `data_func.py` (which contain functions for parsing catalogues into the structure expected by the main code).
The reconstruction code itself is structured into three layers: `MW_dust_ext_calc.py` is the main execution level, `main_func.py` contains the top-level functions to prepare and execute the reconstruction, and `axu_func.py` contains low-level functions.
All the figures in the paper were created with `Plots_main.ipynb` and `lc_description_plots.py`, the former being executed locally and the latter in the Geryon HPC (submitted with `lc_description_plots.sh`).

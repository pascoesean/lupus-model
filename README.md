# lupus-model
### Mathematical Model of Lupus Nephritis

This experiment sets out to create a mathematical model of kidney damage flare ups in ongoing Lupus Nephritis (LN), the condition of kidney damage caused by SLE. Modeling attempts for disease dynamics during SLE in the literature are sparse, and often either too complex to be practical to implement or too broad to capture complex immune system dynamics. To address this, we have reimplemented a relatively simple model introduced by [Budu-Grajdeanu et al](https://doi.org/10.1186/1742-4682-7-14), while adding a stochastic mini-flare element to our dynamical system to model the intricate fluctuations of immune dynamics in SLE without losing generalizability to patient data. In this report, we demonstrate that the general dynamics of the system are reiterated in our model with stochastic flare, then fit these models to new patient data from the ACCESS clinical trial to reveal population level variation in model parameters, using our stochastically driven model to understand patient-to-patient heterogeneity in SLE. 


All code needed to recreate model fits and parameter manipulations is contained in `mathematical_model.ipynb`, with some dependencies in `helpers.py`. Although scipy is used, a sage kernel should be employed as some syntax is sage specific. 


`dataWrangling.R` is used to process data from the access study. Raw data is not uploaded to github for patient privacy but is available from the [Immune Tolerance Network’s Trial Share platform](https://www.itntrialshare.org/project/Studies/ITN034AIPUBLIC/Study Data/begin.view). 

`run_fixed.py` and `run_stochastic.py` can be used with data from `dataWrangling.R` to fit the original model (fixed) or the model with stochastic fluctuations (stochastic). These are designed to be used with SLURM via `jobscript.sh`.

`view_saved.ipynb` and `fitParamViz.R` is then used for interpreting and visualizing model fits that were run in SLURM jobs using `run_fixed.py` and `run_stochastic.py`.
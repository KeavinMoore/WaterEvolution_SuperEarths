# MO-cycling-loss Model

Python model of water exchange on Earth-like planets of mass between 1-10$M_{\oplus}$ orbiting an M-dwarf star.

## Model Requirements

- Python 3.6 
- Numpy
- Scipy
- Matplotlib.pyplot
- A computer :)

## How to run the model

1. Assign the desired model paramaters in the file [MOparams.py](https://github.com/bendavid791/MO-cycling-loss/blob/main/MOparams.py).

2. Open the file [MO_main.ipynb](https://github.com/bendavid791/MO-cycling-loss/blob/main/MO_main.ipynb) and add your path to the current working directory.

3. Depending on whether you have the Magma Ocean data files, set `have_MO_data` to either `False` or `True`. `False` will use existing Magma ocean data files from the current working directory with `M`, `D`, `NUM_EARTH_OCEANS`, and `HABITABLE_ZONE` that you assigned in step 1 in [MOparams.py](https://github.com/bendavid791/MO-cycling-loss/blob/main/MOparams.py). `True` will not use any existing data files, rewriting them as the model runs.

4. Run all cells in the file [MO_main.ipynb](https://github.com/bendavid791/MO-cycling-loss/blob/main/MO_main.ipynb). This will produce 3 output files; cyclingresults and MOresults which are the simulation results, and a plot of the water inventories accross the simulations (using the simulation results data).

# Young Open Cluster Simulation
## How to run
The source code consists of two files: simulate.py and postprocessing.py.
The first runs the simulation and saves some datafiles, the second uses these files to visualise the results.
The files should both be run from the main project folder e.g. `src/simulate.py` (or `python3 src/simulate.py` if the script doesn't detect your python virtual environment automatically).
To run the scripts please use an environment with the following packages installed:
- Amuse
- Matplotlib
- Numpy
- Pandas
- Tqdm

### simulate.py
The main file for simulating an open cluster with stars and gas.
When run without any command line arguments, the file will run a simple simulation containing 10 stars and 100 gas particles for demonstration purposes.
By running `src/simulate.py --help` a promt will show up displaying the different command line arguments that can be passed down to alter the initial conditions of the simulation.

### postprocessing.py
The data visualisation script.
Uses the data saved by `simulate.py` to create several images and/or animations.
The command line argument `--run <name>` can be used to specify which data to use in the visualisation.

## The Report
The LaTeX files for the report can be found in the report directory.
A pdf can be generated using the provided Makefile, by running `make all` from the report directory.
To generate this pdf, `pdflatex` and `biber` should be installed.

---

## Minimum requirements
Simulate a open cluster with the following points implemented:
- Simulate the stars:
  - [x] Generate 1000 stars from a Salpeter IMF.
  - [x] Tune the most massive star to be 30 Msun.
  - [x] Distribute the stars in a fractal distribution.
- Simulate the gas:
  - [x] Distribute the gas in a Plummer sphere.
  - [x] Remove the gas around the stars (create a Swiss cheese), based on the mass of the stars.
- Stellar winds and Supernovae:
  - [x] Add stellar evolution (SEBA)
  - [x] Implement stellar winds module
  - [x] Use 10000 sph particles
  - [ ] Time step 1kyr - 100kyr when the 30MSun star goes supernova
- Cool stuff:
  - [ ] Make movie of simulation

## Research question
> What is the contribution of the supernova to clearing the gas from a newly formed cluster?

## Intermediate errors to fix
- [ ] de Fi error met code -1
- [ ] De initializatie positie van gasdeeltjes door stellar winds
- [ ] gatenkaas beter maken
- [ ] timesteps verbeteren, timesteps van bridge kan max 0.5 keer de totale zijn en stellar evolution heeft eigen timestep, dus die roep je niet aan

## Points of improvement
- [ ] Use MESA for the stellar evolution code
- [ ] Tweak the star formation efficiency 
- [ ] Make the plummer sphere radius much bigger, how does this impact the results?

---
## Contributors:
- Rick Dullaart (s1993879)
- Wouter van Tol (s2041340)
- Rutger Rijnenberg (s1829777)


# Young Open Cluster Simulation
## Minimum requirements
Simulate a open cluster with the following points implemented:
- Simulate the stars:
  - [x] Generate 1000 stars from a Salpeter IMF.
  - [/] Tune the most massive star to be 30 Msun.
  - [x] Distribute the stars in a fractal distribution.
- Simulate the gas:
  - [x] Distribute the gas in a Plummer sphere.
  - [/] Remove the gas around the stars (create a Swiss cheese), based on the mass of the stars.
- Stellar winds and Supernovae:
  - [x] Add stellar evolution (SEBA)
  - [x?] Implement stellar winds module
  - [x] Use 10000 sph particles
  - [ ] Time step 1kyr - 100kyr

## Research question
> What is the contribution of the supernova to clearing the gas from a newly formed cluster?

## Points of improvement
- [ ] Use MESA for the stellar evolution code
- [ ] Tweak the star formation efficiency 
- [ ] Make the plummer sphere radius much bigger, how does this impact the results?

---
## Contributors:
- Rick Dullaart (s1993879)
- Wouter van Tol (s2041340)
- Rutger Rijnenberg (s1829777)

<!--- take a mass function (salpeter) and 1000 stars, maybe 1 30 solar mass starr or tune such that we have at least 1 large star. Take a fractal distribution of stars since it is less bound and stable. Take gas plummer sphere around this distribution of stars. 
As initial conditions eat away the gas around the stars proportional to the mass of the stars, like a swiss cheese. 
we need star evolution, use SEBA, for extra points we can use MESA, use 10000 shp particles.
One of our problems is the bridge timesteps after a supernova starts, then we can use 1000 or 10000 years timesteps.
Make the most massive star 30 solar masses and change it to virial equilibrium.
Is the gas blown away by the stellar winds or by the supernova?
What is the contribution of the supernova to clearing the gas?
present 14th december, deadline is 23rd of december
see if the swiss cheese approach is realistic --->
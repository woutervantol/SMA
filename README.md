# Young Open Cluster Simulation
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
  - [x?] Implement stellar winds module
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

<!--- take a mass function (salpeter) and 1000 stars, maybe 1 30 solar mass starr or tune such that we have at least 1 large star. Take a fractal distribution of stars since it is less bound and stable. Take gas plummer sphere around this distribution of stars. 
As initial conditions eat away the gas around the stars proportional to the mass of the stars, like a swiss cheese. 
we need star evolution, use SEBA, for extra points we can use MESA, use 10000 shp particles.
One of our problems is the bridge timesteps after a supernova starts, then we can use 1000 or 10000 years timesteps.
Make the most massive star 30 solar masses and change it to virial equilibrium.
Is the gas blown away by the stellar winds or by the supernova?
What is the contribution of the supernova to clearing the gas?
present 14th december, deadline is 23rd of december
see if the swiss cheese approach is realistic --->

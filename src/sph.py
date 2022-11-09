#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


print(f"Hi, I'm proccessor {rank} out of {size}")


from amuse.community.fi.interface import Fi
from amuse.units import (units, constants)
from amuse.couple import bridge
from amuse.ext.composition_methods import *
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.ext.protodisk import ProtoPlanetaryDisk
from tqdm import tqdm
import numpy as np

### TEST
from amuse.ext.salpeter import new_salpeter_mass_distribution

n_stars = 1000
alpha_IMF = -2.35

m_stars = new_salpeter_mass_distribution(n_stars, 0.1|units.MSun, 100|units.MSun, alpha_IMF)


#plot the initial mass distribution in loglog
# hist, bins = np.histogram(m_stars.value_in(units.MSun), bins=20)
# logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
# plt.hist(m_stars.value_in(units.MSun), bins=logbins, density=True)
# plt.xscale('log')
# plt.show()


r_cluster = 1.0 | units.parsec

#converter is nodig omdat het anders dimensieloos is, nu kunnen we initial conditions in SI ingeven
from amuse.units import nbody_system
converter=nbody_system.nbody_to_si(m_stars.sum(),r_cluster)

from amuse.community.ph4.interface import ph4
gravity = ph4(converter)
### ENDTEST


hydro = Fi(converter, mode='g6lib')
hydro.parameters.use_hydro_flag = True
hydro.parameters.radiation_flag = False
hydro.parameters.gamma = 1
hydro.parameters.isothermal_flag = True
hydro.parameters.integrate_entropy_flag = False
hydro.parameters.timestep = 0.01 | units.Myr 
hydro.parameters.verbosity = 0
hydro.parameters.eps_is_h_flag = False    # h_smooth is constant
eps = 0.1 | units.au
hydro.parameters.gas_epsilon = eps
hydro.parameters.sph_h_const = eps
print("print 3")

gravhydro = bridge.Bridge(use_threading=False) #, method=SPLIT_4TH_S_M4)
gravhydro.add_system(gravity, (hydro,))
gravhydro.add_system(hydro, (gravity,))
gravhydro.timestep = 0.2 | units.Myr
print("print 4")


def gravity_hydro_bridge(gravity, hydro, gravhydro, bodies,
                         t_end):

    gravity_initial_total_energy = gravity.get_total_energy() + hydro.get_total_energy()
    model_time = 0 | units.Myr
    dt = 0.012|units.yr  #1.0*Pinner
    t_steps = np.arange(model_time.value_in(units.Myr), t_end.value_in(units.Myr), dt.value_in(units.Myr))
    print("print 5")
    for model_time in tqdm(t_steps):
        orbit_planet = orbital_elements_from_binary(bodies[:2], G=constants.G)
        orbit_moon = orbital_elements_from_binary(bodies[1:3], G=constants.G)
        print("Planet:", "ae=", orbit_planet[2].in_(units.AU), orbit_planet[3])
        print("Moon:", "ae=", orbit_moon[2].in_(units.AU), orbit_moon[3])
        
        dE_gravity = gravity_initial_total_energy/(gravity.get_total_energy()+hydro.get_total_energy())
        print("Time:", model_time.in_(units.day), \
              "dE=", dE_gravity)#, dE_hydro

        gravhydro.evolve_model(model_time)
        channel["to_stars"].copy()
        channel["to_disk"].copy()
        channel["to_moon"].copy()
        
        print("S=", bodies[:3])
        print("g=", gravity.particles)
        print(gravity.particles.y.in_(units.au), moon.y.in_(units.au))

    gravity.stop()
    hydro.stop()

t_end = 1.0 | units.yr
gravity_hydro_bridge(gravity, hydro, gravhydro, 
                     bodies, t_end)
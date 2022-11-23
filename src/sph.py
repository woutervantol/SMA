#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


print(f"Hi, I'm proccessor {rank} out of {size}")

import matplotlib.pyplot as plt
import numpy as np
from amuse.community.fi.interface import Fi
from amuse.community.ph4.interface import ph4
from amuse.couple import bridge
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.ext.salpeter import new_salpeter_mass_distribution
from amuse.ext import stellar_wind
from amuse.ic.plummer import new_plummer_sphere
from amuse.lab import (new_plummer_gas_model, new_plummer_sphere)
from amuse.units import (units, constants, nbody_system)
from tqdm import tqdm


#create stars with masses, positions and velocities and put them in the ph4 module
n_stars = 1000
alpha_IMF = -2.35

m_stars = new_salpeter_mass_distribution(n_stars, 0.1|units.MSun, 100|units.MSun, alpha_IMF)
total_mass = np.sum(m_stars)
print(total_mass)
r_cluster = 1.0 | units.parsec
#converter is nodig omdat het anders dimensieloos is, nu kunnen we initial conditions in SI ingeven
converter=nbody_system.nbody_to_si(m_stars.sum(),r_cluster)

bodies=new_plummer_sphere(n_stars, convert_nbody=converter)
bodies.scale_to_standard(converter)
bodies.mass = m_stars
gravity = ph4(converter)
gravity.particles.add_particles(bodies)



## Maak gasdeeltjes
#disc = ProtoPlanetaryDisk(n_stars,
#                              convert_nbody=converter,
#                              Rmin=0.01 | units.parsec,
#                              Rmax=r_cluster,
#                              q_out=10.0,
#                              discfraction=0.01).result

Ngas = 1000
gas = new_plummer_gas_model(Ngas, convert_nbody=converter)
### ENDTEST



#create a hydro code and a gas distribution and put the gas in the hydro code
hydro = Fi(converter, mode='g6lib')
hydro.parameters.use_hydro_flag = True
hydro.parameters.radiation_flag = False
hydro.parameters.gamma = 1
hydro.parameters.isothermal_flag = True
hydro.parameters.integrate_entropy_flag = False
# hydro.parameters.timestep = 0.01 | units.Myr
hydro.parameters.verbosity = 0
hydro.parameters.eps_is_h_flag = False    # h_smooth is constant
eps = 0.1 | units.au
hydro.parameters.gas_epsilon = eps
hydro.parameters.sph_h_const = eps


Ngas = 10000
gas = new_plummer_gas_model(Ngas, convert_nbody=converter) #Note: this is virialized gas, so it has velocities
hydro.particles.add_particles(gas)

supernova_gas_velocity = 12.9 | units.kms
SNstar = bodies[np.argmax(bodies.mass)]
NSNgas = 1000
SNconverter = nbody_system.nbody_to_si(SNstar.mass, 1|units.RSun)
SNgas = new_plummer_gas_model(NSNgas, SNconverter)
#TODO: voeg toe velocities (of u) in de radial direction
print(dir(SNgas))
print(SNgas.u)

# print(bodies.mass)
# SNstar = bodies[np.argmax(bodies.mass)]
# print(SNstar.mass)
# NSNgas = 1000
# SNgas = stellar_wind.new_stellar_wind(SNstar.mass, target_gas=gas, timestep=0)

# print(SNgas)


#bridge the codes
gravhydro = bridge.Bridge(use_threading=False) #, method=SPLIT_4TH_S_M4)
gravhydro.add_system(gravity, (hydro,))
gravhydro.add_system(hydro, (gravity,))
# gravhydro.timestep = 0.2 | units.Myr

channel = {"from_stars": bodies.new_channel_to(gravity.particles),
            "to_stars": gravity.particles.new_channel_to(bodies),
            "from_gas": gas.new_channel_to(hydro.particles),
            "to_gas": hydro.particles.new_channel_to(gas)}
        

def gravity_hydro_bridge(gravity, hydro, gravhydro, bodies,
                         t_end):

    gravity_initial_total_energy = gravity.get_total_energy() + hydro.get_total_energy()
    model_time = 0 | units.Myr
    dt = 0.1|units.Myr  #1.0*Pinner

    t_steps = np.arange(model_time.value_in(units.Myr), t_end.value_in(units.Myr), dt.value_in(units.Myr)) | units.Myr
    fig, ax = plt.subplots(3, 3)
    ax = ax.flatten()
    for i, t in enumerate(tqdm(t_steps)):
        dE_gravity = gravity_initial_total_energy/(gravity.get_total_energy()+hydro.get_total_energy())
        print(dE_gravity, t)
        gravhydro.evolve_model(t)
        channel["to_stars"].copy()
        channel["to_gas"].copy()

        # print("gravitational energy: ", bodies.potential_energy())
        # print("kinetic energy: ", bodies.kinetic_energy())
        print( - bodies.potential_energy() / bodies.kinetic_energy())
        if i < 91 and i%10:
            ax[i//10].scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), s=1)
            ax[i//10].scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), s=1, c=np.log(m_stars.value_in(units.MSun)))
            ax[i//10].scatter(bodies[0].x.value_in(units.parsec), bodies[0].y.value_in(units.parsec), s=5, c="red")
    plt.show()

    gravity.stop()
    hydro.stop()

t_end = 10.0 | units.Myr
gravity_hydro_bridge(gravity, hydro, gravhydro, 
                     bodies, t_end)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


print(f"Hi, I'm proccessor {rank} out of {size}")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
from amuse.community.fractalcluster.interface import new_fractal_cluster_model
from amuse.community.seba.interface import SeBa

# Toevoegen:
# - Channel van stellar winds
# - Channel van stellar evolution
# - 

#seeds for which the highest mass star has mass m with:  29.5Msun < m < 30.5MSun
seeds = [112, 134, 216, 275, 309, 317, 458, 596, 661, 775, 836, 848, 873, 930, 939]
np.random.seed(seeds[np.random.randint(0, len(seeds))]) #take random seed from valid seeds

def create_cheese(gas, stars, r):
    cheesegas = gas.select(lambda gaspos: ((stars.position-gaspos).lengths()<r).any(),["position"])
    return gas.difference(cheesegas)

#create stars with masses, positions and velocities and put them in the ph4 module
n_stars = 1000
alpha_IMF = -2.35
m_stars = new_salpeter_mass_distribution(n_stars, 0.1|units.MSun, 100|units.MSun, alpha_IMF)
total_mass = np.sum(m_stars)
print("max mass star:", m_stars[np.argmax(m_stars)])
r_cluster = 1.0 | units.parsec
converter=nbody_system.nbody_to_si(m_stars.sum(),r_cluster)

bodies=new_fractal_cluster_model(n_stars,fractal_dimension= 1.6, convert_nbody=converter)
bodies.scale_to_standard(converter)
bodies.mass = m_stars
SNstar = bodies[np.argmax(bodies.mass)]
gravity = ph4(converter)
gravity.particles.add_particles(bodies)


## Maak gasdeeltjes
Ngas = 10000
gas = new_plummer_gas_model(Ngas, convert_nbody=converter)

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), gas.z.value_in(units.parsec), s=1)
gas = create_cheese(gas, bodies, 0.6 | units.parsec)
# ax.scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), gas.z.value_in(units.parsec), s=1)
# ax.scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), bodies.z.value_in(units.parsec), s=4, color="black")
# plt.show()





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
Mgas = np.sum(gas.mass)
mgas = Mgas/Ngas


#bridge the codes
gravhydro = bridge.Bridge(use_threading=False) #, method=SPLIT_4TH_S_M4)
gravhydro.add_system(gravity, (hydro,))
gravhydro.add_system(hydro, (gravity,))
# gravhydro.timestep = 0.2 | units.Myr

channel = {"from_stars": bodies.new_channel_to(gravity.particles),
            "to_stars": gravity.particles.new_channel_to(bodies),
            "from_gas": gas.new_channel_to(hydro.particles),
            "to_gas": hydro.particles.new_channel_to(gas),
            #"from_SNgas": SNgas.new_channel_to(hydro.particles),
            #"to_SNgas": hydro.particles.new_channel_to(SNgas)
            }
            
# Stellar evolution
evolution = SeBa()
evolution.particles.add_particles(bodies)
ch_e2g = evolution.particles.new_channel_to(gravity.particles)
ch_g2e = gravity.particles.new_channel_to(bodies)
ch_e2b = evolution.particles.new_channel_to(bodies)
ch_e2b.copy()
        

# Stellar wind  p. 223
from amuse.ext.stellar_wind import new_stellar_wind
dt = 0.1|units.Myr  #1.0*Pinner     # Gekopieerd uit gravity_hydro_bridge()
wind = new_stellar_wind(mgas, target_gas=gas, timestep=dt, derive_from_evolution=True)
wind.particles.add_particles(bodies)
channel_to_wind = bodies.new_channel_to(wind.particles)

# # Stellar wind (p.214-216 in book)
# n_wind = 10000
# from amuse.lab import Particles
# wind_gas = Particles(n_wind)
# wind_gas.mass = 
# # give_particles_some_properties(new_gas)
# gas.add_particles(wind_gas)
# gas.synchronize_to(hydro.gas)


def gravity_hydro_bridge(gravity, hydro, gravhydro, bodies, t_end):

    gravity_initial_total_energy = gravity.get_total_energy() + hydro.get_total_energy()
    model_time = 0 | units.Myr
    dt = 0.1|units.Myr  #1.0*Pinner

    t_steps = np.arange(model_time.value_in(units.Myr), t_end.value_in(units.Myr), dt.value_in(units.Myr)) | units.Myr
    print("t_steps:", t_steps)
    fig, ax = plt.subplots(3, 3)
 #   fig2, ax2 = plt.subplots(3, 3)
    ax = ax.flatten()
 #   ax2 = ax2.flatten()
    for i, t in enumerate(tqdm(t_steps)):
        dE_gravity = gravity_initial_total_energy/(gravity.get_total_energy()+hydro.get_total_energy())
        print("dE:", dE_gravity, "; t=", t)
        evolution.evolve_model(t)
        channel_to_wind.copy()      # wind with hydro and grav: Book 8.1.1 p.323
        wind.evolve_model(t)
        ch_e2g.copy()
        ch_e2b.copy()        
        gravhydro.evolve_model(t)
        channel["to_stars"].copy()
        channel["to_gas"].copy()
        #channel["to_SNgas"].copy()

        current_gasmass = np.sum(gas.mass)
        print("Total mass of gas:", current_gasmass)
        current_gasnumber = current_gasmass/mgas
        print("# of gass particles:", current_gasnumber)
        
        bodies_pd = pd.DataFrame(np.array(bodies.stellar_type.number), columns=["stellar_type"])
        print("\n", bodies_pd.value_counts(), "\n")

        # print("gravitational energy: ", bodies.potential_energy())
        # print("kinetic energy: ", bodies.kinetic_energy())
        print("-Ep/Ek:", - bodies.potential_energy() / bodies.kinetic_energy())
        print("Total mass:", np.sum(bodies.mass) | units.MSun)
        if i < 9:#1:# and i%10==0:
            ax[i].scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), s=1)
            ax[i].scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), s=1)#, c=np.log(m_stars.value_in(units.MSun)))

            # ax2[i].scatter(SNgas.x.value_in(units.parsec), SNgas.y.value_in(units.parsec), s=1)
            ax[i].scatter(bodies[np.argmax(bodies.mass)].x.value_in(units.parsec), bodies[np.argmax(bodies.mass)].y.value_in(units.parsec), s=5, c="red")
        #if i == 8:     #Weggecomment om een supernova te vinden.
            #break
    plt.show()
    
    evolution.stop()
    gravity.stop()
    hydro.stop()

t_end = 50.0 | units.Myr
gravity_hydro_bridge(gravity, hydro, gravhydro, 
                     bodies, t_end)





# Interessant boek? https://misaladino.com/wp-content/uploads/2019/11/Thesis_Martha_Irene.pdf

#### Stellar types ####
""" 
"deeply or fully convective low mass MS star",  # 0
        "Main Sequence star",  # 1
        "Hertzsprung Gap",  # 2
        "First Giant Branch",  # 3
        "Core Helium Burning",  # 4
        "First Asymptotic Giant Branch",  # 5
        "Second Asymptotic Giant Branch",  # 6
        "Main Sequence Naked Helium star",  # 7
        "Hertzsprung Gap Naked Helium star",  # 8
        "Giant Branch Naked Helium star",  # 9
        "Helium White Dwarf",  # 10
        "Carbon/Oxygen White Dwarf",  # 11
        "Oxygen/Neon White Dwarf",  # 12
        "Neutron Star",  # 13
        "Black Hole",  # 14
        "Massless Supernova",  # 15
        "Unknown stellar type",  # 16
        "Pre-main-sequence Star",  # 17
        "Planet",  # 18
        "Brown Dwarf",  # 19
"""        
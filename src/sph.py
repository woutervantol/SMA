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

dt = 0.2 | units.Myr
dt_bridge = 0.1 | units.Myr  #1.0*Pinner

#seeds for which the highest mass star has mass m with:  29.5Msun < m < 30.5MSun
seeds = [112, 134, 216, 275, 309, 317, 458, 596, 661, 775, 836, 848, 873, 930, 939]
np.random.seed(seeds[np.random.randint(0, len(seeds))]) #take random seed from valid seeds

def create_cheese(gas, stars, r):
    cheesegas = gas.select(lambda gaspos: ((stars.position-gaspos).lengths()<r).any(),["position"])
    return gas.difference(cheesegas)    # Testcode hiervan ook opslaan in github! Kunnen ze ook naar kijken. Aanpassen op basis van massa ster (eerst lage massas kazen)

#create stars with masses, positions and velocities and put them in the ph4 module
n_stars = 10    # 1000
alpha_IMF = -2.35
m_stars = new_salpeter_mass_distribution(n_stars, 0.1|units.MSun, 100|units.MSun, alpha_IMF)
total_mass = np.sum(m_stars)

#
m_stars[np.argmax(m_stars)] = 30 | units.MSun   # Added so we have a massive star for other n_stars
#

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
Ngas = 1000  # 10000
gas = new_plummer_gas_model(Ngas, convert_nbody=converter)

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), gas.z.value_in(units.parsec), s=1)
gas = create_cheese(gas, bodies, 0.6 | units.parsec) # Gasdeeltjes weghalen op basis van massa van de ster
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
hydro.parameters.timestep = dt
hydro.parameters.verbosity = 0
hydro.parameters.eps_is_h_flag = False    # h_smooth is constant
eps = 0.1 | units.au
hydro.parameters.gas_epsilon = eps
hydro.parameters.sph_h_const = eps

# Ngas = 10000
# gas = new_plummer_gas_model(Ngas, convert_nbody=converter) #Note: this is virialized gas, so it has velocities
hydro.particles.add_particles(gas)
Mgas = np.sum(gas.mass)
mgas = Mgas/Ngas


#bridge the codes
gravhydro = bridge.Bridge(use_threading=False) #, method=SPLIT_4TH_S_M4)
gravhydro.add_system(gravity, (hydro,))
gravhydro.add_system(hydro, (gravity,))
gravhydro.timestep = dt_bridge # min 2x de output timestep
            
# Stellar evolution
evolution = SeBa()
evolution.particles.add_particles(bodies)
ch_e2g = evolution.particles.new_channel_to(gravity.particles)
ch_g2e = gravity.particles.new_channel_to(bodies) # Unnecessary?
ch_e2b = evolution.particles.new_channel_to(bodies)
ch_e2b.copy()
        

# Stellar wind  p. 223
from amuse.ext.stellar_wind import new_stellar_wind
wind = new_stellar_wind(mgas, target_gas=gas, timestep=dt, derive_from_evolution=True)
wind.particles.add_particles(bodies)
channel_to_wind = bodies.new_channel_to(wind.particles)

channel = {"from_stars": bodies.new_channel_to(gravity.particles),
            "to_stars": gravity.particles.new_channel_to(bodies),
            "from_gas": gas.new_channel_to(hydro.particles),
            "to_gas": hydro.particles.new_channel_to(gas),
            "evo_to_grav": evolution.particles.new_channel_to(gravity.particles),
            "evo_to_stars": evolution.particles.new_channel_to(bodies),
            "stars_to_wind": channel_to_wind
            }

# # Stellar wind (p.214-216 in book)
# n_wind = 10000
# from amuse.lab import Particles
# wind_gas = Particles(n_wind)
# wind_gas.mass = 
# # give_particles_some_properties(new_gas)
# gas.add_particles(wind_gas)
# gas.synchronize_to(hydro.gas)

def ninestepplot(bodies, gas, i, t, maintitle, savename, fig, ax, fig_complete):
    if fig_complete == True:
        fig, ax = plt.subplots(3, 3)
        fig.suptitle(maintitle)
        ax = ax.flatten()
        fig_complete = False
    if i < 9:
        ax[i].scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), s=1)
        ax[i].scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), s=1)#, c=np.log(m_stars.value_in(units.MSun)))
        ax[i].scatter(bodies[np.argmax(bodies.mass)].x.value_in(units.parsec), bodies[np.argmax(bodies.mass)].y.value_in(units.parsec), s=5, c="red")
        ax[i].set_title('t = '+str(t))
    if i == 8:
        fig.savefig(savename)
        fig_complete = True
        plt.show()
    return fig, ax, fig_complete

def simulate(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t):
    evolution.evolve_model(t)    
    channel["evo_to_grav"].copy()
    channel["evo_to_stars"].copy()

    channel["stars_to_wind"].copy()      # wind with hydro and grav: Book 8.1.1 p.323
    wind.evolve_model(t)
    gas.synchronize_to(hydro.particles)  # Deze code zou moeten kloppen, maar zorgt wel voor een uiteindelijke crash: navragen
                                            # Crash: "amuse.support.exceptions.AmuseException: Error when calling 'get_position' of a 'Fi', errorcode is -1"
    # print(dir(gas))
    # print(dir(gas[-1]))
    # print("velocity of gas:", np.sqrt((gas.vx)**2+(gas.vy)**2+(gas.vz)**2))
    # if gas[-1] != gas[99]:
    #    print("Positions:", gas.x[100])

    gravhydro.evolve_model(t)
    channel["to_stars"].copy()
    channel["to_gas"].copy()
    return gravity, hydro, gravhydro, evolution, wind, bodies, gas

def star_control(bodies, n_stars):
        bodies_pd = pd.DataFrame(np.array(bodies.stellar_type.number), columns=["stellar_type"])
        print("\n", bodies_pd.value_counts(), "\n")
        if int(bodies_pd.value_counts().loc[1]) < n_stars:
            if list(bodies_pd.value_counts().index[1]) == [14]: 
                print("Black holes:", int(bodies_pd.value_counts().loc[14]))
            if list(bodies_pd.value_counts().index[1]) == [4]: 
                print("Core helium burning stars:", int(bodies_pd.value_counts().loc[4]))
        return list(bodies_pd.value_counts().index[-1])[0]

def gravity_hydro_bridge(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t_end, dt, dt_bridge, n_stars):
    gravity_initial_total_energy = gravity.get_total_energy() + hydro.get_total_energy()
    model_time = 0 | units.Myr
    dt_SN = 0.001 | units.Myr
    
    t_steps = np.arange(model_time.value_in(units.Myr), t_end.value_in(units.Myr), dt.value_in(units.Myr)) | units.Myr
    t_steps_coderen = np.concatenate((t_steps[:1], t_steps[28:]), axis=None)    # Voor snelheid coderen

    fig_complete = True
    fig, ax = False, False
    for i, t in enumerate(tqdm(t_steps)):
        gravity, hydro, gravhydro, evolution, wind, bodies, gas = simulate(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t)
        if i < 9:
            fig, ax, fig_complete = ninestepplot(bodies, gas, i, t, "Cluster at initialization", "Replace_initialization.png", fig, ax, fig_complete)

        dE_gravity = gravity_initial_total_energy/(gravity.get_total_energy()+hydro.get_total_energy())
        print("dE:", dE_gravity, "; t=", t)
        current_gasmass = np.sum(gas.mass)
        if i == 0:
            start_mass = current_gasmass
        print("Mass change in stars:", current_gasmass-start_mass)

        print("Total mass of gas:", current_gasmass)
        current_gasnumber = current_gasmass/mgas        
        print("# of gass particles:", current_gasnumber)
        print("-Ep/Ek:", - bodies.potential_energy() / bodies.kinetic_energy())
        print("Total mass:", np.sum(bodies.mass) | units.MSun)

        most_advanced_type = star_control(bodies, n_stars)
        print("MAT:", most_advanced_type)
        if most_advanced_type == 14:
            t_SN = t
            break
    
    t_steps_supernova = np.arange(t_SN.value_in(units.Myr), t_end.value_in(units.Myr), dt_SN.value_in(units.Myr)) | units.Myr
    for i, t in enumerate(tqdm(t_steps_supernova)):
        gravity, hydro, gravhydro, evolution, wind, bodies, gas = simulate(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t)
        if i < 9:
            fig, ax, fig_complete = ninestepplot(bodies, gas, i, t, "Cluster after supernova", "Replace_after_supernova.png", fig, ax, fig_complete)
        if (i>2000) & (i<2010):
            fig, ax, fig_complete = ninestepplot(bodies, gas, i-2001, t, "Cluster longer after supernova", "Replace_longer_after_supernova.png", fig, ax, fig_complete)

        dE_gravity = gravity_initial_total_energy/(gravity.get_total_energy()+hydro.get_total_energy())
        print("dE:", dE_gravity, "; t=", t)
        current_gasmass = np.sum(gas.mass)
        print("Total mass of gas:", current_gasmass)
        current_gasnumber = current_gasmass/mgas        
        print("# of gass particles:", current_gasnumber)
        print("-Ep/Ek:", - bodies.potential_energy() / bodies.kinetic_energy())
        print("Total mass:", np.sum(bodies.mass) | units.MSun)

        most_advanced_type = star_control(bodies, n_stars)
        print("MAT:", most_advanced_type)


    plt.show()
    
    evolution.stop()
    gravity.stop()
    hydro.stop()

t_end = 100.0 | units.Myr
gravity_hydro_bridge(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t_end, dt, dt_bridge, n_stars)





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


# 2 timesteps:
# De bridge timestep (minimaal half van de output timestep)
# De output timestep. Hoe vaak wil je wind creeren? Stellar evolution heeft zn eigen timestep en die kan je in principe aanvragen.

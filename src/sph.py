#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from amuse.community.fi.interface import Fi
from amuse.community.ph4.interface import ph4
from amuse.couple import bridge
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.ext.salpeter import new_salpeter_mass_distribution
from amuse.ext.stellar_wind import new_stellar_wind
from amuse.ext import stellar_wind
from amuse.ic.plummer import new_plummer_sphere
from amuse.lab import (new_plummer_gas_model, new_plummer_sphere)
from amuse.units import (units, constants, nbody_system)
from tqdm import tqdm
from amuse.community.fractalcluster.interface import new_fractal_cluster_model
from amuse.community.seba.interface import SeBa

np.random.seed(1)

# x = 1|units.parsec
# t = 0.1|units.Myr
# print((x/t).value_in(units.ms))

dt = 0.1 | units.Myr
dt_winds = 0.05 | units.Myr
dt_hydro = 0.01 | units.Myr
dt_bridge = 0.05 | units.Myr  #1.0*Pinner

def create_swiss_cheese_gas(initial_gas, stars):
    """Create Swiss cheese holes in a gas according to a star distribution.
    
    Removes gas around the position of stars where the mass of the removed gas
    equals the mass of the star. This results in a gas distribution with holes
    or a Swiss cheese like distribution.

    Parameters
    ----------
    initial_gas
        The initial gas distribution in where the holes will be cut into.
    stars
        The distribution of stars used to remove gas from the initial gas
        distribution.

    Returns
    -------
    newgas
        The new swiss cheese gas distribution.
    """
    newgas = initial_gas.copy()
    tryradii = np.logspace(np.log10(0.001), np.log10(1), 50) | units.parsec
    sorted_stars = np.sort(stars.mass)
    for m_star in sorted_stars:
        star = stars[stars.mass == m_star]
        for r in tryradii:
            gas_to_remove = newgas[(star.position - newgas.position).lengths()<r]
            if np.sum(gas_to_remove.mass) > star.mass:
                newgas.remove_particle(gas_to_remove)
                break
    return newgas


#create stars with masses, positions and velocities and put them in the ph4 module
n_stars = 10
alpha_IMF = -2.35
while True:
    m_stars = new_salpeter_mass_distribution(n_stars, 0.1|units.MSun, 100|units.MSun, alpha_IMF)
    masslist = np.sort(m_stars)
    if masslist[-1].value_in(units.MSun) > 27 and masslist[-1].value_in(units.MSun) < 33 and masslist[-2].value_in(units.MSun) < 27:
        break
total_mass = np.sum(m_stars)
print("max mass star:", m_stars[np.argmax(m_stars)])


r_cluster = 1.0 | units.parsec
gravconverter=nbody_system.nbody_to_si(total_mass,r_cluster)

bodies=new_fractal_cluster_model(n_stars,fractal_dimension= 1.6, convert_nbody=gravconverter)
bodies.scale_to_standard(gravconverter)
bodies.mass = m_stars
SNstar = bodies[np.argmax(bodies.mass)]
gravity = ph4(gravconverter)
gravity.particles.add_particles(bodies)



## Maak gasdeeltjes
Ngas = 1000  # 10000
total_gas_mass = 5*total_mass
gasconverter=nbody_system.nbody_to_si(total_gas_mass+total_mass,r_cluster)
gas = new_plummer_gas_model(Ngas, convert_nbody=gasconverter)
# gas.scale_to_standard(gasconverter)

gas = create_swiss_cheese_gas(gas, bodies)

#create a hydro code and a gas distribution and put the gas in the hydro code
hydro = Fi(gravconverter, mode='g6lib')
hydro.parameters.gamma = 1
hydro.parameters.isothermal_flag = True
hydro.parameters.integrate_entropy_flag = False
hydro.parameters.timestep = dt_hydro
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
wind = new_stellar_wind(mgas, target_gas=gas, timestep=dt_winds, derive_from_evolution=True)
wind.particles.add_particles(bodies)


channel = {"from_stars": bodies.new_channel_to(gravity.particles),
            "to_stars": gravity.particles.new_channel_to(bodies),
            "from_gas": gas.new_channel_to(hydro.particles),
            "to_gas": hydro.particles.new_channel_to(gas),
            "evo_to_grav": evolution.particles.new_channel_to(gravity.particles),
            "evo_to_stars": evolution.particles.new_channel_to(bodies),
            "stars_to_wind": bodies.new_channel_to(wind.particles)
            }

def ninestepplot(bodies, gas, i, t, maintitle, savename, fig, ax, fig_complete):
    if fig_complete == True:
        fig, ax = plt.subplots(3, 3)
        fig.suptitle(maintitle)
        ax = ax.flatten()
        fig_complete = False
    if i < 9:
        ax[i].scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), s=1)
        ax[i].scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), s=1)
        ax[i].scatter(bodies[np.argmax(bodies.mass)].x.value_in(units.parsec), bodies[np.argmax(bodies.mass)].y.value_in(units.parsec), s=5, c="red")
        ax[i].set_title('t = '+str(t))
    if i == 8:
        fig.savefig(savename)
        fig_complete = True
        plt.show()
    return fig, ax, fig_complete

def onestepplot():
    plt.scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), s=1, label="gas")
    plt.scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), s=1, label="stars")
    plt.scatter(bodies[np.argmax(bodies.mass)].x.value_in(units.parsec), bodies[np.argmax(bodies.mass)].y.value_in(units.parsec), s=5, c="red")
    plt.legend()
    plt.xlim(-5, 5)
    plt.ylim(-5, 5)
    plt.xlabel("parsec")
    plt.ylabel("parsec")
    plt.show()

def plot_energies(Us, Ks, Ts):
    plt.semilogy(Ts, Us, label="potential energy")
    plt.semilogy(Ts, Ks, label="kinetic energy")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Energy (m^2 kg s^-2)")
    plt.legend()
    plt.show()

def star_control(bodies, n_stars):
        bodies_pd = pd.DataFrame(np.array(bodies.stellar_type.number), columns=["stellar_type"])
        return list(bodies_pd.value_counts().index[-1])[0]

def delete_outofbounds():
    mask = gas.position.lengths() >= 1e18|units.m
    selection = gas[mask]
    gas.remove_particles(selection)




def simulate(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t):
    evolution.evolve_model(t)
    channel["evo_to_grav"].copy()
    channel["evo_to_stars"].copy()
    channel["stars_to_wind"].copy()      # wind with hydro and grav: Book 8.1.1 p.323
    wind.evolve_model(t)
    delete_outofbounds()
    gas.synchronize_to(hydro.particles)
    gravhydro.evolve_model(t)
    channel["to_stars"].copy()
    channel["to_gas"].copy()
    return gravity, hydro, gravhydro, evolution, wind, bodies, gas


def print_info(gravity_initial_total_energy, gravity, hydro, gas, i, start_mass, bodies, t):
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


# print(dir(gas))
def gravity_hydro_bridge(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t_end, dt, dt_bridge, n_stars):
    dt_SN = 0.01 | units.Myr
    t_steps = np.arange(0, t_end.value_in(units.Myr), dt.value_in(units.Myr)) | units.Myr
    # onestepplot()

    Us = []
    Ks = []
    Ts = []

    oldpos = []
    oldvel = []

    for i, t in enumerate(tqdm(t_steps)):
        gravity, hydro, gravhydro, evolution, wind, bodies, gas = simulate(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t)

        U = bodies.potential_energy() + gas.potential_energy()
        K = bodies.kinetic_energy() + gas.kinetic_energy()
        Us.append(abs(U.value_in(units.m**2 * units.kg * units.s**-2)))
        Ks.append(abs(K.value_in(units.m**2 * units.kg * units.s**-2)))
        Ts.append(t.value_in(units.Myr))

        oldpos.append(gas.position)
        oldvel.append(gas.velocity)

        # print_info(gravity_initial_total_energy, gravity, hydro, gas, i, start_mass, bodies, t)

        if star_control(bodies, n_stars) == 4:
            t_SN = t
            print(t_SN)
            hydro.parameters.timestep = dt_SN / 2
            gravhydro.timestep = dt_SN/4
            # plot_energies(Us, Ks, Ts)
            # onestepplot()
            break
    
    t_steps_supernova = np.arange(t_SN.value_in(units.Myr), t_SN.value_in(units.Myr) + 3, dt_SN.value_in(units.Myr)) | units.Myr
    for i, t in enumerate(tqdm(t_steps_supernova)):
        gravity, hydro, gravhydro, evolution, wind, bodies, gas = simulate(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t)

        
        U = bodies.potential_energy() + gas.potential_energy()
        K = bodies.kinetic_energy() + gas.kinetic_energy()
        Us.append(abs(U.value_in(units.m**2 * units.kg * units.s**-2)))
        Ks.append(abs(K.value_in(units.m**2 * units.kg * units.s**-2)))
        Ts.append(t.value_in(units.Myr))

        oldpos.append(gas.position)
        oldvel.append(gas.velocity)

        # print_info(gravity_initial_total_energy, gravity, hydro, gas, i, start_mass, bodies, t)


    
    plot_energies(Us, Ks, Ts)
    # onestepplot()
    

    evolution.stop()
    gravity.stop()
    hydro.stop()

    np.save("./positions.npy", oldpos)
    np.save("./velocities.npy", oldvel)
    np.save("./times.npy", Ts)
    

t_end = 10.0 | units.Myr
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

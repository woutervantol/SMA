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
            }
            
# Stellar evolution
evolution = SeBa()
evolution.particles.add_particles(bodies)
ch_e2g = evolution.particles.new_channel_to(gravity.particles)
ch_g2e = gravity.particles.new_channel_to(bodies) # Unnecessary?
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
    dt = 0.2|units.Myr  #1.0*Pinner

    t_steps = np.arange(model_time.value_in(units.Myr), t_end.value_in(units.Myr), dt.value_in(units.Myr)) | units.Myr
    t_steps_coderen = np.concatenate((t_steps[:1], t_steps[28:]), axis=None)    # Voor snelheid coderen
    print(t_steps_coderen)                                                      #
    print("t_steps:", t_steps)
    fig, ax = plt.subplots(3, 3)
    fig.suptitle('Cluster at initialization')
    fig2, ax2 = plt.subplots(3, 3)
    fig2.suptitle('Cluster during core helium burning') # Kan wsl wel weg.
    fig3, ax3 = plt.subplots(3, 3)
    fig3.suptitle('Cluster right after supernova')
    fig4, ax4 = plt.subplots(3, 3)
    fig4.suptitle('Cluster long after first supernova')
    ax = ax.flatten()
    ax2 = ax2.flatten()
    ax3 = ax3.flatten()
    i_dissapear = False
    i_BH = False
    for i, t in enumerate(tqdm(t_steps_coderen)):
        dE_gravity = gravity_initial_total_energy/(gravity.get_total_energy()+hydro.get_total_energy())
        print("dE:", dE_gravity, "; t=", t)
        evolution.evolve_model(t)
        channel_to_wind.copy()      # wind with hydro and grav: Book 8.1.1 p.323
        wind.evolve_model(t)               
        ch_e2g.copy()
        ch_e2b.copy()
        # Hier de supernova check
        # 1 = MS, 4 = He core burning (duurt +-1 Myr), 14 = black hole.        
        bodies_pd = pd.DataFrame(np.array(bodies.stellar_type.number), columns=["stellar_type"])
        if int(bodies_pd.value_counts().loc[1]) < 1000:
            if i_dissapear == False:
                i_dissapear = i
            print("A main sequence star has evolved!")
            ax2[i-i_dissapear].scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), s=1)
            ax2[i-i_dissapear].scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), s=1) #,c=np.log(m_stars.value_in(units.MSun)))
            ax2[i-i_dissapear].scatter(bodies[np.argmax(bodies.mass)].x.value_in(units.parsec), bodies[np.argmax(bodies.mass)].y.value_in(units.parsec), s=5, c="red")
#            ax2[i].set_title('t = %f Myr' %t)
        
            if list(bodies_pd.value_counts().index[1]) == [4]: 
                print("Core helium burning stars:", int(bodies_pd.value_counts().loc[4]))
                # Zouden hier al kleinere stappen kunnen nemen om dichter bij de echte supernova te komen.

            if list(bodies_pd.value_counts().index[1]) == [14]:
                if i_BH == False:
                    i_BH = i
                print("Black holes:", int(bodies_pd.value_counts().loc[14]))
                dt_supernova = 0.001 | units.Myr
                # dt van stellar_winds moet ook nog aangepast worden! ##############################################
                t_steps_supernova = np.arange(t.value_in(units.Myr), t_end.value_in(units.Myr), dt_supernova.value_in(units.Myr)) | units.Myr
                print(t_steps_supernova)
            
                gravhydro.evolve_model(t)
                channel["to_stars"].copy()
                channel["to_gas"].copy()
                
                for i_SN, t in enumerate(tqdm(t_steps_supernova)):
                    dE_gravity = gravity_initial_total_energy/(gravity.get_total_energy()+hydro.get_total_energy())
                    print("dE:", dE_gravity, "; t=", t)
                    evolution.evolve_model(t)
                    channel_to_wind.copy()      # wind with hydro and grav: Book 8.1.1 p.323
                    wind.evolve_model(t)
                    ch_e2g.copy()
                    ch_e2b.copy()
                    bodies_pd = pd.DataFrame(np.array(bodies.stellar_type.number), columns=["stellar_type"])
                    print("\n", bodies_pd.value_counts(), "\n")
                    gravhydro.evolve_model(t)
                    channel["to_stars"].copy()
                    channel["to_gas"].copy()
                    
                    current_gasmass = np.sum(gas.mass)
                    print("Total mass of gas:", current_gasmass)
                    current_gasnumber = current_gasmass/mgas
                    print("# of gass particles:", current_gasnumber)
                    
                    if i_SN < 9:
                        ax3[i_SN].scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), s=1)
                        ax3[i_SN].scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), s=1) #,c=np.log(m_stars.value_in(units.MSun)))
                        ax3[i_SN].scatter(bodies[np.argmax(bodies.mass)].x.value_in(units.parsec), bodies[np.argmax(bodies.mass)].y.value_in(units.parsec), s=5, c="red")
#                        ax[i].set_title('t = %i Myr' %t)
                    if i_SN == 9:
                        print("Showing supernova picture. Will continue after closing...")
                        plt.show()
                    if i_SN > 21:
                        ax4[i_SN].scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), s=1)
                        ax4[i_SN].scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), s=1) #,c=np.log(m_stars.value_in(units.MSun)))
                        ax4[i_SN].scatter(bodies[np.argmax(bodies.mass)].x.value_in(units.parsec), bodies[np.argmax(bodies.mass)].y.value_in(units.parsec), s=5, c="red")
#                        ax[i].set_title('t = %i Myr' %t)
                    if i_SN == 200:
                        break
                break
            
        gravhydro.evolve_model(t)
        channel["to_stars"].copy()
        channel["to_gas"].copy()

        current_gasmass = np.sum(gas.mass)
        print("Total mass of gas:", current_gasmass)
        current_gasnumber = current_gasmass/mgas
        print("# of gass particles:", current_gasnumber)
        
        # print("gravitational energy: ", bodies.potential_energy())
        # print("kinetic energy: ", bodies.kinetic_energy())
        print("-Ep/Ek:", - bodies.potential_energy() / bodies.kinetic_energy())
        print("Total mass:", np.sum(bodies.mass) | units.MSun)
        if i < 9:#1:# and i%10==0:
            ax[i].scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), s=1)
            ax[i].scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), s=1)#, c=np.log(m_stars.value_in(units.MSun)))
            ax[i].scatter(bodies[np.argmax(bodies.mass)].x.value_in(units.parsec), bodies[np.argmax(bodies.mass)].y.value_in(units.parsec), s=5, c="red")
#            ax[i].set_title('t = %f Myr' %t)
        
        
        bodies_pd = pd.DataFrame(np.array(bodies.stellar_type.number), columns=["stellar_type"])
        print(int(bodies_pd.value_counts().loc[1]))
        print("\n", bodies_pd.value_counts(), "\n")
            
        #if i == 8:     #Weggecomment om een supernova te vinden.
            #break
        if (i_dissapear != False) & (i-i_dissapear == 8):  # Dan zijn er 9 stappen na het verdwijnen van een main sequence ster gebeurd.
            break
    plt.show()
    
    evolution.stop()
    gravity.stop()
    hydro.stop()

t_end = 30.0 | units.Myr
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
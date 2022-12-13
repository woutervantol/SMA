#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Main simulation script.

@author: Rick Dullaart, Rutger Rijnenberg & Wouter van Tol
'''

from os import getcwd, chdir
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm
from amuse.community.fi.interface import Fi
from amuse.community.fractalcluster.interface import new_fractal_cluster_model
from amuse.community.ph4.interface import ph4
from amuse.community.seba.interface import SeBa
from amuse.couple import bridge
from amuse.ext.salpeter import new_salpeter_mass_distribution
from amuse.ext.stellar_wind import new_stellar_wind
from amuse.lab import new_plummer_gas_model
from amuse.units import (units, nbody_system)
from amuse.ext.molecular_cloud import molecular_cloud

# np.random.seed(1)




def fix_cwd():
    """Make sure the script is running from the main project folder."""
    folder = getcwd().split("/")[-1]
    if folder == "src":
        chdir("../")
    folder = getcwd().split("/")[-1]
    if not folder == "SMA":
        print("Please execute this script from the main project folder.")
        sys.exit(1)


class Clustersimulation:
    """Simulation class."""

    def __init__(self, gasmass, run):
        self.load = False
        self.evostars = False
        self.gasmass = gasmass
        self.run = run
        self.dt = 0.1 | units.Myr
        self.dt_winds = 0.05 | units.Myr
        self.dt_hydro = 0.04 | units.Myr
        self.dt_bridge = 0.02 | units.Myr  #1.0*Pinner
        # self.dt_SN = 0.01 | units.Myr
        self.n_stars = 1000
        self.n_gas = 10000
        self.r_cluster = 1.0 | units.parsec
        self.t_end = 50.0 | units.Myr
        self.settle_time = self.t_end - (25.0|units.Myr)

        gravconverter = self.init_stars()
        ## Maak gasdeeltjes
        gasconverter=nbody_system.nbody_to_si(1|units.kg, self.r_cluster*2)
        if not self.load:
            self.gas = self.create_swiss_cheese_gas(gasconverter)
        else:
            self.gas = np.load("./data/initgas.npy", allow_pickle=True)[0]
        # print(type(self.gas))
        # self.gas.remove(self.gas)
        # print(self.gas)
        # # self.lastgas = self.gas[-1]

        plt.scatter(self.gas.x.value_in(units.parsec), self.gas.y.value_in(units.parsec), s=0.5, label="Gas", alpha=0.5)
        plt.scatter(self.bodies.x.value_in(units.parsec), self.bodies.y.value_in(units.parsec), s=3, label="Stars")
        plt.xlabel("Parsec")
        plt.ylabel("Parsec")
        plt.title("Fractal example")
        plt.show()

        if not self.load:
            np.save("./data/initgas.npy", [self.gas])
            np.save("./data/initbodies.npy", [self.bodies])
            addasdas
            

        
        self.hydrocode(gravconverter)
        self.gravhydro = self.bridge()
        self.stellar_evolution()
        self.stellar_wind()

        self.channel = {
            "from_stars": self.bodies.new_channel_to(self.gravity.particles),
            "to_stars": self.gravity.particles.new_channel_to(self.bodies),
            "from_gas": self.gas.new_channel_to(self.hydro.particles),
            "to_gas": self.hydro.particles.new_channel_to(self.gas),
            "evo_to_grav": self.evolution.particles.new_channel_to(self.gravity.particles),
            "evo_to_stars": self.evolution.particles.new_channel_to(self.bodies),
            "stars_to_wind": self.bodies.new_channel_to(self.wind.particles)
        }

    def hydrocode(self, gravconverter):
        '''Create a hydro code and a gas distribution and put the gas in the hydro code.'''
        self.hydro = Fi(gravconverter, mode='g6lib')
        self.hydro.parameters.gamma = 1
        self.hydro.parameters.isothermal_flag = True
        self.hydro.parameters.integrate_entropy_flag = False
        self.hydro.parameters.timestep = self.dt_hydro
        self.hydro.parameters.verbosity = 0
        self.hydro.parameters.eps_is_h_flag = False    # h_smooth is constant
        eps = 0.1 | units.au
        self.hydro.parameters.gas_epsilon = eps
        self.hydro.parameters.sph_h_const = eps
        self.hydro.particles.add_particles(self.gas)

    def bridge(self):
        '''Bridge the codes.'''
        gravhydro = bridge.Bridge(use_threading=False) #, method=SPLIT_4TH_S_M4)
        gravhydro.add_system(self.gravity, (self.hydro,))
        gravhydro.add_system(self.hydro, (self.gravity,))
        gravhydro.timestep = self.dt_bridge # min 2x de output timestep
        return gravhydro

    def stellar_evolution(self):
        ''' Create stellar evolution.'''
        self.evolution = SeBa()
        self.evolution.particles.add_particles(self.bodies)


    def stellar_wind(self):
        '''Create Stellar wind.'''
        # p. 223
        mgas = np.sum(self.gas.mass)/self.n_gas
        self.wind = new_stellar_wind(
            mgas/100, target_gas=self.gas, timestep=self.dt_winds, derive_from_evolution=True
        )
        self.wind.particles.add_particles(self.bodies)

    def init_stars(self, alpha_IMF=-2.35):
        '''Create stars with masses, positions and velocities and put them in the ph4 module.'''
        while True:
            m_stars = new_salpeter_mass_distribution(
                self.n_stars, 0.1|units.MSun, 100|units.MSun, alpha_IMF
            )
            masslist = np.sort(m_stars)
            heaviest_stars = [masslist[-1].value_in(units.MSun), masslist[-2].value_in(units.MSun)]

            if heaviest_stars[0] > 27 and heaviest_stars[0] < 33 and heaviest_stars[1] < 27:
                break
        print("max mass star:", m_stars[np.argmax(m_stars)])
        if not self.load:
            self.total_mass = np.sum(m_stars)
            gravconverter=nbody_system.nbody_to_si(self.total_mass, self.r_cluster)
            self.bodies=new_fractal_cluster_model(self.n_stars, fractal_dimension=1.6, convert_nbody=gravconverter)
            self.bodies.scale_to_standard(gravconverter)
            self.bodies.mass = m_stars
        else:
            self.bodies = np.load("./data/initbodies.npy", allow_pickle=True)[0]
            self.total_mass = np.sum(self.bodies.mass)
            gravconverter=nbody_system.nbody_to_si(self.total_mass, self.r_cluster)
        
        
        self.gravity = ph4(gravconverter)
        self.gravity.particles.add_particles(self.bodies)
        return gravconverter

    def create_swiss_cheese_gas(self, gasconverter):
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
        # gas = new_plummer_gas_model(self.n_gas, convert_nbody=gasconverter)
        gas = molecular_cloud(targetN=self.n_gas, convert_nbody=gasconverter).result

        tryradii = np.logspace(np.log10(0.001), np.log10(1), 50) | units.parsec
        sorted_stars = np.sort(self.bodies.mass)
        for m_star in sorted_stars:
            star = self.bodies[self.bodies.mass == m_star]
            for radius in tryradii:
                gas_to_remove = gas[(star.position - gas.position).lengths()<radius]
                if np.sum(gas_to_remove.mass) > star.mass:
                    gas.remove_particle(gas_to_remove)
                    break
        return gas

    def gravity_hydro_bridge(self):
        '''Run the simulation.'''
        t_steps = np.arange(0, self.t_end.value_in(units.Myr), self.dt.value_in(units.Myr)) | units.Myr

        original_gas = self.gas.copy()
        
        gaslist = []
        bodieslist = []
        timeslist = []
        gas_indices = []

        last_gascount = len(self.gas)
        gascount_at_SN = len(self.gas)

        for timestamp in tqdm(t_steps):
            
            self.simulate_timestamp(timestamp)

            for i in range(1, len(self.gas)):
                if original_gas[-i].key in self.gas.key:
                    gas_indices.append(np.argwhere(self.gas.key == original_gas[-i].key)[0][0])
                    break

            gaslist.append(self.gas.copy())
            bodieslist.append(self.bodies.copy())
            timeslist.append(timestamp)


            #als het aantal deeltjes stijgt totdat het onder het originele niveau
            if len(self.gas) > last_gascount or len(self.gas) > gascount_at_SN:
                # print("star type: {}, time: {}".format(star_control(self.bodies), timestamp))
                self.hydro.parameters.timestep = 0.001 | units.Myr
                self.gravhydro.timestep = 0.0005 | units.Myr
            else:
                # print("star type: {}, time: {}".format(star_control(self.bodies), timestamp))
                self.hydro.parameters.timestep = 0.04 | units.Myr
                self.gravhydro.timestep = 0.02 | units.Myr
                gascount_at_SN = len(self.gas)

            if timestamp > self.settle_time:
                self.evostars = True

            last_gascount = len(self.gas)

        self.evolution.stop()
        self.gravity.stop()
        self.hydro.stop()

        filestring = "_ratio{}_run{}".format(self.gasmass, self.run)

        np.save("./data/gas{}.npy".format(filestring), gaslist)
        np.save("./data/bodies{}.npy".format(filestring), bodieslist)
        np.save("./data/times{}.npy".format(filestring), timeslist)
        np.save("./data/gas_indices{}.npy".format(filestring), gas_indices)


    def simulate_timestamp(self, timestamp):
        '''Run the simulating to the given timestamp.'''
        if self.evostars:
            self.evolution.evolve_model(timestamp - self.settle_time)
        self.channel["evo_to_grav"].copy()
        self.channel["evo_to_stars"].copy()
        self.channel["stars_to_wind"].copy()      # wind with hydro and grav: Book 8.1.1 p.323
        self.wind.evolve_model(timestamp)
        self.gas = delete_outofbounds(self.gas)
        self.gas.synchronize_to(self.hydro.particles)
        self.gravhydro.evolve_model(timestamp)
        self.channel["to_stars"].copy()
        self.channel["to_gas"].copy()



def star_control(bodies):
    bodies_pd = pd.DataFrame(np.array(bodies.stellar_type.number), columns=["stellar_type"])
    return list(bodies_pd.value_counts().index[-1])[0]

def delete_outofbounds(gas):
    '''Delete all the gas that is out of bounds.'''
    mask = gas.position.lengths() >= 1e18|units.m #1e18m is 30 parsec
    selection = gas[mask]
    gas.remove_particles(selection)
    return gas


def print_info(gravity_initial_total_energy, gravity, hydro, gas, i, start_mass, bodies, t,):
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


def main(gasmass, run):
    '''Run the simulation with standard parameters.'''
    my_simulation = Clustersimulation(gasmass, run)
    my_simulation.gravity_hydro_bridge()
    return 0


if __name__ == "__main__":
    fix_cwd()
    runs_per_gasmass = np.arange(0, 10)
    gas_mass_ratios = [5]
    for gasmass in gas_mass_ratios:
        for run in runs_per_gasmass:
            main(gasmass, "{}".format(run))
    # exit(main(5, "finalcomp2"))


# Interessant boek? https://misaladino.com/wp-content/uploads/2019/11/Thesis_Martha_Irene.pdf

#### Stellar types ####
# "deeply or fully convective low mass MS star" #  0
# "Main Sequence star"                          #  1
# "Hertzsprung Gap"                             #  2
# "First Giant Branch"                          #  3
# "Core Helium Burning"                         #  4
# "First Asymptotic Giant Branch"               #  5
# "Second Asymptotic Giant Branch"              #  6
# "Main Sequence Naked Helium star"             #  7
# "Hertzsprung Gap Naked Helium star"           #  8
# "Giant Branch Naked Helium star"              #  9
# "Helium White Dwarf"                          # 10
# "Carbon/Oxygen White Dwarf"                   # 11
# "Oxygen/Neon White Dwarf"                     # 12
# "Neutron Star"                                # 13
# "Black Hole"                                  # 14
# "Massless Supernova"                          # 15
# "Unknown stellar type"                        # 16
# "Pre-main-sequence Star"                      # 17
# "Planet"                                      # 18
# "Brown Dwarf"                                 # 19

# 2 timesteps:
# De bridge timestep (minimaal half van de output timestep)
# De output timestep. Hoe vaak wil je wind creeren? Stellar evolution heeft zn
# eigen timestep en die kan je in principe aanvragen.

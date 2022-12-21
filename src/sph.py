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
from amuse.units import (units, nbody_system, constants)
from amuse.ext.molecular_cloud import molecular_cloud

#file management
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
    def __init__(self, gasmass, run, load=True):
        self.load = load #set False to create initial system or True to load last used system
        self.gasmass = gasmass #ratio between mass in gas and mass in stars
        self.run = run
        self.dt = 0.002 | units.Myr
        self.dt_winds = 0.002 | units.Myr
        self.dt_hydro = 0.001 | units.Myr
        self.dt_bridge = 0.002 | units.Myr  #1.0*Pinner
        self.supernovagasmass = (1|units.MSun)/10
        self.windsgasmass = (1|units.MSun)/10000
        self.n_stars = 10
        self.n_gas = 100
        self.r_cluster = 1.0 | units.parsec
        self.t_end = 10 | units.Myr
        self.settle_time = self.t_end - (10|units.Myr)

        #Create stars and gas
        gravconverter = self.init_stars()
        gasconverter = nbody_system.nbody_to_si((gasmass+1)*self.total_mass, self.r_cluster*2)
        if not self.load:
            self.gas = self.create_swiss_cheese_gas(gasconverter)
        else:
            self.gas = np.load("./data/initgas.npy", allow_pickle=True)[0]
        self.initialised_gas = self.gas.key.copy() #used to keep track of newly formed gas, so we can change their position and velocity by hand

        #Create new initial system and close program
        if not self.load:
            np.save("./data/initgas.npy", [self.gas])
            np.save("./data/initbodies.npy", [self.bodies])
            return
        
        #initialise all codes
        self.hydrocode(gravconverter)
        self.stellar_evolution()
        self.stellar_wind()
        self.gravhydro = self.bridge()

        #Set all needed channels
        self.channel = {
            "to_stars": self.gravity.particles.new_channel_to(self.bodies),
            "to_gas": self.hydro.particles.new_channel_to(self.gas),
            "evo_to_grav": self.evolution.particles.new_channel_to(self.gravity.particles),
            "evo_to_stars": self.evolution.particles.new_channel_to(self.bodies),
            "stars_to_wind": self.bodies.new_channel_to(self.wind.particles)
        }

    def hydrocode(self, gravconverter):
        #Create a hydro code and a gas distribution and put the gas in the hydro code. Parameter values are taken from tutorial
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
        #Bridge gravity and hydro codes
        gravhydro = bridge.Bridge(use_threading=False)
        gravhydro.add_system(self.gravity, (self.hydro,))
        gravhydro.add_system(self.hydro, (self.gravity,))
        gravhydro.timestep = self.dt_bridge
        return gravhydro

    def stellar_evolution(self):
        #Start SeBa code
        self.evolution = SeBa()
        self.evolution.particles.add_particles(self.bodies)


    def stellar_wind(self):
        # Start stellar wind code
        # See p. 223
        self.wind = new_stellar_wind(self.windsgasmass, target_gas=self.gas, timestep=self.dt_winds, derive_from_evolution=True)
        self.wind.particles.add_particles(self.bodies)

    def init_stars(self, alpha_IMF=-2.35):
        #Create a mass Salpeter mass distribution. If it doesn't fit our criteria, reroll the distribution until it does (largest star between 27 and 33 MSun and second largest smaller)
        while True:
            m_stars = new_salpeter_mass_distribution(
                self.n_stars, 0.1|units.MSun, 100|units.MSun, alpha_IMF
            )
            masslist = np.sort(m_stars)
            heaviest_stars = [masslist[-1].value_in(units.MSun), masslist[-2].value_in(units.MSun)]

            if heaviest_stars[0] > 27 and heaviest_stars[0] < 33 and heaviest_stars[1] < 27:
                break
        
        #Create stars or load stars and put them in the ph4 module. 
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
        
        #Start ph4 and return the converter since the hydro code needs it
        self.gravity = ph4(gravconverter)
        self.gravity.particles.add_particles(self.bodies)
        return gravconverter

    def create_swiss_cheese_gas(self, gasconverter):
        # Create uniform cloud of gas with random velocities. Remove gas around the stars with mass similar to the mass of the star.
        gas = molecular_cloud(targetN=self.n_gas, convert_nbody=gasconverter).result

        # Try for a number of radii if the mass of the gas within this radius is similar to the mass of the star
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
        #Run simulation
        t_steps = np.arange(0, self.t_end.value_in(units.Myr), self.dt.value_in(units.Myr)) | units.Myr

        original_gas = self.gas.copy() #used to distinguish the border between original gas and stellar wind gas throughout the simulation
        
        gaslist = []
        bodieslist = []
        timeslist = []
        gas_indices = []

        for timestamp in tqdm(t_steps):
            # Evolve all codes
            self.simulate_timestamp(timestamp)

            # Find border between original gas and stellar wind gas
            for i in range(1, len(self.gas)):
                if original_gas[-i].key in self.gas.key:
                    gas_indices.append(np.argwhere(self.gas.key == original_gas[-i].key)[0][0])
                    break

            # Save data
            gaslist.append(self.gas.copy())
            bodieslist.append(self.bodies.copy())
            timeslist.append(timestamp)

        self.evolution.stop()
        self.gravity.stop()
        self.hydro.stop()

        filestring = "_ratio{}_run{}".format(self.gasmass, self.run)
        np.save("./data/gas{}.npy".format(filestring), gaslist)
        np.save("./data/bodies{}.npy".format(filestring), bodieslist)
        np.save("./data/times{}.npy".format(filestring), timeslist)
        np.save("./data/gas_indices{}.npy".format(filestring), gas_indices)


    def simulate_timestamp(self, timestamp):
        # Evolve all codes
        if timestamp >= self.settle_time:
            self.evolution.evolve_model(timestamp - self.settle_time) # start evolving after settle time
        self.channel["evo_to_grav"].copy()
        self.channel["evo_to_stars"].copy()
        self.channel["stars_to_wind"].copy()      # wind with hydro and grav: Book 8.1.1 p.323
        
        # From supernova onwards we want to increase the mass per particle or else we get millions of particles at the supernova
        if star_control(self.bodies) == 4 or star_control(self.bodies) == 14 or star_control(self.bodies) == 5:
            self.wind.sph_particle_mass = self.supernovagasmass
        else:
            self.wind.sph_particle_mass = self.windsgasmass
        # self.wind.evolve_model(timestamp)
        print(len(self.gas))
        self.wind_init()
        self.gas = delete_outofbounds(self.gas)
        self.gas.synchronize_to(self.hydro.particles)
        self.gravhydro.evolve_model(timestamp)
        self.channel["to_stars"].copy()
        self.channel["to_gas"].copy()
    
    def wind_init(self):
        # Check which gas particles are newly created and initialise them by hand. We place them close to the star and change their velocity
        bigstar = self.bodies[np.argmax(self.bodies.mass)]
        for i in range(1, len(self.gas)):
            if self.gas[-i].key not in self.initialised_gas:
                self.gas[-i].position = bigstar.position + self.gas[-i].position/self.gas[-i].position.length() * (0.001|units.parsec)
                self.gas[-i].velocity = self.gas[-i].velocity/self.gas[-i].velocity.length() * (100|units.kms)
                self.initialised_gas = np.append(self.initialised_gas, self.gas[-i].key)
            else:
                break



def star_control(bodies):
    # Return the evolutionary phase of the largest star. Used to see when it goes supernova
    bodies_pd = pd.DataFrame(np.array(bodies.stellar_type.number), columns=["stellar_type"])
    return list(bodies_pd.value_counts().index[-1])[0]

def delete_outofbounds(gas):
    # Delete all the gas that is out of bounds.
    mask = gas.position.lengths() >= 1e18|units.m #1e18m is around 30 parsec
    selection = gas[mask]
    gas.remove_particles(selection)
    return gas


def print_info(gravity_initial_total_energy, gravity, hydro, gas, i, start_mass, bodies, t,):
    # Used for debugging
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


def main(gasmass, run, load=True):
    # Run the simulation with standard parameters.
    my_simulation = Clustersimulation(gasmass, run, load)
    if load:
        my_simulation.gravity_hydro_bridge()
    return 0


if __name__ == "__main__":
    fix_cwd()
    main(5, "sketchy_final_comp", False)
    exit(main(5, "sketchy_final_comp"))

# 2 timesteps:
# De bridge timestep (minimaal half van de output timestep)
# De output timestep. Hoe vaak wil je wind creeren? Stellar evolution heeft zn
# eigen timestep en die kan je in principe aanvragen.

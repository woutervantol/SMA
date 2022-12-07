#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from os import getcwd, chdir
from sph import *
from matplotlib import animation

fig, ax = plt.subplots(2)

dt_SN = 0.01 | units.Myr

model_time = 0 | units.Myr
t_steps = np.arange(model_time.value_in(units.Myr), t_end.value_in(units.Myr), dt.value_in(units.Myr)) | units.Myr

def makeplot(time, Us, Ks, Ts):
    ax[0].cla()
    ax[0].set_title("{:.2f} Myr".format(time.value_in(units.Myr)))
    ax[0].set_xlim(-5, 5)
    ax[0].set_ylim(-5, 5)
    ax[0].set_xlabel("parsec")
    ax[0].set_ylabel("parsec")
    ax[0].scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), s=0.5, label="gas")
    ax[0].scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), s=3, label="stars")
    ax[0].scatter(bodies[np.argmax(bodies.mass)].x.value_in(units.parsec), bodies[np.argmax(bodies.mass)].y.value_in(units.parsec), s=5, c="red")

    ax[1].cla()
    ax[1].plot(Ts, Us, label="potential energy")
    ax[1].plot(Ts, Ks, label="kinetic energy")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("Time (Myr)")
    ax[1].set_ylabel("Energy (m^2 kg s^-2)")

def fix_cwd():
    folder = getcwd().split("/")[-1]
    if folder == "src":
        chdir("../")
    folder = getcwd().split("/")[-1]
    if not folder == "SMA":
        print("Please execute this script from the main project folder.")
        exit(1)


def calc_time():
    model_time = 0 | units.Myr
    while star_control(bodies, n_stars) != 4:
        model_time = np.round((model_time+dt).value_in(units.Myr), 4) | units.Myr
        yield model_time, dt
    endtime = model_time + (3|units.Myr)
    while model_time < endtime:
        print(star_control(bodies, n_stars))
        model_time = np.round((model_time+dt_SN).value_in(units.Myr), 6) | units.Myr
        yield model_time, dt_SN

Us = []
Ks = []
Ts = []

def update(times, gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t_end, dt, dt_bridge, n_stars):
    time = times[0]
    print(time)
    timestep = times[1]
    U = bodies.potential_energy() + gas.potential_energy()
    K = bodies.kinetic_energy() + gas.kinetic_energy()
    Us.append(abs(U.value_in(units.m**2 * units.kg * units.s**-2)))
    Ks.append(abs(K.value_in(units.m**2 * units.kg * units.s**-2)))
    Ts.append(time.value_in(units.Myr))
    if timestep == dt:
        gravity, hydro, gravhydro, evolution, wind, bodies, gas = simulate(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, time)
        makeplot(time, Us, Ks, Ts)
    else:
        hydro.parameters.timestep = dt_SN / 2
        gravhydro.timestep = dt_SN * 2
        gravity, hydro, gravhydro, evolution, wind, bodies, gas = simulate(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, time)
        makeplot(time, Us, Ks, Ts)

# gravity, hydro, gravhydro, evolution, wind, bodies, gas = simulate(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, 6.400999999999994 | units.Myr)
anim = animation.FuncAnimation(fig, update, frames=calc_time, fargs=(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t_end, dt, dt_bridge, n_stars), save_count=1000)
plt.legend()
writer = animation.FFMpegWriter(fps=len(t_steps)/6.)
fix_cwd()
anim.save("./figures/animation.mp4", writer=writer)
# plt.show()
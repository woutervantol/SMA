import numpy as np
import matplotlib.pyplot as plt
from amuse.units import units
from matplotlib import animation

times = np.load("./data/times_ratio1_run0.npy", allow_pickle=True)
avg_pos = np.zeros(len(test))
avg_vel = np.zeros(len(test))
runs_per_gasmass = np.arange(0, 50)
    gas_mass_ratios = [1, 5, 10, 15, 20, 25]
    for gasmass in gas_mass_ratios:
        for run in runs_per_gasmass:
            filestring = "_ratio{}_run{}".format(gasmass, run)

            gas = np.load("./data/gas{}.npy".format(filestring), allow_pickle=True)
            bodies = np.load("./data/bodies{}.npy".format(filestring), allow_pickle=True)
            # times = np.load("./data/times{}.npy".format(filestring), allow_pickle=True)
            gas_indices = np.load("./data/gas_indices{}.npy".format(filestring), allow_pickle=True)

            for i in range(len(times)):
                avg_pos[i] += np.average(gas[i].position.lengths())
                avg_vel[i] += np.average(gas[i].velocity.lengths())
avg_pos /= 50*6
avg_vel /= 50*6

plt.plot(times, avg_pos)
plt.xlabel("Time (Myr)")
plt.ylabel("Average distance of gas from center")
plt.show()

plt.plot(times, avg_vel)
plt.xlabel("Time (Myr)")
plt.ylabel("Average velocity of gas")
plt.show()












filestring = "_ratio{}_run{}".format(5, "test")

gas = np.load("./data/gas{}.npy".format(filestring), allow_pickle=True)
bodies = np.load("./data/bodies{}.npy".format(filestring), allow_pickle=True)
times = np.load("./data/times{}.npy".format(filestring), allow_pickle=True)
gas_indices = np.load("./data/gas_indices{}.npy".format(filestring), allow_pickle=True)
for i in range(len(times)):
    times[i] = times[i].value_in(units.Myr)

#find time at which supernova starts
SNtime = 0 | units.Myr
for i in range(1, len(times)):
    if len(gas[i]) > len(gas[i-1]):
        SNtime = times[i]
        break
print(SNtime)

def plot_avg_velocities():
    plt.figure()
    vels = []
    for t in range(len(times)):
        vel = np.sqrt(gas[t].velocity[:gas_indices[i], 0].value_in(units.ms)**2 + gas[t].velocity[:gas_indices[i], 1].value_in(units.ms)**2 + gas[t].velocity[:gas_indices[i], 2].value_in(units.ms)**2)
        vels.append(np.average(vel))
    plt.title("Gas velocity")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Average gas velocity (m/s)")
    plt.vlines(SNtime, np.min(vels), np.max(vels), color="black", linestyle="dashed")
    plt.plot(times, vels)
    plt.savefig("./figures/avg_vel.png")

def plot_avg_positions():
    plt.figure()
    poss = []
    for t in range(len(times)):
        pos = np.sqrt(gas[t].position[:gas_indices[i], 0].value_in(units.parsec)**2 + gas[t].position[:gas_indices[i], 1].value_in(units.parsec)**2 + gas[t].position[:gas_indices[i], 2].value_in(units.parsec)**2)
        poss.append(np.average(pos))
    plt.title("Gas position")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Average gas position (parsec)")
    plt.vlines(SNtime, np.min(poss), np.max(poss), color="black", linestyle="dashed")
    plt.plot(times, poss)
    plt.savefig("./figures/avg_pos.png")

def plot_energies():
    Us = []
    Ks = []
    for t in range(len(times)):
        print(t)
        Us.append(-gas[t][:gas_indices[t]].potential_energy().value_in(units.m**2 * units.kg * units.s**-2) - bodies[t].potential_energy().value_in(units.m**2 * units.kg * units.s**-2))
        Ks.append(gas[t][:gas_indices[t]].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2) + bodies[t].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2))
    plt.figure()
    plt.semilogy(times, Us, label="potential energy")
    plt.semilogy(times, Ks, label="kinetic energy")
    plt.title("Energy")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Energy (m^2 kg s^-2)")
    plt.vlines(SNtime, np.min(Ks), np.max(Ks), color="black", linestyle="dashed")
    plt.legend()
    plt.savefig("./figures/energy.png")

    plt.figure()
    data = []
    for t in range(len(times)):
        data.append(Ks[t]/Us[t])
    plt.semilogy(times, data)
    plt.vlines(SNtime, np.min(data), np.max(data), color="black", linestyle="dashed")
    plt.title("boundness energy")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Energy ratio between kinetic and potential energy")
    plt.savefig("./figures/energy_rel.png")

def make_animation(frame):
    print(frame)
    ax.cla()
    ax.set_title("Time: {:.2f} Myr".format(times[frame]))
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.set_xlabel("parsec")
    ax.set_ylabel("parsec")
    
    ax.scatter(gas[frame].x[:gas_indices[frame]].value_in(units.parsec), gas[frame].y[:gas_indices[frame]].value_in(units.parsec), s=0.5, label="Gas")
    ax.scatter(gas[frame].x[gas_indices[frame]:].value_in(units.parsec), gas[frame].y[gas_indices[frame]:].value_in(units.parsec), s=2, label="Supernova gas", color="orange")
    ax.scatter(bodies[frame].x.value_in(units.parsec), bodies[frame].y.value_in(units.parsec), s=3, label="Stars")
    if times[frame] >= SNtime:
        ax.scatter(bodies[frame][np.argmax(bodies[frame].mass)].x.value_in(units.parsec), bodies[frame][np.argmax(bodies[frame].mass)].y.value_in(units.parsec), s=10, c="red", label="Supernova")
    ax.legend()


def make_columndensity(frame):
    print(frame)
    ax.cla()
    ax.set_title("Time: {:.2f} Myr".format(times[frame]))
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.set_xlabel("parsec")
    ax.set_ylabel("parsec")

    plt.hist2d(gas[frame].x[:gas_indices[frame]].value_in(units.parsec), gas[frame].y[:gas_indices[frame]].value_in(units.parsec), bins=(50, 50), range=((-10, 10), (-10, 10)), cmap="OrRd")
    ax.scatter(bodies[frame].x.value_in(units.parsec), bodies[frame].y.value_in(units.parsec), s=3, label="Stars")


# fig, ax = plt.subplots()
# anim = animation.FuncAnimation(fig, make_animation, frames=len(times))
# writer = animation.FFMpegWriter(fps=len(times)/20.)
# anim.save("./figures/animation.mp4", writer=writer)

fig, ax = plt.subplots()
anim = animation.FuncAnimation(fig, make_columndensity, frames=len(times))
writer = animation.FFMpegWriter(fps=len(times)/20.)
anim.save("./figures/columndensity.mp4", writer=writer)

# fig, ax = plt.subplots()
# anim = animation.FuncAnimation(fig, make_pos_hist, frames=len(times))
# writer = animation.FFMpegWriter(fps=len(gasvelocities)/20.)
# anim.save("./figures/pos_hist.mp4", writer=writer)

# fig, ax = plt.subplots()
# anim = animation.FuncAnimation(fig, make_vel_hist, frames=len(times))
# writer = animation.FFMpegWriter(fps=len(gasvelocities)/20.)
# anim.save("./figures/vel_hist.mp4", writer=writer)

plot_avg_positions()
plot_avg_velocities()
plot_energies()










# def make_vel_hist(frame):
#     print(frame)
#     ax.cla()
#     ax.set_title("Time: {:2f} Myr".format(times[frame]))
#     ax.set_xlabel("Absolute speed (m/s)")
#     ax.set_ylabel("Nr of gas particles")
#     ax.set_xlim(0, 10000)
#     ax.set_ylim(0, 200)
#     ax.hist(np.sqrt(gasvelocities[frame][:, 0].value_in(units.ms)**2 + gasvelocities[frame][:, 1].value_in(units.ms)**2 + gasvelocities[frame][:, 2].value_in(units.ms)**2), bins=100, range=(0, 10000))

# def make_pos_hist(frame):
#     print(frame)
#     ax.cla()
#     ax.set_title("Time: {:2f} Myr".format(times[frame]))
#     ax.set_xlabel("Position from center (parsec)")
#     ax.set_ylabel("Nr of gas particles")
#     ax.set_xlim(0, 10)
#     ax.set_ylim(0, 200)
#     ax.hist(np.sqrt(gaspositions[frame][:, 0].value_in(units.parsec)**2 + gaspositions[frame][:, 1].value_in(units.parsec)**2 + gaspositions[frame][:, 2].value_in(units.parsec)**2), bins=100, range=(0, 10))

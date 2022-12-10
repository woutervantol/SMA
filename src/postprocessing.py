import numpy as np
import matplotlib.pyplot as plt
from amuse.units import units
from matplotlib import animation

filestring = ""
filestring = "_ratio{}_run{}".format(5, "test")

kinetic_energies = np.load("./data/kinetic{}.npy".format(filestring), allow_pickle=True)
potential_energies = np.load("./data/potential{}.npy".format(filestring), allow_pickle=True)
times = np.load("./data/times{}.npy".format(filestring), allow_pickle=True)
gaspositions = np.load("./data/gaspositions{}.npy".format(filestring), allow_pickle=True)
gasvelocities = np.load("./data/gasvelocities{}.npy".format(filestring), allow_pickle=True)
starpositions = np.load("./data/starpositions{}.npy".format(filestring), allow_pickle=True)
starvelocities = np.load("./data/starvelocities{}.npy".format(filestring), allow_pickle=True)

#find time at which supernova starts
SNtime = 0 | units.Myr
for i in range(len(gasvelocities)):
    # print(len(gasvelocities[i]))
    if len(gasvelocities[i]) > len(gasvelocities[0]):
        SNtime = times[i]
        break

def plot_avg_velocities():
    plt.figure()
    vels = []
    for frame in range(len(gasvelocities)):
        temp = np.sqrt(gasvelocities[frame][:, 0].value_in(units.ms)**2 + gasvelocities[frame][:, 1].value_in(units.ms)**2 + gasvelocities[frame][:, 2].value_in(units.ms)**2)
        relevant = temp[temp < 4000]
        vels.append(np.average(relevant))
    plt.title("Gas velocity")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Average gas velocity (m/s)")
    plt.vlines(SNtime, np.min(vels), np.max(vels), color="black", linestyle="dashed")
    plt.plot(times, vels)
    plt.savefig("./figures/avg_vel.png")

def plot_avg_positions():
    plt.figure()
    poss = []
    for frame in range(len(gasvelocities)):
        temp = np.sqrt(gaspositions[frame][:, 0].value_in(units.parsec)**2 + gaspositions[frame][:, 1].value_in(units.parsec)**2 + gaspositions[frame][:, 2].value_in(units.parsec)**2)
        relevant = temp[temp < 10]
        poss.append(np.average(relevant))
    plt.title("Gas position")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Average gas position (m/s)")
    plt.vlines(SNtime, np.min(poss), np.max(poss), color="black", linestyle="dashed")
    plt.plot(times, poss)
    plt.savefig("./figures/avg_pos.png")

def plot_energies():
    plt.figure()
    plt.semilogy(times, potential_energies, label="potential energy")
    plt.semilogy(times, kinetic_energies, label="kinetic energy")
    plt.title("Energy")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Energy (m^2 kg s^-2)")
    plt.vlines(SNtime, np.min(kinetic_energies), np.max(kinetic_energies), color="black", linestyle="dashed")
    plt.legend()
    plt.savefig("./figures/energy.png")

def plot_relative_energy():
    plt.figure()
    data = kinetic_energies / potential_energies
    plt.semilogy(times, data)
    plt.vlines(SNtime, np.min(data), np.max(data), color="black", linestyle="dashed")
    plt.title("boundness energy")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Energy ratio")
    plt.legend()
    plt.savefig("./figures/energy_rel.png")

def make_vel_hist(frame):
    print(frame)
    ax.cla()
    ax.set_title("Time: {:2f} Myr".format(times[frame]))
    ax.set_xlabel("Absolute speed (m/s)")
    ax.set_ylabel("Nr of gas particles")
    ax.set_xlim(0, 10000)
    ax.set_ylim(0, 200)
    ax.hist(np.sqrt(gasvelocities[frame][:, 0].value_in(units.ms)**2 + gasvelocities[frame][:, 1].value_in(units.ms)**2 + gasvelocities[frame][:, 2].value_in(units.ms)**2), bins=100, range=(0, 10000))

def make_pos_hist(frame):
    print(frame)
    ax.cla()
    ax.set_title("Time: {:2f} Myr".format(times[frame]))
    ax.set_xlabel("Position from center (parsec)")
    ax.set_ylabel("Nr of gas particles")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 200)
    ax.hist(np.sqrt(gaspositions[frame][:, 0].value_in(units.parsec)**2 + gaspositions[frame][:, 1].value_in(units.parsec)**2 + gaspositions[frame][:, 2].value_in(units.parsec)**2), bins=100, range=(0, 10))

def make_animation(frame):
    print(frame)
    ax.cla()
    ax.set_title("Time: {:.2f} Myr".format(times[frame]))
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)
    ax.set_xlabel("parsec")
    ax.set_ylabel("parsec")
    
    #for some reason array has quantity elements instead of being a quantity vector so we have to iterate by hand to convert units
    starpos = np.zeros((len(starpositions[0]), 3))
    for i in range(len(starpositions[frame])):
        starpos[i,0] = starpositions[frame][i,0].value_in(units.parsec)
        starpos[i,1] = starpositions[frame][i,1].value_in(units.parsec)
        starpos[i,2] = starpositions[frame][i,2].value_in(units.parsec)
    ax.scatter(gaspositions[frame][:,0].value_in(units.parsec), gaspositions[frame][:,1].value_in(units.parsec), s=0.5, label="gas")
    ax.scatter(starpos[:,0], starpos[:,1], s=3, label="stars")
    # ax.scatter(bodies[np.argmax(bodies.mass)].x.value_in(units.parsec), bodies[np.argmax(bodies.mass)].y.value_in(units.parsec), s=5, c="red")

fig, ax = plt.subplots()
anim = animation.FuncAnimation(fig, make_animation, frames=len(times))
writer = animation.FFMpegWriter(fps=len(gasvelocities)/20.)
anim.save("./figures/animation.mp4", writer=writer)

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
plot_relative_energy()

import numpy as np
import matplotlib.pyplot as plt
from amuse.units import units
from matplotlib import animation

# times = np.load("./data/veelsims/times_ratio1_run0.npy", allow_pickle=True)
# for i in range(len(times)):
#     times[i] = times[i].value_in(units.Myr)

# runs_per_gasmass = [6, 7, 8, 12, 20, 24, 25, 27, 28]
# gas_mass_ratios = [5]
# for gasmass in gas_mass_ratios:
#     avg_pos = np.zeros(len(times))
#     avg_vel = np.zeros(len(times))
#     for run in runs_per_gasmass:
#         filestring = "_ratio{}_run{}".format(gasmass, run)

#         gas = np.load("./data/veelsims/gas{}.npy".format(filestring), allow_pickle=True)
#         bodies = np.load("./data/veelsims/bodies{}.npy".format(filestring), allow_pickle=True)
#         # times = np.load("./data/times{}.npy".format(filestring), allow_pickle=True)
#         gascomp = np.load("./data/veelsims/compgas{}.npy".format(filestring), allow_pickle=True)
#         bodiescomp = np.load("./data/veelsims/compbodies{}.npy".format(filestring), allow_pickle=True)

#         gas_indices = np.load("./data/veelsims/gas_indices{}.npy".format(filestring), allow_pickle=True)

#         SNtime = 0 | units.Myr
#         for i in range(1, len(times)):
#             print(len(gas[i]), times[i])
#             if len(gas[i]) > len(gas[i-1]):
#                 SNtime = times[i]
#                 # break
#         print(SNtime)

#         Us = []
#         Uscomp = []
#         Ks = []
#         Kscomp = []
#         for t in range(len(times)):
#             print(t)
#             #stars
#             # Us.append(-bodies[t].potential_energy().value_in(units.m**2 * units.kg * units.s**-2))
#             # Uscomp.append(-bodiescomp[t].potential_energy().value_in(units.m**2 * units.kg * units.s**-2))
#             # Ks.append(bodies[t].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2))
#             # Kscomp.append(bodiescomp[t].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2))
#             #gas
#             # Us.append(-gas[t][:gas_indices[t]].potential_energy().value_in(units.m**2 * units.kg * units.s**-2))
#             # Uscomp.append(-gascomp[t].potential_energy().value_in(units.m**2 * units.kg * units.s**-2))
#             # Ks.append(gas[t][:gas_indices[t]].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2))
#             # Kscomp.append(gascomp[t].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2))
#             #both
#             Us.append(-gas[t][:gas_indices[t]].potential_energy().value_in(units.m**2 * units.kg * units.s**-2) - bodies[t].potential_energy().value_in(units.m**2 * units.kg * units.s**-2))
#             # Uscomp.append(-gascomp[t].potential_energy().value_in(units.m**2 * units.kg * units.s**-2) - bodiescomp[t].potential_energy().value_in(units.m**2 * units.kg * units.s**-2))
#             Ks.append(gas[t][:gas_indices[t]].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2) + bodies[t].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2))
#             # Kscomp.append(gascomp[t].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2) + bodiescomp[t].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2))
#         plt.figure()
#         plt.semilogy(times, Us, label="Potential energy")
#         # plt.semilogy(times, Uscomp, label="Potential energy without winds")
#         plt.semilogy(times, Ks, label="Kinetic energy")
#         # plt.semilogy(times, Kscomp, label="Kinetic energy without winds")
#         # plt.plot(times, np.array(Ks) - np.array(Us))

#         # plt.title("Energy comparison{} {}".format(gasmass, run))
#         plt.xlabel("Time (Myr)")
#         plt.ylabel("Energy (m^2 kg s^-2)")
#         plt.axvline(SNtime, 0, 1, color="black", linestyle="dashed", label="Start supernova")
#         plt.axvline(20, 0, 1, color="gray", linestyle="dashed", label="Start SeBa")
#         plt.legend()
#         plt.show()


        # for i in range(len(times)):
        #     avg_pos[i] += np.average(gas[i].position.lengths().value_in(units.parsec))
        #     avg_vel[i] += np.average(gas[i].velocity.lengths().value_in(units.ms))
#     print(gasmass)
#     avg_pos /= len(runs_per_gasmass)
#     avg_vel /= len(runs_per_gasmass)
#     plt.figure(1)
#     plt.xlabel("Time (Myr)")
#     plt.ylabel("Average distance of gas from center (parsec)")
#     plt.plot(times, avg_pos, label="gas mass ratio: {}".format(gasmass))
#     plt.legend()
#     plt.figure(2)
#     plt.xlabel("Time (Myr)")
#     plt.ylabel("Average velocity of gas")
#     plt.plot(times, avg_vel, label="gas mass ratio: {}".format(gasmass))
#     plt.legend()
# plt.show()












filestring = "_ratio{}_run{}".format(5, "final")
filestringcomp = "_ratio{}_run{}".format(5, "finalcomp")

gas = np.load("./data/gas{}.npy".format(filestring), allow_pickle=True)
bodies = np.load("./data/bodies{}.npy".format(filestring), allow_pickle=True)
gascomp = np.load("./data/gas{}.npy".format(filestringcomp), allow_pickle=True)
bodiescomp = np.load("./data/bodies{}.npy".format(filestringcomp), allow_pickle=True)


    

times = np.load("./data/times{}.npy".format(filestring), allow_pickle=True)
gas_indices = np.load("./data/gas_indices{}.npy".format(filestring), allow_pickle=True)
for i in range(len(times)):
    times[i] = times[i].value_in(units.Myr)

#find time at which supernova starts
SNtime = 0 | units.Myr
for i in range(1, len(times)):
    # print(len(gas[i]), times[i])
    if len(gas[i]) > len(gas[i-1]):
        SNtime = times[i]
        # break
print(SNtime)

def plot_avg_velocities():
    plt.figure()
    vels = []
    velscomp = []
    for t in range(len(times)):
        vels.append(np.average(gas[t][:gas_indices[t]].velocity.length().value_in(units.ms)))
        velscomp.append(np.average(gascomp[t].velocity.length().value_in(units.ms)))
    plt.title("Gas velocity comparison")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Average gas velocity (m/s)")
    plt.axvline(SNtime, 0, 1, color="black", linestyle="dashed", label="Start supernova")
    plt.axvline(25, 0, 1, color="gray", linestyle="dashed", label="Start SeBa")
    plt.plot(times, vels, label="With wind")
    plt.plot(times, velscomp, label="Without wind")
    plt.legend()
    plt.savefig("./figures/avg_vel.png")

def plot_avg_positions():
    plt.figure()
    poss = []
    posscomp = []
    for t in range(len(times)):
        poss.append(np.average(gas[t][:gas_indices[t]].position.lengths().value_in(units.parsec)))
        posscomp.append(np.average(gascomp[t].position.lengths().value_in(units.parsec)))
    plt.title("Gas position comparison")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Average gas position (parsec)")
    plt.axvline(SNtime, 0, 1, color="black", linestyle="dashed", label="Start supernova")
    plt.axvline(25, 0, 1, color="gray", linestyle="dashed", label="Start SeBa")
    plt.plot(times, poss, label="With wind")
    plt.plot(times, posscomp, label="Without wind")
    plt.legend()
    plt.savefig("./figures/avg_pos.png")

def plot_energies():
    Us = []
    Uscomp = []
    Ks = []
    Kscomp = []
    for t in range(len(times)):
        print(t)
        #stars
        # Us.append(-bodies[t].potential_energy().value_in(units.m**2 * units.kg * units.s**-2))
        # Uscomp.append(-bodiescomp[t].potential_energy().value_in(units.m**2 * units.kg * units.s**-2))
        # Ks.append(bodies[t].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2))
        # Kscomp.append(bodiescomp[t].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2))
        #gas
        Us.append(-gas[t][:gas_indices[t]].potential_energy().value_in(units.m**2 * units.kg * units.s**-2))
        Uscomp.append(-gascomp[t].potential_energy().value_in(units.m**2 * units.kg * units.s**-2))
        Ks.append(gas[t][:gas_indices[t]].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2))
        Kscomp.append(gascomp[t].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2))
        #both
        # Us.append(-gas[t][:gas_indices[t]].potential_energy().value_in(units.m**2 * units.kg * units.s**-2) - bodies[t].potential_energy().value_in(units.m**2 * units.kg * units.s**-2))
        # Uscomp.append(-gascomp[t].potential_energy().value_in(units.m**2 * units.kg * units.s**-2) - bodiescomp[t].potential_energy().value_in(units.m**2 * units.kg * units.s**-2))
        # Ks.append(gas[t][:gas_indices[t]].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2) + bodies[t].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2))
        # Kscomp.append(gascomp[t].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2) + bodiescomp[t].kinetic_energy().value_in(units.m**2 * units.kg * units.s**-2))
    plt.figure()
    plt.semilogy(times, Us, label="Potential energy")
    plt.semilogy(times, Uscomp, label="Potential energy without winds")
    plt.semilogy(times, Ks, label="Kinetic energy")
    plt.semilogy(times, Kscomp, label="Kinetic energy without winds")
    plt.title("Energy comparison")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Energy (m^2 kg s^-2)")
    plt.axvline(SNtime, 0, 1, color="black", linestyle="dashed", label="Start supernova")
    plt.axvline(25, 0, 1, color="gray", linestyle="dashed", label="Start SeBa")
    plt.legend()
    plt.savefig("./figures/energy.png")


def make_animation(frame):
    print(frame)
    ax[0].cla()
    ax[0].set_title("Time: {:.2f} Myr".format(times[frame]))
    ax[0].set_xlim(-10, 10)
    ax[0].set_ylim(-10, 10)
    ax[0].set_xlabel("parsec")
    ax[0].set_ylabel("parsec")
    
    ax[0].scatter(gas[frame].x[:gas_indices[frame]+1].value_in(units.parsec), gas[frame].y[:gas_indices[frame]+1].value_in(units.parsec), s=0.5, label="Gas", alpha=0.5)
    ax[0].scatter(gas[frame].x[gas_indices[frame]+1:].value_in(units.parsec), gas[frame].y[gas_indices[frame]+1:].value_in(units.parsec), s=1, label="Supernova gas", color="yellow")
    ax[0].scatter(bodies[frame].x.value_in(units.parsec), bodies[frame].y.value_in(units.parsec), s=3, label="Stars")
    # if times[frame] >= SNtime:
    #     ax.scatter(bodies[frame][np.argmax(bodies[frame].mass)].x.value_in(units.parsec), bodies[frame][np.argmax(bodies[frame].mass)].y.value_in(units.parsec), s=10, c="red", label="Supernova")
    # ax.legend()
    ax[1].cla()
    ax[1].set_xlim(-10, 10)
    ax[1].set_ylim(-10, 10)
    ax[1].set_xlabel("parsec")
    ax[1].set_ylabel("parsec")

    ax[1].scatter(gascomp[frame].x.value_in(units.parsec), gascomp[frame].y.value_in(units.parsec), s=0.5, label="Gas", alpha=0.5)
    ax[1].scatter(bodiescomp[frame].x.value_in(units.parsec), bodiescomp[frame].y.value_in(units.parsec), s=3, label="Stars")
    



def make_columndensity(frame):
    print(frame)
    ax.cla()
    ax.set_title("Time: {:.2f} Myr".format(times[frame]))
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.set_xlabel("parsec")
    ax.set_ylabel("parsec")

    plt.hist2d(gas[frame].x[:gas_indices[frame]].value_in(units.parsec), gas[frame].y[:gas_indices[frame]].value_in(units.parsec), bins=(100, 100), range=((-10, 10), (-10, 10)), cmap="OrRd")
    ax.scatter(bodies[frame].x.value_in(units.parsec), bodies[frame].y.value_in(units.parsec), s=3, label="Stars")


# fig, ax = plt.subplots(2)
# ax[0].set_box_aspect(1)
# ax[1].set_box_aspect(1)
# anim = animation.FuncAnimation(fig, make_animation, frames=len(times))
# writer = animation.FFMpegWriter(fps=len(times)/20.)
# anim.save("./figures/animation.mp4", writer=writer)

# fig, ax = plt.subplots()
# anim = animation.FuncAnimation(fig, make_columndensity, frames=len(times))
# writer = animation.FFMpegWriter(fps=len(times)/20.)
# anim.save("./figures/columndensity.mp4", writer=writer)

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

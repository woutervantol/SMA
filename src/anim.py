from os import getcwd, chdir
from sph import *
from matplotlib import animation

fig, ax = plt.subplots()

dt_SN = 0.01 | units.Myr

model_time = 0 | units.Myr
t_steps = np.arange(model_time.value_in(units.Myr), t_end.value_in(units.Myr), dt.value_in(units.Myr)) | units.Myr

def fix_cwd():
    folder = getcwd().split("/")[-1]
    print(folder)
    if folder == "src":
        chdir("../")
    folder = getcwd().split("/")[-1]
    print(folder)

def makeplot(time):
    ax.cla()
    ax.set_title("{:.2f} Myr".format(time.value_in(units.Myr)))
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)
    ax.set_xlabel("parsec")
    ax.set_ylabel("parsec")
    ax.scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), s=0.5, label="gas")
    ax.scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), s=3, label="stars")
    ax.scatter(bodies[np.argmax(bodies.mass)].x.value_in(units.parsec), bodies[np.argmax(bodies.mass)].y.value_in(units.parsec), s=5, c="red")


def calc_time():
    model_time = 0 | units.Myr
    while star_control(bodies, n_stars) != 4:
        model_time = np.round((model_time+dt).value_in(units.Myr), 4) | units.Myr
        yield model_time, dt
    endtime = model_time + (1|units.Myr)
    while model_time < endtime:
        print(star_control(bodies, n_stars))
        model_time = np.round((model_time+dt_SN).value_in(units.Myr), 6) | units.Myr
        yield model_time, dt_SN



def update(times, gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t_end, dt, dt_bridge, n_stars):
    time = times[0]
    print(time)
    timestep = times[1]
    if timestep == dt:
        gravity, hydro, gravhydro, evolution, wind, bodies, gas = simulate(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, time)
        makeplot(time)
    else:
        hydro.parameters.timestep = dt_SN / 2
        gravhydro.timestep = dt_SN * 2
        gravity, hydro, gravhydro, evolution, wind, bodies, gas = simulate(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, time)
        makeplot(time)

# gravity, hydro, gravhydro, evolution, wind, bodies, gas = simulate(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, 6.400999999999994 | units.Myr)
anim = animation.FuncAnimation(fig, update, frames=calc_time, fargs=(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t_end, dt, dt_bridge, n_stars), save_count=1000)
plt.legend()
writer = animation.FFMpegWriter(fps=len(t_steps)/6.)
fix_cwd()
anim.save("./figures/animation.mp4", writer=writer)
# plt.show()
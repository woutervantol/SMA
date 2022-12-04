from sph import *
from matplotlib import animation

fig, ax = plt.subplots()



model_time = 0 | units.Myr
t_steps = np.arange(model_time.value_in(units.Myr), t_end.value_in(units.Myr), dt.value_in(units.Myr)) | units.Myr

# for i, t in enumerate(tqdm(t_steps)):
def makeplot(time):
    ax.cla()
    ax.set_title("{:.1f} Myr".format(time.value_in(units.Myr)))
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)
    ax.set_xlabel("parsec")
    ax.set_ylabel("parsec")
    ax.scatter(gas.x.value_in(units.parsec), gas.y.value_in(units.parsec), s=0.5, label="gas")
    ax.scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), s=3, label="stars")
    ax.scatter(bodies[np.argmax(bodies.mass)].x.value_in(units.parsec), bodies[np.argmax(bodies.mass)].y.value_in(units.parsec), s=5, c="red")
    
    
    

def update(time, gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t_end, dt, dt_bridge, n_stars):
    if star_control(bodies, n_stars) != 14:
        gravity, hydro, gravhydro, evolution, wind, bodies, gas = simulate(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, time)
        print(time, "of", t_end)
        makeplot(time)
    else:
        pass

anim = animation.FuncAnimation(fig, update, frames=t_steps, fargs=(gravity, hydro, gravhydro, evolution, wind, channel, bodies, gas, t_end, dt, dt_bridge, n_stars))
plt.legend()
writer = animation.FFMpegWriter(fps=len(t_steps)/6.)
anim.save("./animation.mp4", writer=writer)
# plt.show()
import numpy as np
from amuse.units import (units, constants)
from matplotlib import pyplot as plt
from tqdm import tqdm



n_stars = 1000


alpha_IMF = -2.35
from amuse.ext.salpeter import new_salpeter_mass_distribution
m_stars = new_salpeter_mass_distribution(n_stars, 0.1|units.MSun, 100|units.MSun, alpha_IMF)
total_mass = np.sum(m_stars)

r_cluster = 1.0 | units.parsec

#converter is nodig omdat het anders dimensieloos is, nu kunnen we initial conditions in SI ingeven
from amuse.units import nbody_system
converter=nbody_system.nbody_to_si(m_stars.sum(),r_cluster)
from amuse.ic.plummer import new_plummer_sphere
bodies=new_plummer_sphere(n_stars, convert_nbody=converter)
bodies.scale_to_standard(converter)

# print(dir(bodies))


from amuse.community.ph4.interface import ph4
gravity = ph4(converter)
gravity.particles.add_particles(bodies)
print("gravitational energy: ", gravity.get_total_energy())
print("kinetic energy: ", bodies.kinetic_energy)

channel = gravity.particles.new_channel_to(bodies)

times = np.arange(0, 9, 1) | units.Myr
figs, axes = plt.subplots(3, 3, sharex=True, sharey=True)
axes = axes.flatten()
import seaborn as sns
for i, t in enumerate(tqdm(times)):
    gravity.evolve_model(t)
    channel.copy()
    # sns.kdeplot(x=bodies.x.value_in(units.parsec), y=bodies.y.value_in(units.parsec), ax=axes[i])
    axes[i].scatter(bodies.x.value_in(units.parsec), bodies.y.value_in(units.parsec), s=1, c=np.log(m_stars.value_in(units.MSun)))
    # plt.axis("equal")
plt.show()
gravity.stop()

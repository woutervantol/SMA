import numpy as np
import matplotlib.pyplot as plt

positions = np.load("./positions.npy", allow_pickle=True)
velocities = np.load("./velocities.npy", allow_pickle=True)
times = np.load("./times.npy", allow_pickle=True)

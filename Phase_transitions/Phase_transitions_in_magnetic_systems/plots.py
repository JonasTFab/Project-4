import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("MC_cycles.txt", skiprows=1)
sorted_data = np.transpose(data)
print(sorted_data)
T = sorted_data[0]
num_cycles = sorted_data[1]
ave_energy = sorted_data[2]
ave_mag = sorted_data[3]
ave_energy_squared = sorted_data[4]
ave_mag_squared = sorted_data[5]

plt.plot(num_cycles,ave_energy)
plt.show()
plt.plot(num_cycles,ave_energy_squared)
plt.show()
plt.plot(num_cycles,ave_mag)
plt.show()
plt.plot(num_cycles,ave_mag_squared)
plt.show()

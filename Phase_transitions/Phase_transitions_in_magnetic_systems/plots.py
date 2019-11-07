import matplotlib.pyplot as plt
import numpy as np

def c_plots(random_file,ordered_file):
    #split data random matrix in files into seperate types
    data_random_matrix = np.loadtxt(random_file, skiprows=1)
    data_ordered_matrix = np.loadtxt(ordered_file, skiprows=1)
    sorted_data_random = np.transpose(data_random_matrix)
    sorted_data_ordered = np.transpose(data_ordered_matrix)
    #split data random matrix in files into seperate types
    T = sorted_data_random[0]
    num_cycles = sorted_data_random[1]
    ave_energy_r = sorted_data_random[2]
    ave_mag_r = sorted_data_random[3]
    ave_energy_squared_r = sorted_data_random[4]
    ave_mag_squared_r = sorted_data_random[5]
    accepted_configs_r = sorted_data_random[6]

    #split data ordered matrix in files into seperate types
    ave_energy_o = sorted_data_ordered[2]
    ave_mag_o = sorted_data_ordered[3]
    ave_energy_squared_o = sorted_data_ordered[4]
    ave_mag_squared_o = sorted_data_ordered[5]
    accepted_configs_o = sorted_data_ordered[6]

    plt.figure(1)#plt.subplot(211)
    plt.title("Average energy (T = %.1f)" %((T[0])))
    plt.plot(num_cycles, ave_energy_o,label = "Ordered matrix",color = 'g')
    #plt.ylabel("E")
    #plt.grid()
    #plt.legend()
    #plt.subplot(212)
    plt.plot(num_cycles,ave_energy_r, label = "Random matrix",color = 'r')
    plt.xlabel('Monte Carlo Cycles')
    plt.ylabel("E")
    plt.grid()
    plt.legend()
    #plt.show()

    plt.figure(2)#plt.subplot(211)
    plt.title("Average magnetization(T = %.1f)" %((T[0])))
    plt.plot(num_cycles, ave_mag_o,label = "Ordered matrix",color = 'g')
    #plt.ylabel("M")
    #plt.grid()
    #plt.legend()
    #plt.subplot(212)
    plt.plot(num_cycles,ave_mag_r,label = "Random matrix",color = 'r')
    plt.xlabel("Monte Carlo Cycles")
    plt.ylabel("M")
    plt.grid()
    plt.legend()
    #plt.show()

    plt.figure(3)#plt.subplot(211)
    plt.title("Accepted configurations[Ac](T = %.1f)" %((T[0])))
    plt.plot(num_cycles, accepted_configs_o,label = "Ordered matrix",color = 'g')
    #plt.ylabel("Ac")
    #plt.grid()
    #plt.legend()
    #plt.subplot(212)
    plt.plot(num_cycles,accepted_configs_r,label = "Random matrix",color = 'r')
    plt.xlabel("Monte Carlo Cycles")
    plt.ylabel("Ac")
    plt.grid()
    plt.legend()
    plt.show()

def d_plots(filename):
    counted_energies = np.loadtxt(filename)
    unique, counts = np.unique(counted_energies, return_counts=True)
    probabilities = counts/len(counted_energies)
    print(counts)
    print(unique)

    #plt.hist(unique)
    #plt.show()
    plt.plot(unique,probabilities)
    plt.show()

d_plots("4d_counted_energies.txt")

#c_plots("MC_cycles_random_T24.txt","MC_cycles_ordered_T24.txt")
#c_plots("MC_cycles_random_T1.txt","MC_cycles_ordered_T1.txt")

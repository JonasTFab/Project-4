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
    plt.figure()
    counted_energies = np.loadtxt(filename)
    mean = counted_energies[0]
    variance = counted_energies[1]
    sigma = np.sqrt(variance)
    unique, counts = np.unique(counted_energies[4000:-1], return_counts=True)
    probabilities = counts/np.sum(counts)
    x = list(filename.split('_')[3])
    b = x
    if b[0] == "r":
        if b[6] == '1':
            T = 1
        else:
            T = 2.4
        title = ("P(E) for random matrix [T = %.1f]" %T)
    else:
        if b[7] == '1':
            T = 1
        else:
            T = 2.4
        title = ("P(E) for ordered matrix [T = %.1f]" %T)
    plt.title(title)
    plt.bar(unique, height = probabilities, width = 5)
    plt.plot(np.linspace(mean,mean),np.linspace(0,max(probabilities)),"gold",label = "Average energy")
    plt.plot(np.linspace(mean-sigma,mean+sigma,2),np.linspace(max(probabilities)/2,max(probabilities)/2,2),markersize = 20,marker ="|",color = "r",label = "Standard deviation")
    plt.legend()
    plt.ylabel('P')
    plt.xlabel('Energy')
    "Remember to write that the expectation value calculated in main is at the top of this"

#average energy for this system -499.874
"First two elements of the following files are 1. average energy, 2. variance"
#d_plots("4d_counted_energies_random1.txt")
#d_plots("4d_counted_energies_ordered1.txt")
#d_plots("4d_counted_energies_random2.txt")
#d_plots("4d_counted_energies_ordered2.txt")
#plt.show()

def e_plots(filename):
    plt.subplot(211)
    dat = np.transpose(np.loadtxt(filename))
    Temp,avg_e,h_capacity,susceptibility,avg_m = dat[0],dat[1],dat[2],dat[3],dat[4]
    plt.plot(Temp,avg_e,"r")
    plt.title("Average energy")
    plt.grid()
    plt.ylabel("E")
    plt.show()
    plt.subplot(212)
    plt.plot(Temp,avg_m,"g")
    plt.label("Average magnetization")
    plt.xlabel("T")
    plt.ylabel("M")
    plt.legend()
    plt.grid()
    plt.show()

    plt.subplot(211)
    plt.plot(Temp,h_capacity, "b",label = ("Heat capacity"))
    plt.xlabel("T")
    plt.ylabel("C_v")
    plt.grid()
    plt.legend()
    plt.subplot(212)
    plt.plot(Temp,susceptibility,"black",label = ("Susceptibility"))
    plt.xlabel("T")
    plt.ylabel("\u03A7")
    plt.grid()
    plt.legend()
    plt.show()

e_plots("4e_L40.txt")

#c_plots("MC_cycles_random_T24.txt","MC_cycles_ordered_T24.txt")
#c_plots("MC_cycles_random_T1.txt","MC_cycles_ordered_T1.txt")

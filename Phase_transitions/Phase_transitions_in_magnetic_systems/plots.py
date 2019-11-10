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

#c_plots("MC_cycles_random_T24.txt","MC_cycles_ordered_T24.txt")
#c_plots("MC_cycles_random_T1.txt","MC_cycles_ordered_T1.txt")

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

def e_plots(filename,L,color):
    plt.figure(1)
    ax1 = plt.subplot(211)
    dat = np.transpose(np.loadtxt(filename))
    Temp,avg_e,h_capacity,avg_m,susceptibility = dat[0],dat[1],dat[2],dat[3],dat[4]
    plt.plot(Temp,avg_e,color)
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0+box.height*0.04, box.width, box.height])
    plt.grid()
    plt.legend()
    plt.ylabel("E")
    ax2 = plt.subplot(212)
    line = ax2.plot(Temp,avg_m,color,label= "L=%s"%L)#,label =("Average magnetization"))
    plt.xlabel("T")
    plt.ylabel("M")
    plt.grid()
    #plt.show()
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0+box.height*0.04, box.width, box.height])

    # Put a legend to the right of the current axis
    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),ncol = 5,fancybox=True)
    #plt.show()

    plt.figure(2)
    plt.subplot(211)
    plt.plot(Temp,h_capacity, "b")#,label = ("Heat capacity"))
    plt.xlabel("T")
    plt.ylabel("C_v")
    plt.grid()
    plt.legend()
    plt.subplot(212)
    plt.plot(Temp,susceptibility,"black")#,label = ("Susceptibility"))
    plt.xlabel("T")
    plt.ylabel("\u03A7")
    plt.grid()
    plt.legend()



e_plots("4e_L40.txt","40","r")
e_plots("4e_L60.txt","60","b")
e_plots("4e_L80.txt","80","g")
e_plots("4e_L100.txt","100","black")
plt.show()

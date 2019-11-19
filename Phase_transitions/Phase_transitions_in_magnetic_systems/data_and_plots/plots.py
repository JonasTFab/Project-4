import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def c_plots(random_file,ordered_file):
    num_particles = 20*20
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
    plt.plot(num_cycles, ave_energy_o/num_particles,label = "Ordered matrix",color='g',linewidth=1)
    #plt.ylabel("E")
    #plt.grid()
    #plt.legend()
    #plt.subplot(212)
    plt.plot(num_cycles,ave_energy_r/num_particles, label = "Random matrix",color = 'r',linewidth=1)
    plt.xlabel('Monte Carlo Cycles')
    plt.ylabel("E")
    plt.grid()
    plt.legend(fontsize=10)
    #plt.show()

    plt.figure(2)#plt.subplot(211)
    plt.title("Average magnetization(T = %.1f)" %((T[0])))
    plt.plot(num_cycles, ave_mag_o/num_particles,label = "Ordered matrix",color = 'g',linewidth=1)
    #plt.ylabel("M")
    #plt.grid()
    #plt.legend()
    #plt.subplot(212)
    plt.plot(num_cycles,ave_mag_r/num_particles,label = "Random matrix",color = 'r',linewidth=1)
    plt.xlabel("Monte Carlo Cycles")
    plt.ylabel("M")
    plt.grid()
    plt.legend(fontsize=10)
    #plt.show()

    plt.figure(3)#plt.subplot(211)
    plt.title("Accepted configurations[Ac](T = %.1f)" %((T[0])))
    plt.plot(num_cycles, accepted_configs_o,label = "Ordered matrix",color = 'g',linewidth=1)
    #plt.ylabel("Ac")
    #plt.grid()
    #plt.legend()
    #plt.subplot(212)
    plt.plot(num_cycles,accepted_configs_r,label = "Random matrix",color = 'r',linewidth=1)
    plt.xlabel("Monte Carlo Cycles")
    plt.ylabel("Accepts")
    plt.grid()
    plt.legend(fontsize=10)
    plt.show()

c_plots("MC_cycles_random_T24.txt","MC_cycles_ordered_T24.txt")
c_plots("MC_cycles_random_T1.txt","MC_cycles_ordered_T1.txt")

def d_plots(filename,color):
    plt.figure(6)
    counted_energies = np.loadtxt(filename)
    mean = counted_energies[0]
    variance = counted_energies[1]
    sigma = np.sqrt(variance)
    unique, counts = np.unique(counted_energies[4000:-1], return_counts=True)
    probabilities = counts/np.sum(counts)
    min_idx = np.argmin(abs(mean-sigma - unique))
    max_idx = np.argmin(abs(mean+sigma -unique))
    #print(np.sum(probabilities[min_idx:max_idx]))
    x = list(filename.split('_')[3])
    b = x
    if b[0] == "r":
        label = "Random"
        sub_num = 1
        if b[6] == '1':
            T = 1
        else:
            T = 2.4
        title = ("P(E) for random matrix [T = %.1f]" %T)
    else:
        sub_num = 2
        label = "Odered"
        if b[7] == '1':
            T = 1
        else:
            T = 2.4
        title = ("P(E) for ordered matrix [T = %.1f]" %T)
    plt.title("P(E) for ordered and random lattice [T = %.1f]" %T,fontsize=16)
    plt.subplot(2,1,sub_num)
    plt.bar(unique, height = probabilities, width = 5,color = color,label = label)
    plt.plot(np.linspace(mean,mean),np.linspace(0,max(probabilities)),"gold",label = "Average energy",linewidth = 3)
    plt.plot(np.linspace(mean-sigma,mean+sigma,2),np.linspace(max(probabilities)/2,max(probabilities)/2,2),markersize = 20,marker ="|",color = "black",label = "Standard deviation",linewidth = 3)
    plt.legend(fontsize=10)
    plt.ylabel('P',fontsize=15)
    "Remember to write that the expectation value calculated in main is at the top of this"

#average energy for this system -499.874
"First two elements of the following files are 1. average energy, 2. variance"
"""
d_plots("4d_counted_energies_random1.txt",'b')
d_plots("4d_counted_energies_ordered1.txt",'r')
plt.xlabel('Energy')
plt.show()
d_plots("4d_counted_energies_random2.txt",'b')
d_plots("4d_counted_energies_ordered2.txt",'r')
plt.xlabel('Energy',fontsize=15)
plt.show()
"""

def e_plots(filename,L,color):
    particles = (float(L)*float(L))
    nu = 1
    a = 1

    dat = np.transpose(np.loadtxt(filename))
    Temp,avg_e,h_capacity,avg_m,susceptibility = dat[0],dat[1],dat[2],dat[3],dat[4]
    #smoothed_h_capacity = pd.rolling_mean(h_capacity,30) #take a moving average over 30 measurments
    smoothed_h_capacity = pd.Series(h_capacity).rolling(window=20,center = True).mean()
    smoothed_susceptibility = pd.Series(susceptibility).rolling(window=20, center = True).mean()
    idx = np.argmax(smoothed_h_capacity)
    idx1 = np.argmax(smoothed_susceptibility)
    critical_T = Temp[idx]
    #print(Temp[idx],'hcap',L)
    #print(Temp[idx1],'sus',L)

    plt.figure(1)
    ax1 = plt.subplot(211)
    ax1.plot(Temp,avg_e/particles,color)
    plt.title("Average energy and magnetization",fontsize=15)
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0+box.height*0.04, box.width, box.height])
    plt.grid()
    plt.ylabel("E",fontsize=15)
    ax2 = plt.subplot(212)
    line = ax2.plot(Temp,avg_m/particles,color,label= "L=%s"%L)#,label =("Average magnetization"))
    plt.xlabel("T",fontsize=15)
    plt.ylabel("M",fontsize=15)
    plt.grid()
    #plt.show()
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0+box.height*0.04, box.width, box.height])

    # Put a legend to the right of the current axis
    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),ncol = 5,fancybox=True,fontsize=15)


    plt.figure(2)
    ax3 = plt.subplot(211)
    plt.title("Heat capacity and susceptibility",fontsize=15)
    ax3.plot(Temp,h_capacity/particles,color = color,alpha=0.4)#,label = ("Heat capacity"))
    plt.plot(Temp,smoothed_h_capacity/particles,color = color)
    plt.ylabel("$C_v$",fontsize=15)
    plt.grid()
    ax4 = plt.subplot(212)
    ax4.plot(Temp,susceptibility/particles,color = color,alpha = 0.4)#,label = ("Susceptibility"))
    ax4.plot(Temp,smoothed_susceptibility/particles,color = color,label= "L=%s"%L)
    plt.xlabel("T",fontsize=15)
    plt.ylabel("$\chi$",fontsize=15)#("\u03A7")
    plt.grid()
    box = ax4.get_position()
    ax4.set_position([box.x0, box.y0+box.height*0.04, box.width, box.height])
    # Put a legend to the right of the current axis
    ax4.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),ncol = 5,fancybox=True,fontsize=15)
    return(critical_T)

"""
T_40 = e_plots("4e_L40.txt","40","r")
T_60 = e_plots("4e_L60.txt","60","b")
T_80 = e_plots("4e_L80.txt","80","g")
T_100 = e_plots("4e_L100.txt","100","black")
plt.show()
"""


def lin_reg():
    nu = 1.
    L_size = np.array([40.,60.,80.,100.])
    Tc_vals = np.array([T_40,T_60,T_80,T_100])

    lin_fit = np.polyfit(Tc_vals/L_size,Tc_vals,1)
    a = lin_fit[0]
    print(lin_fit)
    x = np.linspace(-0.01,0.055,100)
    critical_T = a*L_size**(-nu)+Tc_vals
    print(critical_T)
    plt.plot(Tc_vals/L_size,Tc_vals,"o",color = "r")
    plt.plot(x,lin_fit[0]*x +lin_fit[1])
    plt.grid()
    plt.show()
#lin_reg()

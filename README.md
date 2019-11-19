# Phase Transitions in Magnetic Systems

The Ising model is the model to study and the main focus in this project is how a phase transition in a magnetic system develops as the temperature increases. The algorithm is based on the Markov Chain Monte Carlo method which is called the Metropolis algorithm. Our code is tested for relatively small Monte Carlo cycles before going up to 100 000 cycles. Because of large amount of iterations, parallelization are introduced for optimization.

The operations are computed using C++ as the tool of choice with MPI compiler for parallelization. This coise was made as C++ tends to faster than Python. The code including the algorithm is found in the folder Phase_transitions/Phase_transitions_in_magnetic_systems/ as main.cpp. Because of the heavy domputations, data are stored in seperate text files which is found in the subfolder data_and_plots. As the name implies, the data are plotted with Python using MatPlotLib. The .py file that plots any stored values named plots.py is also found in this folder.

Various functions are implemented in the code (creating initial matrices) but the focous is the Ising_model function. Here lies the Metropolis algorithm. The code requires four arguments: 1.(lattice length) 2.(Temperature) 3.(random/ordered) 4.(Cycles)

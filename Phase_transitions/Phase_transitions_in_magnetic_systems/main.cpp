#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cstdlib>
#include <random>
#include <sstream>
#include <string>
#include "mpi.h"
//For debugging:
// compile with: g++ main.cpp -o main1 -larmadillo -llapack -lblas
// execute with ./main1 length temperature ordered/random

//For paralellization
//mpicxx  -o main_mpi.x  main.cpp -std=c++11
//mpiexec -n 2 ./main_mpi.x

double J = 1;
double k_b = 1;
//std::mt19937 generator (time(NULL)); //seed rng with time now
//output files
std::ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}


int spin(std::mt19937 &generator){ //generate random spins up or down with mersenne twister
  std::uniform_real_distribution<double> dis(0.0, 1.0);
  double ran_nr = dis(generator);
  double divide = 0.5;
  //double inverse_period = RAND_MAX;
  //double ran_nr = rand()/inverse_period;
  if (ran_nr < divide){
        return -1;
    }
    else {
        return 1;
      }
} // end of function spin()

arma::Mat<double> spin_system(int L,std::string ordering, std::mt19937 &generator){ //set up the lattice of spins with random spins up or down
  arma::Mat<double> spin_matrix = arma::mat(L, L,arma::fill::ones);
  if (ordering == "random"){
  for (int i = 0; i < L; i++){
    for(int j = 0; j < L; j++){
      spin_matrix(i,j) = spin(generator);
    }
    }
  }
  //std::cout << spin_matrix << std::endl;
  return spin_matrix;
} // end of function spin_system()

arma::Mat<double> ising_model(int L, double T, arma::mat spin_matrix, int MC_cycles, arma::vec save_energies,
            double &ave_energy, double &spec_heat_cap, double &ave_mag, double &susceptibility, std::mt19937 generator){
  std::uniform_real_distribution<double> dis(-1, L); //chose a random spin to flip
  std::uniform_real_distribution<double> r_dis(0.0, 1.0);
  double energy = 0;
  int Magnetization = 0;
  int accepted_configs = 0;           // number of accepted cofigurations
  int N = L*L;                        // total number of spins
  int M = pow(2,N);                   // total number of possible states
  double beta = 1/((double) k_b*T);
  arma::Mat<double> w = arma::vec(17);
  for (int de = -8; de <= 8; de += 4){
    w(de+8) = exp(-de/T);
  }
  //calculate initial energy
  for (int x = 0; x < L; x++){
    for(int y = 0; y < L; y++){
    energy -= spin_matrix(x,y) * (spin_matrix(periodic(x,L,-1),y) + spin_matrix(x,periodic(y,L,-1)));
    Magnetization += spin_matrix(x,y);
    }
  }

  double ave_energy_squared = 0;
  double ave_mag_squared = 0;

  //The Metropolis Algorithm
  for (int i = 0; i < MC_cycles; i++){
    for (int s=0; s < N; s++){
      int random_x = dis(generator);//random i index to flip
      int random_y = dis(generator);//random j index to flip

      int delta_energy = 2*spin_matrix(random_x,random_y) *
      (spin_matrix(periodic(random_x,L,-1),random_y) +
      spin_matrix(periodic(random_x,L,1),random_y) +
      spin_matrix(random_x,periodic(random_y,L,-1)) +
      spin_matrix(random_x,periodic(random_y,L,1)));

      if (r_dis(generator) <= w(delta_energy + 8)) {
        spin_matrix(random_x,random_y) *= (-1);
        energy += delta_energy;
        Magnetization += 2*spin_matrix(random_x,random_y);
        accepted_configs += 1;
      }
    }
    save_energies(i) = energy;
    // Updating the expectation values
    ave_energy += energy;
    ave_energy_squared += energy*energy;
    ave_mag += abs(Magnetization);
    ave_mag_squared += Magnetization*Magnetization;
  }

  // Normalize
  ave_energy /= (double) MC_cycles;
  ave_energy_squared /= (double) MC_cycles;
  ave_mag /= (double) MC_cycles;
  ave_mag_squared /= (double) MC_cycles;

  spec_heat_cap = (ave_energy_squared - ave_energy*ave_energy)/(k_b*T*T);
  susceptibility = (ave_mag_squared - ave_mag*ave_mag)/(k_b*T);
  save_energies(0) = ave_energy;
  save_energies(1) = spec_heat_cap*(k_b*T*T);
  /*
  std::cout << "Average energy:                            " << ave_energy << std::endl;
  std::cout << "Average energy squared:                    " << ave_energy_squared << std::endl;
  std::cout << "Specific heat capacity:                    " << spec_heat_cap << std::endl;
  std::cout << "Average magnetization:                     " << ave_mag << std::endl;
  std::cout << "Average magnetization squared:             " << ave_mag_squared << std::endl;
  std::cout << "Susceptibility:                            " << susceptibility << "\n\n";
  */

  if (L==2){
    // Analytic solution of mean energy, mean energy squared, mean
    // magnetization and mean magnetization squared for the special case L=2
    double Z = 2*exp(-8*beta*J) + 2*exp(8*beta*J) + 12;
    double an_ave_energy = -(16*J/Z) * (exp(8*beta*J) - exp(-8*beta*J));
    double an_ave_energy_squared = (128*J*J/Z) * (exp(8*beta*J) + exp(-8*beta*J));
    double an_ave_mag = (8*exp(8*beta*J) + 16) / Z;
    double an_ave_mag_squared = (32*exp(8*beta*J) + 32) / Z;
    double an_spec_heat_cap = (1/(k_b*T*T)) * (an_ave_energy_squared - an_ave_energy*an_ave_energy);
    double an_susceptibility = (1/(k_b*T)) * (an_ave_mag_squared - an_ave_mag*an_ave_mag);
    /*
    std::cout << "Analytic average energy (L=2):                   " << an_ave_energy << std::endl;
    std::cout << "Analytic average energy squared (L=2):           " << an_ave_energy_squared << std::endl;
    std::cout << "Analytic specific heat capacity (L=2):           " << an_spec_heat_cap << std::endl;
    std::cout << "Analytic average magnetization (L=2):            " << an_ave_mag << std::endl;
    std::cout << "Analytic average magnetization squared (L=2):    " << an_ave_mag_squared << std::endl;
    std::cout << "Analytic susceptibility (L=2):                   " << an_susceptibility << "\n\n";
    */

    // Unit testing: tests if the numerical calculation are more or less
    // equal the analytical solution. The error is set relatively high because
    // we are working we random number distribution which basically never gets
    // the exact same results as the analytical solution unless lattice is very high.
    double p = 10;
    double eps_e = fabs(an_ave_energy/100);
    double eps_m = fabs(an_ave_mag/100);
    double eps_c = p*fabs(an_spec_heat_cap/100);
    double eps_x = p*fabs(an_susceptibility/100);
    double diff_energy = fabs(fabs(ave_energy)-fabs(an_ave_energy));
    double diff_mag = fabs(fabs(ave_mag)-fabs(an_ave_mag));
    double diff_shc = fabs(fabs(spec_heat_cap)-fabs(an_spec_heat_cap));
    double diff_sus = fabs(fabs(susceptibility)-fabs(an_susceptibility));
    if (diff_energy > eps_e){
      std::string msg = "The difference in analytical and numerical average energy is larger than: ";
      std::cout << msg << eps_e << " Difference is: " << diff_energy << std::endl;
      }
    if (diff_mag > eps_m){
      std::string msg = "The difference in analytical and numerical average magnetization is larger than: ";
      std::cout << msg << eps_m << " Difference is: " << diff_mag << std::endl;
      }
    if (diff_shc > eps_c){
      std::string msg = "The difference in analytical and numerical specific heat capacity is larger than: ";
      std::cout << msg << eps_c << " Difference is: " << diff_shc << std::endl;
      }
    if (diff_sus > eps_x){
      std::string msg = "The difference in analytical and numerical susceptibility is larger than: ";
      std::cout << msg << eps_x << " Difference is: " << diff_sus << std::endl;
      }
    if (diff_sus < eps_x && diff_mag < eps_m && diff_shc < eps_c && diff_sus < eps_x){
      std::string msg = "Unit test complete! The algorithm works just fine for L=2, Thumbs up! We may now assume it works for any L.";
      std::cout << msg << std::endl;
    }

    }



  /*
  // Writes to file. Used for storing the computed data for
  // increasing value of Monte Carlo cycles
  ofile << std::setw(15) << std::setprecision(10) << T;
  ofile << std::setw(15) << std::setprecision(10) << MC_cycles;
  ofile << std::setw(15) << std::setprecision(10) << ave_energy;
  ofile << std::setw(15) << std::setprecision(10) << ave_mag;
  ofile << std::setw(15) << std::setprecision(10) << ave_energy_squared;
  ofile << std::setw(15) << std::setprecision(10) << ave_mag_squared;
  ofile << std::setw(15) << std::setprecision(10) << accepted_configs;
  ofile << "\n";*/

  return save_energies;
} // end of function ising_model()





int main(int argc, char* argv[]){
  int N;

  if (argc != 5){
    std::cout << "Bad usage! Enter on command line: 1.(./filename) 2.(Lattice_length) 3.(Temperature) 4.(random/ordered) 5. number of monte carlo cycles";
    exit(1);
  }
  int L = atoi(argv[1]);
  double Temp = atof(argv[2]);
  std::string ordering = argv[3];
  int MC_cycles = atoi(argv[4]);
  std::string fileout;

  double ave_energy;
  double spec_heat_cap;
  double ave_mag;
  double susceptibility;

  if (Temp != 0){
    std::mt19937 generator (time(NULL));   //seed rng with time now
    //define filename of the utput file
    if (ordering=="random"){
    fileout = "MC_cycles_random.txt";
    }
    else {
    fileout = "MC_cycles_ordered.txt";
    }

    /*
    std::string argument = std::to_string(MC_cycles);
    fileout.append(argument);
    fileout.append("_T_");
    std::string argument2 = std::to_string(Temp);
    fileout.append(argument2);
    fileout.append(ordering);
    */


    ofile.open(fileout);
    arma::Mat<double> matrix = spin_system(L,ordering,generator);
    /*
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << std::setw(15)  << "Temp:";
    ofile << std::setw(15)  << "MC_cycles:";
    ofile << std::setw(15)  << "Avg_E:";
    ofile << std::setw(15)  << "Avg_E^2:";
    ofile << std::setw(15)  << "Avg_mag:";
    ofile << std::setw(15)  << "Avg_mag^2:";
    ofile << std::setw(15)  << "accepted conf:\n";
    */

    // Testing our algorithm for increasing size of Monte Carlo cycles.
    // Then store in a text file for plotting using Python
    //for (int num_cycles = 100; num_cycles <= MC_cycles; num_cycles += 100){
    //  arma::Mat<double> save_energies = arma::vec(MC_cycles);
    //  ising_model(L,Temp,matrix,MC_cycles,save_energies,
    //          ave_energy, spec_heat_cap, ave_mag, susceptibility, generator);
    //}


    //ofile.close();
    arma::Mat<double> save_energies = arma::vec(MC_cycles); //vector contains all energies and will be used to calculate probabilities
    // Starting the algorithm with a unit test (L=2)
    arma::Mat<double> output_save_energies_test = ising_model(2,Temp,matrix,MC_cycles,save_energies,
                                  ave_energy, spec_heat_cap, ave_mag, susceptibility, generator);
    arma::Mat<double> output_save_energies = ising_model(L,Temp,matrix,MC_cycles,save_energies,
                                  ave_energy, spec_heat_cap, ave_mag, susceptibility, generator);

    //std::string output_file = "4d_counted_energies_";
    //output_file.append(ordering);
    //output_file.append(std::to_string(int(Temp)));
    //output_file.append(".txt");
    //std::cout << output_file << std::endl;
    //ofile.open(output_file);
    //ofile <<  output_save_energies;
    //ofile.close();



  }



  else{
    int num_procs, proc_rank;
    double T_max = 2.4;
    double T_min = 2.0;
    int steps = 200;
    double delta_T = (T_max-T_min)/(steps-1);

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank (MPI_COMM_WORLD, &proc_rank);

    MPI_Bcast(&MC_cycles,1,MPI_INT,0, MPI_COMM_WORLD);
    MPI_Bcast(&T_min,1,MPI_INT,0, MPI_COMM_WORLD);
    MPI_Bcast(&T_max,1,MPI_INT,0, MPI_COMM_WORLD);
    MPI_Bcast(&L,1,MPI_INT,0, MPI_COMM_WORLD);
    MPI_Bcast(&delta_T,1,MPI_DOUBLE,0, MPI_COMM_WORLD);


    arma::Col <double> Temp_vec = arma::vec(steps);
    for (int i=0; i<steps; i++){
      Temp_vec(i) = T_min + i*delta_T;
    }

    // MPI initializations
    std::random_device rd;
    std::mt19937 generator (time(NULL) << proc_rank); //seed for different ranks
    int runs = steps/num_procs;
    std::string test_file = "test_file.txt";
    ofile.open(test_file);

    int n_array = 4;

    double *quantity_vec = new double [n_array];
    double *new_quantity_vec = new double [n_array];
    //double rank_T = Temp_vec(proc_rank+num_procs*i);
    double rank_T;

    int old_work_done = 0;
    double time_start = MPI_Wtime();

    if (proc_rank==0){        // unit test: testing if the algorithm works for L=2
      int test_L = 2;
      double test_temp = 2;
      arma::Mat<double> matrix = spin_system(test_L,ordering,generator);
      arma::Mat<double> output_save_energies = ising_model(test_L,test_temp,matrix,MC_cycles,arma::vec(MC_cycles),
                                    ave_energy, spec_heat_cap, ave_mag, susceptibility, generator);
    }

    for (double i=T_min; i<= T_max; i += delta_T){
      arma::Mat<double> matrix = spin_system(L,ordering,generator);

      ising_model(L, i, matrix, MC_cycles, arma::vec(MC_cycles),
                  ave_energy, spec_heat_cap, ave_mag, susceptibility, generator);

      quantity_vec[0] = ave_energy;
      quantity_vec[1] = spec_heat_cap;
      quantity_vec[2] = ave_mag;
      quantity_vec[3] = susceptibility;

      MPI_Reduce(&quantity_vec[0], &new_quantity_vec[0], 4,MPI_DOUBLE,
      MPI_SUM, 0, MPI_COMM_WORLD);
      if(proc_rank == 0){
        ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
        ofile << std::setw(15)  << i;
        ofile << std::setw(15)  << new_quantity_vec[0]/num_procs;
        ofile << std::setw(15)  << new_quantity_vec[1]/num_procs;
        ofile << std::setw(15)  << new_quantity_vec[2]/num_procs;
        ofile << std::setw(15)  << new_quantity_vec[3]/num_procs << std::endl;

        int work_done = 100*(i-T_min) / (T_max-T_min);
        if (work_done != old_work_done){
          std::cout << work_done << "% done" << std::endl;
          old_work_done = work_done;
        }
      }
    }
    if (proc_rank == 0){
      int time_taken = MPI_Wtime() - time_start;
      int hours = time_taken/3600;
      int minutes = time_taken/60 - hours*60;
      int seconds = time_taken - minutes*60 - hours*3600;
      std::cout << "Time taken: " << time_taken << " seconds" << std::endl;
      std::cout << hours << " hour(s) " << minutes << " minute(s) " << seconds << " second(s)" << std::endl;
    }

    MPI_Finalize ();

    ofile.open(fileout);
    arma::Mat<double> matrix = spin_system(L,ordering,generator);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << std::setw(15)  << "accepted conf:\n";


  }
} // end of function main()

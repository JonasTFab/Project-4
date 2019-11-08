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
//mpiexec -n 2 ./main_mpi.x 8

double J = 1;
double k_b = 1;//.38064852e-23;
std::mt19937 generator (time(NULL)); //seed rng with time now
//output files
std::ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}


int spin(){ //generate random spins up or down with mersenne twister
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

arma::Mat<double> spin_system(int L,std::string ordering){ //set up the lattice of spins with random spins up or down
  arma::Mat<double> spin_matrix = arma::mat(L, L,arma::fill::ones);
  if (ordering == "random"){
  for (int i = 0; i < L; i++){
    for(int j = 0; j < L; j++){
      spin_matrix(i,j) = spin();
    }
    }
  }
  //std::cout << spin_matrix << std::endl;
  return spin_matrix;
} // end of function spin_system()

arma::Mat<double> ising_model(int L, double T, arma::mat spin_matrix, int MC_cycles, arma::vec save_energies){
  std::uniform_real_distribution<double> dis(-1, L); //chose a random spin to flip
  std::uniform_real_distribution<double> r_dis(0.0, 1.0);
  double energy = 0;
  int Magnetization = 0;
  int accepted_configs = 0; //number of accepted cofigurations
  int N = L*L;
  int M = pow(2,N);
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

  double ave_energy = 0;
  double ave_energy_squared = 0;
  double ave_mag = 0;
  double ave_mag_squared = 0;

  //The Metropolis Algorithm
  //for (int i = 0; i < MC_cycles; i++){      // without MPI
  for (int i = 0; i < MC_cycles; i++){      // with MPI
    for (int s=0; s < N; s++){
      int random_x = dis(generator);//random i index to flip
      int random_y = dis(generator);//random j index to flip

      int delta_energy = 2*spin_matrix(random_x,random_y) *
      (spin_matrix(periodic(random_x,L,-1),random_y) +
      spin_matrix(periodic(random_x,L,1),random_y) +
      spin_matrix(random_x,periodic(random_y,L,-1)) +
      spin_matrix(random_x,periodic(random_y,L,1)));

      //std::cout << delta_energy << std::endl;
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


  double spec_heat_cap = (ave_energy_squared - ave_energy*ave_energy)/(k_b*T*T);
  double susceptibility = (ave_mag_squared - ave_mag*ave_mag)/(k_b*T);
  save_energies(0) = ave_energy;
  save_energies(1) = spec_heat_cap*(k_b*T*T);
  /*
  std::cout << "Average energy:                            " << ave_energy << std::endl;
  std::cout << "Average energy squared:                    " << ave_energy_squared << std::endl;
  std::cout << "Specific heat capacity:                    " << spec_heat_cap << std::endl;
  std::cout << "Average magnetization:                     " << ave_mag << std::endl;
  std::cout << "Average magnetization squared:             " << ave_mag_squared << std::endl;
  std::cout << "Susceptibility:                            " << susceptibility << "\n\n";*/

  // Analytic solution of mean energy, mean
  // energy squared, mean magnetization and
  // mean magnetization squared
  double Z = 2*exp(-8*beta*J) + 2*exp(8*beta*J) + 12;
  double an_ave_energy = -(16*J/Z) * (exp(8*beta*J) - exp(-8*beta*J));
  double an_ave_energy_squared = (128*J*J/Z) * (exp(8*beta*J) + exp(-8*beta*J));
  double an_ave_mag = (8*exp(8*beta*J) + 16) / Z;
  double an_ave_mag_squared = (32*exp(8*beta*J) + 32) / Z;
  double an_spec_heat_cap = (1/(k_b*T*T)) * (an_ave_energy_squared - an_ave_energy*an_ave_energy);
  double an_susceptibility = (1/(k_b*T)) * (an_ave_mag_squared - an_ave_mag*an_ave_mag);
  /*
  std::cout << "Analytic average energy:                   " << an_ave_energy << std::endl;
  std::cout << "Analytic average energy squared:           " << an_ave_energy_squared << std::endl;
  std::cout << "Analytic specific heat capacity:           " << an_spec_heat_cap << std::endl;
  std::cout << "Analytic average magnetization:            " << an_ave_mag << std::endl;
  std::cout << "Analytic average magnetization squared:    " << an_ave_mag_squared << std::endl;
  std::cout << "Analytic susceptibility:                   " << an_susceptibility << "\n\n";*/


  /*
  ofile << std::setw(15) << std::setprecision(10) << T;
  ofile << std::setw(15) << std::setprecision(10) << MC_cycles;
  ofile << std::setw(15) << std::setprecision(10) << ave_energy;
  ofile << std::setw(15) << std::setprecision(10) << ave_mag;
  ofile << std::setw(15) << std::setprecision(10) << ave_energy_squared;
  ofile << std::setw(15) << std::setprecision(10) << ave_mag_squared;
  ofile << std::setw(15) << std::setprecision(10) << accepted_configs;
  ofile << "\n";*/
  //std::cout << energy << std::endl;
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

  if (Temp != 0){
    //std::mt19937 generator (time(NULL));   //seed rng with time now
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
    arma::Mat<double> matrix = spin_system(L,ordering);
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << std::setw(15)  << "Temp:";
    ofile << std::setw(15)  << "MC_cycles:";
    ofile << std::setw(15)  << "Avg_E:";
    ofile << std::setw(15)  << "Avg_E^2:";
    ofile << std::setw(15)  << "Avg_mag:";
    ofile << std::setw(15)  << "Avg_mag^2:";
    ofile << std::setw(15)  << "accepted conf:\n";

    //for (int num_cycles = 100; num_cycles <= MC_cycles; num_cycles += 100){
    //  ising_model(L,Temp,matrix,num_cycles,arma::vec(num_cycles));
    //}

    /*
    ofile.close();
    arma::Mat<double> save_energies = arma::vec(MC_cycles); //vector will contain all energies and will be used to calculate probabilities
    arma::Mat<double> output_save_energies = ising_model(L,Temp,matrix,MC_cycles,save_energies);
    std::string output_file = "4d_counted_energies_";
    output_file.append(ordering);
    output_file.append(std::to_string(int(Temp)));
    output_file.append(".txt");
    std::cout << output_file << std::endl;
    ofile.open(output_file);
    ofile <<  output_save_energies;
    ofile.close();
    */

  }



  else{
    int num_procs, proc_rank;
    // MPI initializations
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank (MPI_COMM_WORLD, &proc_rank);
    //std::cout << "Rank " << proc_rank << " out of " << num_procs << std::endl;
    std::mt19937 generator (time(NULL) << proc_rank); //seed for different ranks

    int local_cycles = MC_cycles/num_procs;

    double* ave_energy;
    double* ave_energy_squared;
    double* ave_mag;
    double* ave_mag_squared;

    arma::Mat<double> matrix = spin_system(L,ordering);
    double Temp = 1;
    ising_model(L, Temp, matrix, MC_cycles/num_procs, arma::vec(MC_cycles/num_procs),
                &ave_energy, &ave_energy_squared, &ave_mag, &ave_mag_squared);

    std::cout << ave_energy << std::endl;

    double sum_ave_energy = 0;


    MPI_Finalize ();

    // Finding expectation values as a function of T
    double T_max = 2.3;
    double T_min = 2.0;
    int steps = 8;
    double delta_T = (T_max-T_min)/steps;
    //double delta_T = 0.05;
    //int steps = (T_max-T_min)/delta_T + 2;



    arma::Col <double> Temp_vec = arma::vec(steps);
    for (int i=0; i<steps; i++){
      Temp_vec(i) = T_min + i*delta_T;
    }
  }

} // end of function main()

















/* JUNK

arma::Mat<double> new_spin_matrix = spin_matrix;
new_spin_matrix(random_x,random_y) *= (-1);//new lattice with one randomly flipped spin
  new_energy = new_spin_matrix(random_x,random_y)*
  (new_spin_matrix(periodic(random_x,L,-1),random_y) +
  new_spin_matrix(periodic(random_x,L,1),random_y) +
  new_spin_matrix(random_x,periodic(random_y,L,-1)) +
  new_spin_matrix(random_x,periodic(random_y,L,1)));


  // Sjekk Coding energy differences på side 16 i statphys
  for (int i=0; i<=M; i++){
    //w(i) = -8*J + 16*i/M;
  }

  //std::cout << energy << std::endl;
  if (delta_energy <= 0){
    spin_matrix = new_spin_matrix;
    energy = new_energy;
  }
  else{
    w(i) = exp(-beta*delta_energy); //this can be precalculated, how?
    double r = r_dis(generator);
    if (r <= w(i)){
      spin_matrix = new_spin_matrix;
      energy = new_energy;
    }
  }


  else{
    w(i) = exp(-beta*delta_energy); //this can be precalculated, how?
    double r = r_dis(generator);
    if (r <= w(i)){
      spin_matrix(random_x,random_y) *= (-1);
     }
std::cout << energy << std::endl;

*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cstdlib>
#include <random>
#include <sstream>
#include <string>
//For debugging:
// compile with: g++ main.cpp -o main1 -larmadillo -llapack -lblas
// execute with ./main1 length temperature ordered/random

//For paralellization
//mpicxx  -o main_mpi.x  main.cpp -std=c++11
//mpiexec -n 2 ./main_mpi.x 8

double J = 1;
double k_b = 1;//.38064852e-23;
std::mt19937 generator (time(NULL)); //seed rng with time now
//output file
std::ofstream ofile;


// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}
//inline function for computing heat capacity at constant volume
//inline double C_V(double T,mean_E){
//return 1/(k_b*T*T)*(mean_E*mean_E)
//}

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
  std::cout << spin_matrix << std::endl;
  return spin_matrix;
} // end of function spin_system()

int ising_model(int L, double T, arma::mat spin_matrix, int MC_cycles){
  std::uniform_real_distribution<double> dis(-1, L); //chose a random spin to flip
  std::uniform_real_distribution<double> r_dis(0.0, 1.0);
  double energy = 0;
  int Magnetization = 0;
  int N = L*L;
  int M = pow(2,N);
  double beta = 1/(k_b*T);
  arma::Mat<double> w = arma::vec(17);
  for (int de = -8; de <= 8; de += 4){
    w(de+8) = exp(-de/T);
  }
  //calculate initial energy
  for (int x = 0; x < L; x++){
    for(int y = 0; y < L; y++){
    energy -= spin_matrix(x,y) * (spin_matrix(periodic(x,L,-1),y) + spin_matrix(x,periodic(y,L,-1))); //*(spin_matrix(i-1,j)+spin_matrix(i,j-1)+spin_matrix(i+1,j)+spin_matrix(i,j+1));//J = 1
    Magnetization += spin_matrix(x,y);
    }
  }

  double ave_energy = 0;
  double ave_energy_squared = 0;
  double ave_mag = 0;
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

    //std::cout << delta_energy << std::endl;
    if (r_dis(generator) <= w(delta_energy + 8)) {
      spin_matrix(random_x,random_y) *= (-1);
      energy += delta_energy;
      Magnetization += 2*spin_matrix(random_x,random_y);
       }
      }
     // Updating the expectation values
     ave_energy += energy;
     ave_energy_squared += energy*energy;
     ave_mag += Magnetization;
     ave_mag_squared += Magnetization*Magnetization;
     }

     // Normalize
     ave_energy /= MC_cycles;
     ave_energy_squared /= MC_cycles;
     ave_mag /= MC_cycles;
     ave_mag_squared /= MC_cycles;

     /*
     std::cout << "Average energy:                  " << ave_energy << std::endl;
     std::cout << "Average energy squared:          " << ave_energy_squared << std::endl;
     std::cout << "Average magnetization:           " << ave_mag << std::endl;
     std::cout << "Average magnetization squared:   " << ave_mag_squared << std::endl;
     */

     // Mean absolute value of the magnetization,
     // Spesific heat capacity C_V and Suscecbtibility
     // using analytical solutions
      double mean_mag = arma::accu(abs(spin_matrix))/N;
      double C_v = J / (k_b*T*T);
      double chi = 1 / (k_b*T);
      //std::cout << "Analytic mean magnetization:      " << mean_mag << std::endl;
      //std::cout << "Analytic specific heat capacity:  " << C_v << std::endl;
      //std::cout << "Analytic suscecbtibility:         " << chi << std::endl;
      ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
      ofile << std::setw(15) << std::setprecision(10) << T;
      ofile << std::setw(15) << std::setprecision(10) << MC_cycles;
      ofile << std::setw(15) << std::setprecision(10) << ave_energy;
      ofile << std::setw(15) << std::setprecision(10) << ave_mag;
      ofile << std::setw(15) << std::setprecision(10) << ave_energy_squared;
      ofile << std::setw(15) << std::setprecision(10) << ave_mag_squared;
      ofile << "\n";

  //std::cout << energy << std::endl;
  return 0;
} // end of function ising_model()

int main(int argc, char* argv[]){
  int L = atoi(argv[1]);
  double Temp = atof(argv[2]);
  std::string ordering = argv[3];
  int MC_cycles = atoi(argv[4]);
  if (argc != 5){
    std::cout << "Bad usage! Enter on command line: 1.(./filename) 2.(Lattice_length) 3.(Temperature) 4.(random/ordered) 5. number of monte carlo cycles";
    exit(1);
  }

  //define filename of the utput file
  std::string fileout = "MC_cycles.txt";
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
  ofile << std::setw(15)  << "Avg_mag^2:\n";
  for (int num_cycles = 100; num_cycles <= MC_cycles; num_cycles += 100){
    ising_model(L,Temp,matrix,num_cycles);
  }

  ofile.close();
  }

















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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cstdlib>
#include <random>
//For debugging:
// compile with: g++ main.cpp -o main1 -larmadillo -llapack -lblas
// execute with ./main1 length temperature

//For paralellization
//mpicxx  -o main_mpi.x  main.cpp -std=c++11
//mpiexec -n 2 ./main_mpi.x 8
double J = 1;
double k_b = 1;//.38064852e-23;
std::mt19937 generator (time(NULL)); //seed rng with time now

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
} // end of function ran()

arma::Mat<double> spin_system(int L){ //set up the lattice of spins with random spins up or down
  arma::Mat<double> spin_matrix = arma::mat(L, L);
  for (int i = 0; i < L; i++){
    for(int j = 0; j < L; j++){
      spin_matrix(i,j) = spin();
    }
  }
  std::cout << spin_matrix << "\n";
  return spin_matrix;
}

int ising_model(int L, double T, arma::mat spin_matrix){
  std::uniform_real_distribution<double> dis(-1, L); //chose a random spin to flip
  double energy = 0;
  int N = L*L;
  int M = pow(2,N);
  double beta = 1/(k_b*T);
  arma::Mat<double> w = arma::vec(5);

  for (int x = 0; x < L; x++){
    for(int y = 0; y < L; y++){
    energy -= spin_matrix(x,y)*(spin_matrix(periodic(x,L,-1),y) + spin_matrix(x,periodic(y,L,-1))); //*(spin_matrix(i-1,j)+spin_matrix(i,j-1)+spin_matrix(i+1,j)+spin_matrix(i,j+1));//J = 1
    }
  }
  std::cout << energy << std::endl;
  double mean_energy = energy/N;
  double Z = (2*exp(-8)+exp(8) + 12);
  for (int i = 0; i < 5; i++){
    double new_energy = 0;
    int random_x = dis(generator);//random i index to flip
    int random_y = dis(generator);//random j index to flip
    arma::Mat<double> new_spin_matrix = spin_matrix;
    new_spin_matrix(random_x,random_y) *= (-1);//new lattice with one randomly flipped spin
    for (int x = 0; x < L; x++){
      for(int y = 0; y < L; y++){
      new_energy -= new_spin_matrix(x,y)*(new_spin_matrix(periodic(x,L,-1),y) + new_spin_matrix(x,periodic(y,L,-1))); //*(spin_matrix(i-1,j)+spin_matrix(i,j-1)+spin_matrix(i+1,j)+spin_matrix(i,j+1));//J = 1
      }
    }
    double delta_energy = new_energy-energy;
    if (delta_energy<= 0){
      energy = new_energy;
    }
    else{
      w(i) = exp(-beta*delta_energy);
    }
    std::cout << energy << std::endl;
  }

  //std::cout << energy << std::endl;
  return mean_energy;
}

int main(int argc, char* argv[]){
  int N;

  int L = atoi(argv[1]);
  double Temp = atof(argv[2]);
  if (argc != 3){
    std::cout << "Bad usage! Enter on command line: 1.(./filename) 2.(Lattice_length) 3.(Temperature)";
    exit(1);
  }
  //std::cout << "chose lattice length L. N = LXL: L= ";
  //std::cout << "chose temperature T. N = LXL: L= ";

  arma::Mat<double> matrix = spin_system(L);
  ising_model(L,Temp,matrix);
  }

















/* JUNK

  if(std::floor(N)!=std::ceil(N)){
    std::cout << "Bad Usage: " << N <<
    "Not an integer" << std::endl;
    exit(1);
    }
  else{
    ;
  }
*/

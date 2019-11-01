#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cstdlib>
#include <random>
double J = 1;
std::mt19937 generator (time(NULL)); //seed rng with time now

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

//we want a sample mean
//sample eror
double system(int L){
  arma::Mat<double> spin_matrix = arma::mat(L, L);
  for (int i = 0; i < L; i++){
    for(int j = 0; j < L; j++){
      spin_matrix(i,j) = spin();
    }
  }
  std::cout << spin_matrix << "\n";
  return 0;
}

int ising_model(int L, arma::mat spin_matrix){
  double energy = 0;
  for (int i = 0; i < L; i++){
    for(int j = 0; j < L; j++){
    energy += spin_matrix(i,j)*(spin_matrix(i-1,j)+spin_matrix(i,j-1)+spin_matrix(i+1,j)+spin_matrix(i,j+1));

    }
  }
  energy *= -J;
  return 0;
}

int main(int argc, char* argv[]){
  int N;
  int L;
  std::cout << "chose lattice length L. N = LXL: L= ";
  std::cin >> L;
  system(L);
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

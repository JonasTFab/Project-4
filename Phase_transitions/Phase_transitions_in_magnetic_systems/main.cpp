#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cstdlib>
#include <random>
double J = 1;
std::mt19937 generator (time(NULL));

int spin(){ //generate random spins up or down
  //srand(time(NULL));// seed random number generator with the time now
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


int ising_model(int L, int N){

  double energy = 0;
  for(int j = 1; j < N; j++){

  }
  return 0;
}

int main(int argc, char* argv[]){

  }
  }

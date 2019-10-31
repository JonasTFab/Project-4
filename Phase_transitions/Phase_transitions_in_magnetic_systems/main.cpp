#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cstdlib>
double J = 1;

int spin(){ //generate random spins up or down
    double divide = 0.5;
    double ran_nr = rand()/RAND_MAX;
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

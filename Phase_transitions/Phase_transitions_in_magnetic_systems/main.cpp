# include <iostream>
# include <armadillo>

double ran(){ //Just a function to make the random initializer look a bit better
    double ran_nr;
    double invers_period = RAND_MAX;
    ran_nr = rand()/invers_period;
        return ran_nr;
} // end of function ran()




int ising_model(int L, int dim, double T){
    double E,M_mean,Cv,chi;
    int N = pow(L,dim);
    int M = pow(dim,N);
    double Z;





  return 0;
}


int main(){
    double T = 1;
    ising_model(2,2,T);

    return 0;
}

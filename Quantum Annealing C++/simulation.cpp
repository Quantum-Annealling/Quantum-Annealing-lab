#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <cstring>
using namespace std;

static const double fraction = 1.0 / (RAND_MAX + 1.0);

void printMatrix(vector<vector<double> > matrix){
    for (int i=0; i<matrix.size(); i++){
        for (int j=0; j<matrix[i].size(); j++){
            cout<< matrix[i][j] << "\t"; 
        }  
        cout << endl;
    }
}
void printMatrix(vector<vector<int> > matrix){
    for (int i=0; i<matrix.size(); i++){
        for (int j=0; j<matrix[i].size(); j++){
            cout<< matrix[i][j] << "\t"; 
        }  
        cout << endl;
    }
}
int matrixSum(vector<vector<int> > matrix){
    int sum = 0;
    for (int i=0; i<matrix.size(); i++)
        for (int j=0; j<matrix[i].size(); j++)
           sum += matrix[i][j];
    return sum;
}
void writeVector(vector<double> vect, string file){
    ofstream ofile;
    int n = file.length(); 
    char filec[n+1]; 
    strcpy(filec, file.c_str()); 
    ofile.open(filec);
    for(int i=0;i<vect.size();i++)
        ofile<<vect[i]<<endl;
    ofile.close();
}
int getRandomInt(int min, int max){ 
    //max is exclusive
    return min + static_cast<int>((max - min) * (rand() * fraction));
}
double getUniform(){
    return rand() * fraction;
}
double energyAfterNN(vector<vector<int> > &spins, int i, int k, vector<vector<double> > &J, int N, int M, double beta, double gamma){
    double energy = 0 ;
    if(i == 0 && k != 0){
        energy = spins[i][k]*(J[i][N-1]*spins[N-1][k] + J[i][(i + 1)%N]*spins[(i+1)%N][k])/M;
        energy += 0.5/beta*log(1.0/tanh(beta*gamma/M))*spins[i][k]*(spins[i][k-1] + spins[i][(k+1)%M]);
    }
    if(k == 0 && i != 0){
        energy = spins[i][k]*(J[i][i-1]*spins[i-1][k] + J[i][(i + 1)%N]*spins[(i+1)%N][k])/M;
        energy += 0.5/beta*log(1.0/tanh(beta*gamma/M))*spins[i][k]*(spins[i][N-1] + spins[i][(k+1)%M]);
    }
    if(i == 0 && k == 0){
        energy = spins[i][k]*(J[i][N-1]*spins[N-1][k] + J[i][(i + 1)%N]*spins[(i+1)%N][k])/M;
        energy += 0.5/beta*log(1.0/tanh(beta*gamma/M))*spins[i][k]*(spins[i][N-1] + spins[i][(k+1)%M]);
    }
    else{
        energy = spins[i][k]*(J[i][i-1]*spins[i-1][k] + J[i][(i + 1)%N]*spins[(i+1)%N][k])/M;
        energy += 0.5/beta*log(1.0/tanh(beta*gamma/M))*spins[i][k]*(spins[i][k-1] + spins[i][(k+1)%M]);
    }
    return energy;
}
void flipSpin(vector<vector<int> > &spins, vector<vector<double> > &J, int N, int M, double beta, double gamma, double &magnetization){
    //choose random spin on NxM grid
    int i = getRandomInt(0,N);
    int k = getRandomInt(0,M);
    //calculates the energy difference for the Metropolis algorithm
    double energyAfter  = energyAfterNN(spins, i, k, J, N, M, beta, gamma);
    double deltaEnergy =  2 * energyAfter; 
    //Metropolis algorithm
    if (deltaEnergy < 0){
        spins[i][k] = -1 * spins[i][k];
        //keeping track of magnetization
        magnetization = magnetization + 2 * spins[i][k];
    }
    else{
        if(getUniform() < exp(-1.0*deltaEnergy*beta)){
            spins[i][k]  = spins[i][k]* -1;
            //keeping track of magnetization
            magnetization = magnetization + 2 * spins[i][k];
        }
    }
}
//initialization of grid of spins
vector<vector<int> > spinsInit(int N, int M){
    vector<vector<int> > spins(N, vector<int>(M));
    for(int i=0; i<N; i++)
         for (int j=0; j<M; j++)
                spins[i][j] = getRandomInt(0,2)*2-1;
    return spins;
}
vector<vector<double> > JInit(int N){
    vector<vector<double> > J(N, vector<double>(N));
    for (int i=0; i<N; i++){
        J[i][(i+1)%N] = 1;
        J[(i+1)%N][i] = 1;
    }
    J[N-1][0] = 1;
    J[0][N-1] = 1;
    return J;
} 
double avrMag(vector<double> MagnetizationHistory, int N, int M){
    double integral = 0;
    for(int i = 0; i < MagnetizationHistory.size()/2; i++)
        integral += MagnetizationHistory[ MagnetizationHistory.size()/2 + i - 1];
    return integral/(MagnetizationHistory.size()/2);
}
void simulation(int NT, float snapNT, vector<vector<int> > &spins, vector<vector<double> > &J, int N, int M, double beta, double gamma, double &magnetization, vector<double> &magnetizationHistory){
    for(int i = 0; i < NT; i++){
        flipSpin(spins,J,N,M,beta,gamma,magnetization);
        magnetization = matrixSum(spins);
        if(fmod(i,snapNT) == 0){
            magnetizationHistory.push_back(magnetization/(N*M));
        }
    }
}
void simulationGammaRange(int NT, float snapNT, vector<vector<int> > &spins, vector<vector<double> > &J, int N, int M, double beta, double gamma, int Nsteps, double &magnetization, vector<double> &magnetizationHistory, vector<double> &gammaHistory, vector<double> &magnetizationGammaHistory){
    for(int i = 0; i < Nsteps; i++){
        spins = spinsInit(N,M);
        J = JInit(N);
        magnetization = matrixSum(spins);
        gammaHistory.push_back(gamma/Nsteps*i);
        cout << "Gamma = " << gammaHistory[i] << endl;
        simulation(NT,snapNT,spins,J,N,M,beta,gammaHistory[i],magnetization,magnetizationHistory);  
        magnetizationGammaHistory.push_back(fabs(avrMag(magnetizationHistory,N,M)));
        cout << (float)(i+1)/Nsteps*100 << "% done" << endl;
    } 
}
int main(int nargs, const char* argv[]){
    srand(time(NULL));
    //constants and simulation parameters 
    double T = 0.001; //termodynamic temperature 
    double gamma = 2; //transverse field strenght factor
    double kB = 1.0; //Boltzman constant
    double beta = 1.0/kB/T;  
    int N = 10; //N spinów kwantowych
    int M = 10; //additional dimension of spins
    int NT = 2000000; //liczba kroków czasowych, w których losowany jest jeden spin
    float snapNT = (float)(NT/10000); //a value of magnetization is saved every 1000 steps of the simulation
    int Gsteps = 32; //number of evaluation points of gamma for the |<s>| = f(gamma) plot
    int Tsteps = 10; //number of times temperature is decreased during annealing
   
    //initialization for first simulation
    vector<double> magnetizationHistory;
    double magnetization = 0;
    vector<vector<int> > spins = spinsInit(N,M);
    vector<vector<double> > J = JInit(N);
    vector<double> magnetizationGammaHistory;
    vector<double> gammaHistory;
    simulationGammaRange(NT, snapNT, spins, J, N, M, beta, gamma, Gsteps, magnetization, magnetizationHistory, gammaHistory, magnetizationGammaHistory);
    //simulation(NT,snapNT,spins,J,N,M,beta,gamma,magnetization,magnetizationHistory);
    writeVector(magnetizationHistory,"outputMagnetizationHistory.dat"); 
    writeVector(gammaHistory,"outputGammaHistory.dat");
    writeVector(magnetizationGammaHistory,"outputMagnetizationGammaHistory.dat");
    
    return 0;
}
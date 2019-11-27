#ifndef SUPERNETWORK_H
#define SUPERNETWORK_H
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstring>
using namespace std;
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file
#include "network.h" // Network class definition

class superNetwork {

private:
	int netID; // Network ID number to indicate the type of the system with specifying the number of buses
	int solverChoice; // 1 for GUROBI, 2 for CVXGEN, 3 for GUROBI Centralized
	int setRhoTuning; // parameter to select adaptive rho, fixed rho, and type of adaptive rho
	int numberOfCont; // gets the number of contingency scenarios in the variable numberOfCont
	int numberOfGenerators; // get the number of generators in the system
	int iterCountAPP; // Iteration counter for APP coarse grain decomposition algorithm
	double alphaAPP = 5.0; // APP Parameter/Path-length
	int consLagDim; // Dimension of the vectors of APP Lagrange Multipliers and Power Generation Consensus	
	double finTol; // Initial Guess of the Final tolerance of the APP iteration/Stopping criterion
	int intervalCount; // count of the dispatch interval to which the particular network instance for the coarse grain belongs
	int lastInterval; // Flas to indicate if the network belongs to last interval: 0=not last interval; 1=last interval
	vector< Network* > contNetVector; // Vector of base-case and contingency scenario network objects
	vector<double> singleNetTimeVec;
	vector<double> largestNetTimeVec;
	double virtualNetExecTime;

public:
	superNetwork( int, int, int, int, int, int ); // constructor
	~superNetwork(); // destructor	
	void runSimulation(int, double[], double[], double[], double[], double[], int, GRBEnv*); // runs the distributed SCOPF simulations using ADMM-PMP with CVXGEN custom solver
	int getGenNumber(); // returns the number of Generators in the network
	double *getPowPrev(); // returns what I think about previous dispatch interval generators
	double *getPowNext(); // returns what I think about next door fellow
	double *getPowSelf(); // returns the values of what I think about myself
	double getvirtualNetExecTime();
}; // end class superNetwork

#endif // SUPERNETWORK_H

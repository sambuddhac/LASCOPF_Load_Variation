// Network class definition.
#ifndef NETWORK_H
#define NETWORK_H
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>
// include the definitions of classes, generator, load, transmission line, node
#include "generator.h"
#include "load.h"
#include "transl.h"
#include "node.h"
#include <vector>
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file
using namespace std;

class Network {

public:
	Network( int, int, int, int, int, int, int, int ); // constructor
	~Network(); // destructor
	void setNetworkVariables( int, int ); // sets variables of the network
	void runSimulation(int, double[], double[], int, int, double[], double[], double[], double[], double[], GRBEnv*); // runs the SCOPF simulation
	void runSimAPPGurobiBase(int, double[], double[], int, double[], double[], double[], double[], double[], GRBEnv*); // runs the APP coarse grain Gurobi OPF for base case
	void runSimAPPGurobiCont(int, double[], double[], int, double[], double[], GRBEnv*); // runs the APP coarse grain Gurobi OPF for contingency scenarios
	void runSimulationCentral(int, double[], double[], double[], double[], double[], GRBEnv*); // runs the centralized GUROBI based OPF simulation
	int retContCount(); // returns the number of contingency scenarios
	int indexOfLineOut( int ); // gets the serial number of transmission line outaged in this scenario 
	int getGenNumber(); // returns the number of Generators in the network
	double *getPowSelf(); // returns the values of what I think about my generation values
	double *getPowPrev(); // returns the values of what I think about my previous generation values
	double *getPowNext(); // returns the values of what I think about my next generation values
	double *getPowSelfGUROBI(); // returns the values of what I think about my generation values
	double *getPowNextGUROBI(); // returns the values of what I think about my next generation values
	double *getPowPrevGUROBI(); // returns the values of what I think about my previous generation values
	double returnVirtualExecTime();

private:
	// Define the Network
	int networkID; // ID number of the network instance
	int scenarioIndex; // ID number of the contingency scenario, corresponding to which the network instance is generated, 0 for base case
	int prePostContScen; // Takes values 1 or 0, to indicate if the network instance belongs to the post- or pre-contingency scenario, respectively
	int genNumber, genFields; // Number of Generators & Fields
	int loadNumber, loadFields; // Number of Loads & Fields
	int translNumber, translFields; // Number of Transmission Lines & Fields
	int deviceTermCount; // Number of device terminals
	int nodeNumber; // Number of Nodes	
	double Rho; // ADMM Tuning Parameter
	int intervalID; // dispatch interval to which this base case or contingency scenario network belongs
	int lastFlag; // 0 is not the last interval and 1 if the last interval
	vector <int> outagedLine; // vector of serial numbers of the lines that has been outaged in a particular contingency scenario
	int OutagedLine; // serial number of the line that has been outaged in a particular contingency scenario
	int contingencyCount; // Total number of contingency scenarios
	int solverChoice; // 1 for GUROBI, 2 for CVXGEN, 3 for centralized
	bool Verbose; // Parameter to decide whether to display intermediate results or not
	vector< double > pSelfBeleif; // Belief about the generator MW output of the generators in this coarse grain from the previous APP iteration
	vector< double > pSelfBeleifInner; // Belief about the generator MW output of the generators in this coarse grain from the previous inner APP iteration
	vector< double > pPrevBeleif; // Belief about the generator MW output of the generators in this coarse grain from the previous APP iteration
	vector< double > pNextBeleif; // Belief about the generator MW output of the generators in this coarse grain from the previous APP iteration
	vector< int > connNodeNumList; // List of identifiers of nodes for the 300 bus system
	vector< int > nodeValList; // List of assigned numerical ranking/order/serial numbers to the nodes for the 300 bus system
	int assignedNodeSer=0; // assigned node serial for nodes of 330 bus system, initialized to zero
	double pSelfBuffer[100]; // Buffer for the generator MW output of the generators in this coarse grain from the previous APP iteration
	double pPrevBuffer[100]; // Buffer for the generator MW output of the generators in this coarse grain from the previous APP iteration
	double pNextBuffer[100]; // Buffer for the generator MW output of the generators in this coarse grain from the previous APP iteration
	double pSelfBufferGUROBI[100]; // Buffer for the generator MW output of the generators in this coarse grain from the previous APP iteration
	double pNextBufferGUROBI[100]; // Buffer for the generator MW output of the generators in next coarse grain from the previous APP iteration
	double pPrevBufferGUROBI[100]; // Buffer for the generator MW output of the generators in previous coarse grain from the previous APP iteration
	// Names of output files string variables
	string matrixResultString;
	string devProdString;
	string iterationResultString; 
	string lmpResultString;
	string objectiveResultString;
	string primalResultString;
	string dualResultString;
	vector<double> genSingleTimeVec;
	vector<double> genADMMMaxTimeVec;
	double virtualExecTime;
	double divConvMWPU; // Divisor, which is set to 100 for all other systems, except two bus system, for which it is set to 1
	// Create vectors of Generators, Loads, Transmission lines and Nodes
	vector< Generator > genObject;
	vector< Load > loadObject;
	vector< transmissionLine > translObject;
	vector< Node > nodeObject;

}; // end class Network

#endif // NETWORK_H

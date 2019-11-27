// Generator class definition.
// Member functions defined in generator.cpp
#ifndef GENERATOR_H
#define GENERATOR_H
#include "gensolverFirstBase.h" // definition of Gensolver class for base case scenario
#include "gensolverIntermediateBase.h" // definition of Gensolver class for base case scenario
#include "gensolverLastBase.h" // definition of Gensolver class for base case scenario
#include "gensolverCont.h" // definition of Gensolver class for contingency scenario
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file
// Include definition of Node class 
//#include "node.h"

class Node; // forward declaration

class Generator {

public:
	Generator( int, int, int, int, int, Node *, GensolverFirstBase &, GensolverInterBase &, GensolverLastBase &, GensolverCont &, int, int ); // constructor 
	~Generator(); // destructor
	int getGenID(); // returns the ID of the Generator
	int getGenNodeID(); // returns the ID of the Node to which the Generator is connected
	void gpowerangleMessage(int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double[], double[]); //const; // Real power and angle iterate
	void gpowerangleMessageGUROBI(int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double[], double[], GRBEnv*); //const; // Real power and angle iterate
	void setGenData(); // Function to set values of generator data
	double genPower(); //const; // Function to return the value of the power iterate
	double genPowerPrev(); //const; // Function to return the value of the power iterate
	double genPowerNext(); //const; // Function to return the value of the power iterate
	double calcPtilde(); //const; // Calculates the difference between Power iterate and average power
	double calcPavInit() const; // gets the Ptilde before iterations start from node
	double calcvtilde() const; // Calculates the difference between v and average v
	double objectiveGen(); // Calculates the objective function value after each iteration of ADMM
	double objectiveGenGUROBI(); // Objective from GUROBI ADMM 
	double getu() const; // Gets the value of price of real power from node
	double calcThetatilde(); //const; // calculates the difference between voltage angle and average voltage angle
	double getv(); // Gets the value of price of voltage angle
	double getPgenPrev();
	double getPgenNext();
	double getRMax();
	double getRMin();
	double getPMax();
	double getPMin();
	double getQuadCoeff();
	double getLinCoeff();
	double getConstCoeff();

private:
	int genID; // Generator object id number
	int numberOfGenerators; // Total number of generators belonging to the network in any particular scenario
	int baseContScenario; // ID number of the scenario to which the generator belongs, 0 for base-case and 1 for contingency-scenarios
	int dispatchInterval;
	int flagLast;
	int contCountGen; // Total number of contingency scenarios
	double Pg, Thetag, PgenPrev, PgenNext; // Power and angle iterates
	Node *connNodegPtr; // connection node object
	double v; // Voltage angle constraint Lagrange Multiplier
	GensolverFirstBase genSolverFirstBase; // Generator Solver object for base case scenario
	GensolverInterBase genSolverInterBase; // Generator Solver object for base case scenario
	GensolverLastBase genSolverLastBase; // Generator Solver object for base case scenario
	GensolverCont genSolverCont; // Generator solver object for contingency scenarios
	int scenarioContCount; // Contingency scenario index; 0 for base case
	double objOpt; // Optimal value of objective
}; //end class Generator

#endif // GENERATOR_H
	

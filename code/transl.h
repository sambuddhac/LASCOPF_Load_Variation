// transmissionLine class definition.
// Member functions defined in transl.cpp
#ifndef TRANSL_H
#define TRANSL_H
// Include definition of Node class 
//#include "node.h"

class Node; // forward declaration

class transmissionLine {

public:
	transmissionLine( int, Node *, Node *, double, double, double, int ); // constructor
	~transmissionLine(); // destructor
	int getTranslID(); // returns the ID of the transmission line
	int getTranslNodeID1(); // returns ID number of node at end-1 to which the transmission line is connected
	int getTranslNodeID2(); // returns ID number of node at end-2 to which the transmission line is connected
	void tpowerangleMessage( double, double, double, double, double, double, double, double, double, double, double, double, double ); //const; // Real power and Voltage angle iterate
	void setTranData(); // Function to set transmission line data values
	double translPower1(); //const; // Function to return the value of the power iterate
	double translPower2(); //const; // Function to return the value of the power iterate
	double calcPtilde1(); //const; // Calculates the difference between Power iterate and average power for end 1
	double calcPavInit1() const; // gets the Ptilde before iterations start from node at end 1
	double calcvtilde1() const; // Calculates the difference between v and average v for end 1
	double getu1() const; // Gets the value of price of real power from node for end 1
	double calcThetatilde1(); //const; // calculates the difference between voltage angle and average voltage angle for end 1
	double getv1(); //const; // Gets the value of price of voltage angle for end 1
	double calcPtilde2(); //const; // Calculates the difference between Power iterate and average power for end 2
	double calcPavInit2() const; // gets the Ptilde before iterations start from node at end 2
	double calcvtilde2() const; // Calculates the difference between v and average v for end 2
	double getu2() const; // Gets the value of price of real power from node for end 2
	double calcThetatilde2(); //const; // calculates the difference between voltage angle and average voltage angle for end 2
	double getv2(); //const; // Gets the value of price of voltage angle for end 1
	double getReactance(); // Gets the reactance of the transmission line
	double getFlowLimit(); // Gets the value of power flow line limit
	int getOutageScenario(); // returns scenario in which the line is outaged

private:
	int translID; // transmissionLine object id number
	int contScenTracker; // Count of contingency scenario when this articlar line is outaged; only for centralized GUROBI
	double ptMax; // Maximum MW transfer capability
	double Pt1, Pt2, Thetat1, Thetat2; // Iterates for angle and power
	double resT, reacT; // Resistance and Reactance of the transmission line
	Node *connNodet1Ptr, *connNodet2Ptr; // pointers to the node objects at the two ends of the transmission line
	double v1, v2; // Voltage angle constraint Lagrange Multiplier

}; //end class transmissionLine

#endif // TRANSL_H
	

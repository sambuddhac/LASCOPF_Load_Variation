// Node class definition.
// Member functions defined in node.cpp
#ifndef NODE_H
#define NODE_H

#include <vector>

using namespace std;

class Node {

public:
	Node( int, int ); // constructor
	~Node(); // destructor
	int getNodeID(); // returns ID of the node to the caller
	double npinitMessage( double ); // calculates initial average power
	double npinitMessageCont( double, int ); // calculates initial average power at contingency
	double devpinitMessage() const; // returns the initial average node power imbalance to the connected devices before the iterations begin
	double devpinitMessageCont( int ) const; // returns the initial average node power imbalance to the connected devices before the iterations begin in the contingency cases
	void powerangleMessage( double, double *, double *, int, int = 0 ); // gets the power, angle, angle price, contingency count, and device identity from the devices
	double PavMessageCont( int = 0 ) const; // Calculate average power after present iteration/node at contingency
	double ThetaavMessageCont( int = 0 ) const; // Calculate average angle after present iteration/node at contingency
	double ThetaavPrevMessageCont( int = 0 ) const; // Return average angle from previous iteration
	double uMessage( int = 0 ); // Real power balance lagrange multiplier iterate
	double vavMessage() const; // Function to return the average price of the voltage angle constraint
	void setgConn(); // Function to set number of generators connected
	void settConn( int ); // Function to set number of transmission lines connected
	void setlConn( int ); // Function to set number of loads connected
	void reset(); // resets the P_avg, Theta_avg, v_avg to zero after each iteration
	void redContNodeCount( int ); // stores the contingency scenario number corresponding to an outaged line in a vector

private:
	int nodeID; // Node object id number
	int gConnNumber, tConnNumber, lConnNumber; // number of generators, transmission lines and loads connected to the node
	double *P_avg; // average values of power
	double *Theta_avg; // average values of voltage angle
	double *Theta_prev_avg; // average values of voltage angle from previous iteration
	double *u; // Lagrange multiplier corresponding to power balance
	double v_avg; // average Lagrange multiplier corresponding to voltage angle constraint
	int continCount; // Number of contingency scenarios
	int PDevCount; // Number of devices connected to the particular node object
	double Pinitavg; // initial average power
	double Pinitavgc; // initial average power at contingency
	int *nodeFlag; // node flag to indicate whether u has been calculated in a particular iteration
	int tranIDPrev, loadIDPrev; // ID numbers of previous line or load connected to the particular node (in order to avoid duplication of contingency connections)
	vector< int > outageSerial; // vector consisting of the contingency indices of the outaged lines connected to a particular node
	vector< double > PinitAvgCont; // vector consisting of the average initial power injections corresponding to contingency scenarios

}; //end class Node

#endif // NODE_H
	

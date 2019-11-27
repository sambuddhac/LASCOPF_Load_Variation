// transmissionLine class definition.
// Member functions defined in transl.cpp
#ifndef TRANSL_H
#define TRANSL_H

class Node; // forward declaration

class transmissionLine {

public:
	transmissionLine( int, Node *, Node *, double, double, double ); // constructor
	~transmissionLine(); // destructor
	int getTranslID(); // returns the ID of the transmission line
	int getTranslNodeID1(); // returns ID number of node at end-1 to which the transmission line is connected
	int getTranslNodeID2(); // returns ID number of node at end-2 to which the transmission line is connected
	void tpowerangleMessage( double, double, double, double, double, double, double, double, double, double, double, double, double, int = 0 ); // Real power and Voltage angle iterate
	void setTranData(); // Function to set transmission line data values
	double translPower1(); // Function to return the value of the power iterate
	double translPower2(); // Function to return the value of the power iterate
	double calcPtilde1( int = 0 ); // Calculates the difference between Power iterate and average power for end 1
	double calcPavInit1() const; // gets the Ptilde before iterations start from node at end 1
	double calcvtilde1() const; // Calculates the difference between v and average v for end 1
	double getu1( int = 0 ) const; // Gets the value of price of real power from node for end 1
	double calcThetatilde1( int = 0 ); // calculates the difference between voltage angle and average voltage angle for end 1
	double getv1( int = 0 ); // Gets the value of price of voltage angle for end 1
	double calcPtilde2( int = 0 ); // Calculates the difference between Power iterate and average power for end 2
	double calcPavInit2() const; // gets the Ptilde before iterations start from node at end 2
	double calcvtilde2() const; // Calculates the difference between v and average v for end 2
	double getu2( int = 0 ) const; // Gets the value of price of real power from node for end 2
	double calcThetatilde2( int = 0 ); // calculates the difference between voltage angle and average voltage angle for end 2
	double getv2( int = 0 ); // Gets the value of price of voltage angle for end 1

protected:
	int translID; // transmissionLine object id number
	double ptMax; // Maximum MW transfer capability
	int deviceNature; // Type of device; 1 for Generating device, 0 otherwise
	double Pt1, Pt2, Thetat1, Thetat2; // Iterates for angle and power
	double resT, reacT; // Resistance and Reactance of the transmission line
	Node *connNodet1Ptr, *connNodet2Ptr; // pointers to the node objects at the two ends of the transmission line
	double v1, v2; // Voltage angle constraint Lagrange Multiplier

}; //end class transmissionLine

#endif // TRANSL_H
	

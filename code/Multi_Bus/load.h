// Load class definition.
// Member functions defined in load.cpp
#ifndef LOAD_H
#define LOAD_H

class Node; // forward declaration

class Load {

public:
	Load( int, Node *, double ); // constructor
	~Load(); // destructor
	int getLoadID(); // returns the ID of the Load
	int getLoadNodeID(); // returns the ID number of the node to which the load is connected
	double pinitMessage(); // Calculated initial Pav
	void lpowerangleMessage( double, double, double, double, int = 0 ); // Real power and Voltage angle iterate
	void setLoadData(); // Function to set load data values
	double calcPtilde( int = 0 ); // Calculates the difference between Power iterate and average power
	double calcPavInit() const; // gets the Ptilde before iterations start from node
	double calcvtilde() const; // Calculates the difference between v and average v
	double getu( int = 0 ) const; // Gets the value of price of real power from node
	double calcThetatilde( int = 0 ); // calculates the difference between voltage angle and average voltage angle
	double getv( int = 0 ); // Gets the value of price of voltage angle

protected:
	int loadID; // Load object id number
	int deviceNature; // 1 for Generating device, 0 otherwise
	double Pl; // MW consumption
	double Thetal; // Voltage angle
	Node *connNodelPtr; // connection node object
	double v; // Voltage angle constraint Lagrange Multiplier

}; //end class Load

#endif // LOAD_H
	

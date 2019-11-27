// loadContingency derived class definition
// Member functions defined in loadcont.cpp
#ifndef LOADCONT_H
#define LOADCONT_H

#include "load.h" // Load class definition

class Node; // forward declaration

class loadContingency : public Load {

private:
	int loadCID; // ID number of the contingency scenario corresponding to the load for the particular contingency scenario

public:
	loadContingency( int, Node *, double, int ); // constructor
	int getLoadContCounter() const; // gets the contingency scenario count of the load
	void setLoadContCounter( int ); // set the contingency count
	double pinitMessageCont(); // return the average power injection at a contingency scenario before the start of iterations
	double calcPavInitc( int ) const; // returns Ptilde before iterations start
}; // end class loadContingency

#endif // LOADCONT_H

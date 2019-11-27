// tlineContingency derived class definition
// Member functions defined in tlinecont.cpp
#ifndef TLINECONT_H
#define TLINECONT_H

#include "transl.h" // transmissionLine class definition

class Node; // forward declaration

class tlineContingency : public transmissionLine {

private:
	int translCID; // ID number of the contingency scenario corresponding to the outaged line
	double translContInd; // 1.0 for line to be analyzed for contingency, 0.0 for not analyzing
	int contingencyFlag; // flag when set to 1 indicates the line is outaged, 0 otherwise

public:
	tlineContingency( int, Node *, Node *, double, double, double, double, int ); // constructor
	int getTranslContCounter() const; // gets the contingency scenario count of the line
	void setTranslCont( double, int ); // sets the contingency scenario count and index of the line
	void reduceNodeCount(); // Reduce the node count by one corresponding to the outaged transmission line
	int getContFlag() const; // Return the flag of contingencyFlag
	double calcPavInitc1( int ) const; // return the Ptilde from sending end before iterations start
	double calcPavInitc2( int ) const; // return the Ptilde from receiving end before iterations start
}; // end class tlineContingency

#endif // TLINECONT_H

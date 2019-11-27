// Member functions for class tlineContingency
// tlineContingency class inherits from class transmissionLine
#include <iostream>
#include "node.h"
// include tlineContingency class definition from tlinecont.h
#include "tlinecont.h"

using namespace std;

// default constructor
tlineContingency::tlineContingency( int idOfTransl, Node *nodeConnt1, Node *nodeConnt2, double PowertMax, double Reactance, double Resistance, double indexOfContingency, int countOfContingency )
		   : transmissionLine( idOfTransl, nodeConnt1, nodeConnt2, PowertMax, Reactance, Resistance ) // call the base-class constructor
{
	setTranslCont( indexOfContingency, countOfContingency );

} // end tlineContingency constructor

// set the contingency count and index
void tlineContingency::setTranslCont( double indexOfContingency, int countOfContingency )
{
	translCID = ( countOfContingency <= 0 ? 0 : countOfContingency );
	translContInd = indexOfContingency;
	contingencyFlag = 0; 

} // end of setTranslCont

// get the contingency count
int tlineContingency::getTranslContCounter() const
{
	return translCID; 

} // end function getTranslContCounter

// Reduce the node count by one corresponding to the outaged line
void tlineContingency::reduceNodeCount()
{
	contingencyFlag = 1; // sets the contingencyFlag to indicate the outaged line
	// Passes to the connected nodes the contingency number
	connNodet1Ptr->redContNodeCount( translCID ); 
	connNodet2Ptr->redContNodeCount( translCID ); 

} // end function tlineContingency

// return the value of contingencyFlag
int tlineContingency::getContFlag() const
{
	return contingencyFlag; 
} // end function getContFlag

double tlineContingency::calcPavInitc1( int contin ) const // function calcPavInitc1 begins
{
	return connNodet1Ptr->devpinitMessageCont( contin ); // seeks the initial Ptilde from the node under contingency
} // function calcPavInitc ends

double tlineContingency::calcPavInitc2( int contin ) const // function calcPavInitc2 begins
{
	return connNodet2Ptr->devpinitMessageCont( contin ); // seeks the initial Ptilde from the node under contingency
} // function calcPavInitc ends

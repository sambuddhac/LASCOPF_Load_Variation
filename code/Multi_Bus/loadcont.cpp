// Member functions for class loadContingency
// loadContingency class inherits from class Load
#include <iostream>
#include "node.h"
// include loadContingency class definition from loadcont.h
#include "loadcont.h"

using namespace std;

// default constructor
loadContingency::loadContingency( int idOfLoad, Node *nodeConnl , double Load_P, int countOfContingency )
		   : Load( idOfLoad, nodeConnl, Load_P ) // call the base-class constructor
{
	setLoadContCounter( countOfContingency );

} // end loadContingency constructor

// set the contingency count
void loadContingency::setLoadContCounter( int countOfContingency )
{
	loadCID = ( countOfContingency <= 0 ? 0 : countOfContingency );

} // end of setLoadContCounter

// get the contingency count
int loadContingency::getLoadContCounter() const
{
	return loadCID; 

} // end function getLoadContCounter

double loadContingency::pinitMessageCont() // function pinitMessageCont begins
{
	double pinit = 0.0; // declare and initialize pinit
	pinit = pinit + connNodelPtr->npinitMessageCont( Pl, loadCID ); // passes to node object the power value to calculate the initial Pav
	return pinit; // return pinit
} // function pinitMessageCont ends

double loadContingency::calcPavInitc( int contin ) const // function calcPavInitc begins
{
	return ( Pl - connNodelPtr->devpinitMessageCont( contin ) ); // seeks the initial Ptilde from the node under contingency
} // function calcPavInitc ends


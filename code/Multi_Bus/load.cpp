// Member functions for class Load.
#include <iostream>

// include Load class definition from load.h and node.h
#include "load.h"
#include "node.h"

using namespace std;

Load::Load( int idOfLoad, Node *nodeConnl, double Load_P ) // constructor begins
	: loadID( idOfLoad ),
	  Pl( Load_P ),
	  connNodelPtr( nodeConnl ),
	  Thetal( 0.0 ),
	  deviceNature( 0 )
{
	//cout << "\nInitializing the parameters of the load with ID: " << loadID << endl;
	connNodelPtr->setlConn( idOfLoad ); // increments the load connection variable to node
	setLoadData(); // calls setLoadData member function to set the parameter values

} // constructor ends

Load::~Load() // destructor
{
	//cout << "\nThe load object having ID " << loadID << " have been destroyed.\n";

} // end of destructor

int Load::getLoadID() // function getLoadID begins
{
	return loadID; // returns the ID of the load object
} // end of getLoadID function

int Load::getLoadNodeID() // function getLoadNodeID begins
{
	return connNodelPtr->getNodeID(); // returns the ID number of the node to which the load object is connected
} // end of getLoadNodeID function

void Load::setLoadData() // start setLoadLimits function
{
	
	v = 0.0; // Initialize the Lagrange multiplier corresponding voltage angle constraint to zero
} // end setLoadLimits function

double Load::pinitMessage() //const // function pinitMessage begins
{
	double pinit = 0.0; // declare and initialize pinit
	pinit = pinit + connNodelPtr->npinitMessage( Pl ); // passes to node object the power value to calculate the initial Pav
	return pinit; // return pinit
} // function pinitMessage ends

void Load::lpowerangleMessage( double lRho, double vprevavg, double Aprevavg, double vprev, int contCount ) //const // function lpowerangleMessage begins
{
	Thetal = vprevavg + Aprevavg - vprev;
	//*cout << "\nFrom load " << loadID << " in the scenario " << contCount;
	//*cout << "\nTheta from load solver: " << Thetal << endl;
	connNodelPtr->powerangleMessage( Pl, &v, &Thetal, deviceNature, contCount ); // passes to node object the corresponding iterates of power, angle and v
} // function lpowerangleMessage ends

double Load::calcPtilde( int contin ) //const // function calcPtilde begins
{
	double P_avg = connNodelPtr->PavMessageCont( contin ); // Gets average power from the corresponding node object
	double Ptilde = Pl - P_avg; // calculates the difference between power iterate and average
	return Ptilde; // returns the difference
} // function calcPtilde ends

double Load::calcPavInit() const // function calcPavInit begins
{
	return ( Pl - connNodelPtr->devpinitMessage() ); // seeks the initial Ptilde from the node
} // function calcPavInit ends

double Load::getu( int contin ) const // function getu begins
{
	double u = connNodelPtr->uMessage( contin ); // gets the value of the price corresponding to power balance from node
	//cout << "u: " << u << endl;
	return u; // returns the price
} // function getu ends

double Load::calcThetatilde( int contin ) //const // function calcThetatilde begins
{
	double Theta_avg = connNodelPtr->ThetaavMessageCont( contin ); // get the average voltage angle at the particular node
	double Theta_tilde = Thetal - Theta_avg; // claculate the deviation between the voltage angle of the device and the average
	//*cout << "\nLoad ID # " << loadID << " contingency scenario# " << contin << " angle = " << Thetal << " , average angle = " << Theta_avg << endl;
	return Theta_tilde; // return the deviation
} // function calcThetatilde ends

double Load::calcvtilde() const // function calcvtilde begins
{
	double v_avg = connNodelPtr->vavMessage(); // get the average of the Lagrange multiplier corresponding to voltage angle balance 
	double v_tilde = v - v_avg; // calculate the deviation of the node Lagrange multiplier to the average
	return v_tilde; // return the deviation
} // function calcvtilde ends

double Load::getv( int contin ) //const // function getv begins
{
	//cout << "v_initial: " << v << endl;
	v = v + calcThetatilde( contin ); // Calculate the value of the Lagrange multiplier corresponding to angle constraint
	//cout << "v_final: " << v << endl;
	return v; // Calculate the value of the Lagrange multiplier corresponding to angle constraint
} // function getv ends	
	

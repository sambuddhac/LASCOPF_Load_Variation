// Member functions for class Node.
#include <iostream>
// include Node class definition from node.h
#include <vector>
#include <algorithm>
#include "node.h"

using namespace std;

Node::Node( int idOfNode, int countOfContin ) // constructor begins
	: nodeID( idOfNode ),
	  continCount( countOfContin )
{
	//cout << "\nInitializing the parameters of the node with ID: " << nodeID << endl;

	// initialize the connected devices to zero for node
	gConnNumber = 0; // number of generators connected to a particular node
	tConnNumber = 0; // number of transmission lines connected to a particular node
	lConnNumber = 0; // number of loads connected to a particular node

	nodeFlag = new int[ ( continCount + 1 ) ]; // flag to indicate if a particular node has been accounted for by any one device connected to it for calculation of u
	for (int i = 0; i < ( continCount + 1 ); ++i ) {
		nodeFlag[ i ] = 0; // Initialize to zero
	}
 
	tranIDPrev = -1; // Value to initiate the ID of transmission line
	loadIDPrev = -1; // Value to initiate the ID of load 
	PDevCount = 0; // initialize number of devices connectedto a node to zero
	P_avg = new double[ ( continCount + 1 ) ]; // Initialize average power to zero
	for (int i = 0; i < ( continCount + 1 ); ++i ) {
		P_avg[ i ] = 0.0; // Initialize the average power to zero
	}
	
	Theta_avg = new double[ ( continCount + 1 ) ]; // initialize average angle to zero
	for (int i = 0; i < ( continCount + 1 ); ++i ) {
		Theta_avg[ i ] = 0.0; // Initialize the voltage angle to zero
	}

	Theta_prev_avg = new double[ ( continCount + 1 ) ]; // initialize average angle to zero
	for (int i = 0; i < ( continCount + 1 ); ++i ) {
		Theta_prev_avg[ i ] = 0.0; // Initialize the voltage angle to zero
	}
	
	u = new double[ ( continCount + 1 ) ]; // initialize power balance price to zero
	for (int i = 0; i < ( continCount + 1 ); ++i ) {
		u[ i ] = 0.0; // initialize power balance price to zero
	} 

	v_avg = 0.0; // initialize average value of voltage angle price to zero
	Pinitavg = 0.0; // initialize initial average power to zero
	Pinitavgc = 0.0; // initialize initial average power at contingency to zero

} // constructor ends

Node::~Node() // destructor
{
	//cout << "\nThe node object having ID " << nodeID << " have been destroyed.\n";

} // end of destructor

int Node::getNodeID() // function getNodeID begins
{
	return nodeID; // returns node ID to the caller
} // end of function getNodeID

void Node::setgConn()
{
	++gConnNumber; // increment the number of generators connected by one whenever a generator is connected to the node

}

void Node::settConn( int tranID )
{
	if ( tranID != tranIDPrev ) {
		++tConnNumber; // increment the number of txr lines connected by one whenever a txr line is connected to the node
		tranIDPrev = tranID; // Update current value of transmission line ID
	}

}

void Node::setlConn( int lID )
{
	if ( loadIDPrev != lID ) {
		++lConnNumber; // increment the number of loads connected by one whenever a load is connected to the node
		loadIDPrev = lID; // Update current value of load ID
	}

}

// function redContNodeCount begins
void Node::redContNodeCount( int translCID ) 
{
	outageSerial.push_back( translCID ); // store the value of the scenario number corresponding to the outaged line in the vector
} // function redContNodeCount ends

double Node::npinitMessage( double Pload ) // function npinitMessage begins
{
	Pinitavg = Pinitavg + Pload / ( gConnNumber + tConnNumber + lConnNumber ); // calculate average power 
	return Pinitavg; // return initial average power
} // function npinitMessage ends

double Node::npinitMessageCont( double Pload, int contin ) // function npinitMessageCont begins
{
	if ( find( outageSerial.begin(), outageSerial.end(), contin ) != outageSerial.end() ) {
		Pinitavgc = Pload / ( gConnNumber + ( tConnNumber - 1 ) + lConnNumber ); // calculate average power under outage
	}
	else {
		Pinitavgc = Pload / ( gConnNumber + tConnNumber + lConnNumber ); // calculate average power
	}
	PinitAvgCont.push_back( Pinitavgc ); // store the value of the average power injection at the particular contingency scenario in the vector 
	return Pinitavgc; // return initial average power
} // function npinitMessageCont ends


double Node::devpinitMessage() const// function devpinitMessage begins
{
	return Pinitavg; // return the initial average node power imbalance to the devices
} // function devpinitMessage ends

double Node::devpinitMessageCont( int contin ) const// function devpinitMessageCont begins
{
	if ( lConnNumber == 0 )
		return Pinitavgc; // if no load is connected to the node, return zero
	else
		return PinitAvgCont[ contin - 1 ]; // return the initial average node power imbalance to the devices under contingency
} // function devpinitMessageCont ends


void Node::powerangleMessage( double Power, double* AngPrice, double* Angle, int devIdentity, int contCountGen ) //const // function powerangleMessage begins
{
	if ( devIdentity == 1 ) {
		P_avg[ 0 ] = P_avg[ 0 ] + Power / ( gConnNumber + tConnNumber + lConnNumber ); // calculate average power under base case
		Theta_avg[ 0 ] = Theta_avg[ 0 ] + *( Angle ) / ( gConnNumber + tConnNumber + lConnNumber ); // calculate average voltage angle under base case
		//*cout << "\nP_avg from node solver: " << P_avg[ 0 ] << "\tTheta_avg from node solver: " << Theta_avg[ 0 ] << endl;
		for ( int i = 1; i <= contCountGen; ++i ) {
			if ( find( outageSerial.begin(), outageSerial.end(), i ) != outageSerial.end() ) {
				P_avg[ i ] = P_avg[ i ] + Power / ( gConnNumber + ( tConnNumber - 1 ) + lConnNumber ); // calculate average power under outage
				Theta_avg[ i ] = Theta_avg[ i ] + *( Angle + i ) / ( gConnNumber + ( tConnNumber - 1 ) + lConnNumber ); // calculate average angle under outage
			}
			else {
				P_avg[ i ] = P_avg[ i ] + Power / ( gConnNumber + tConnNumber + lConnNumber ); // calculate average power
				Theta_avg[ i ] = Theta_avg[ i ] + *( Angle + i ) / ( gConnNumber + tConnNumber + lConnNumber ); // calculate average voltage
			}
			//*cout << "\nP_avg from node solver: " << P_avg[ i ] << "\tTheta_avg from node solver: " << Theta_avg[ i ] << endl;
		}
		//*cout << "\nContingency scenario number from Generator " << contCountGen << endl; 
	}
	else {
		if ( find( outageSerial.begin(), outageSerial.end(), contCountGen ) != outageSerial.end() ) {
			P_avg[ contCountGen ] = P_avg[ contCountGen ] + Power / ( gConnNumber + ( tConnNumber - 1 ) + lConnNumber ); // calculate average power under outage
			Theta_avg[ contCountGen ] = Theta_avg[ contCountGen ] + *( Angle ) / ( gConnNumber + ( tConnNumber - 1 ) + lConnNumber ); // calculate average angle under outage
			//*cout << "\nContingency scenario number from calculation for outaged condition " << contCountGen << endl; 
		}
		else {
			P_avg[ contCountGen ] = P_avg[ contCountGen ] + Power / ( gConnNumber + tConnNumber + lConnNumber ); // calculate average power
			Theta_avg[ contCountGen ] = Theta_avg[ contCountGen ] + *( Angle ) / ( gConnNumber + tConnNumber + lConnNumber ); // calculate average voltage
			//*cout << "\nContingency scenario number from calculation for healthy condition " << contCountGen << endl; 
		}
		//*cout << "\nP_avg from node solver: " << P_avg[ contCountGen ] << "\tTheta_avg from node solver: " << Theta_avg[ contCountGen ] << endl;
	}
	//**v_avg = v_avg + AngPrice / ( gConnNumber + tConnNumber + lConnNumber ); // calculate average voltage angle price
	if ( ( devIdentity == 1 ) || ( ( devIdentity == 0 ) && ( contCountGen == 0 ) ) ) {
		//*cout << "\nContingency scenario number for base-case condition " << contCountGen << endl; 
		PDevCount++; // increment device count by one indicating that a particular device connected to the node has been taken into account
	}
} // function powerangleMessage ends

double Node::PavMessageCont( int cont ) const // function PavMessage begins
{
	if ( PDevCount == ( gConnNumber + tConnNumber + lConnNumber ) )  // if all the devices are taken care of return the average power
		return P_avg[ cont ];
} // function PavMessage ends

double Node::uMessage( int cont ) // function uMessage begins
{
	if ( nodeFlag[ cont ] != 0 ) {
		//cout << nodeFlag << endl;
		return u[ cont ];
	} 
	else {
		if ( PDevCount == ( gConnNumber + tConnNumber + lConnNumber ) ) {
			u[ cont ] += P_avg[ cont ]; 		
			//cout << nodeFlag << endl;
			++nodeFlag[ cont ]; // this node has already been accounted for
			return u[ cont ];  // if all the devices are taken care of calculate and return the power price
		}
	}

} // function uMessage ends

double Node::ThetaavPrevMessageCont( int cont ) const // function ThetaavMessage begins
{
	if ( PDevCount == ( gConnNumber + tConnNumber + lConnNumber ) ) 
		return Theta_prev_avg[ cont ];  // if all the devices are taken care of return the average angle
} // function ThetaavMessage ends

double Node::ThetaavMessageCont( int cont ) const // function ThetaavMessage begins
{
	if ( PDevCount == ( gConnNumber + tConnNumber + lConnNumber ) ) 
		return Theta_avg[ cont ];  // if all the devices are taken care of return the average angle
} // function ThetaavMessage ends

double Node::vavMessage() const // function vavMessage begins
{
	if ( PDevCount == ( gConnNumber + tConnNumber + lConnNumber ) )  
		return v_avg;  // if all the devices are taken care of return the average angle price
} // function vavMessage ends

	
void Node::reset() // function reset begins
{
	PDevCount = 0;
	for (int i = 0; i < ( continCount + 1 ); ++i ) {
		P_avg[ i ] = 0.0; // reset the average power to zero
	}
	v_avg = 0.0;
	for (int i = 0; i < ( continCount + 1 ); ++i ) {
		Theta_prev_avg[ i ] = Theta_avg[ i ]; // Store the previous average node voltage angle from the previous iteration
		Theta_avg[ i ] = 0.0; // reset the voltage angle to zero
	}
	for (int i = 0; i < ( continCount + 1 ); ++i ) {
		nodeFlag[ i ] = 0; // reset to zero
	}
} // function reset ends

		


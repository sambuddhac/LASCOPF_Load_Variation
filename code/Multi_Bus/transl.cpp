// Member functions for class transmissionLine.
#include <iostream>
#include <iomanip>
#include <cmath>
// include transmissionLine class definition from transl.h, Node class definition from node.h
#include "transl.h"
#include "node.h"

using namespace std;

transmissionLine::transmissionLine( int idOfTransl, Node *nodeConnt1, Node *nodeConnt2, double PowertMax, double Reactance, double Resistance ) // constructor begins
	: translID( idOfTransl ),
	  connNodet1Ptr( nodeConnt1 ),
	  connNodet2Ptr( nodeConnt2 ),
	  ptMax( PowertMax ),
	  reacT( Reactance ),
	  resT( Resistance ), 
	  deviceNature( 0 )
{
	//*cout << "\nInitializing the parameters of the transmission line with ID: " << translID << endl;
	connNodet1Ptr->settConn( idOfTransl ); // increments the txr line connection variable to node 1
	connNodet2Ptr->settConn( idOfTransl ); // increments the txr line connection variable to node 2
	setTranData(); // calls setTranData member function to set the parameter values

} // constructor ends

transmissionLine::~transmissionLine() // destructor
{
	//cout << "\nThe transmission line object having ID " << translID << " have been destroyed.\n";

} // end of destructor

int transmissionLine::getTranslID() // function gettranslID begins
{
	return translID; // returns the ID of the generator object
} // end of gettranslID function

int transmissionLine::getTranslNodeID1() // function getGenNodeID begins
{
	return connNodet1Ptr->getNodeID(); // returns the ID number of the node to which the generator object is connected
} // end of getGenNodeID function

int transmissionLine::getTranslNodeID2() // function getGenNodeID begins
{
	return connNodet2Ptr->getNodeID(); // returns the ID number of the node to which the generator object is connected
} // end of getGenNodeID function

void transmissionLine::setTranData() // member function to set parameter values of transmission lines
{
	Thetat1 = 0.0; // Initialize the angle iterate at end-1
	Thetat2 = 0.0; // Initialize the angle iterate at end-2
	Pt1 = 0.0; // Initialize the power iterate at end-1
	Pt2 = 0.0; // Initialize the power iterate at end-2
	v1 = 0.0; // Initialize the Lagrange multiplier corresponding to end-1 voltage angle constraint to zero
	v2 = 0.0; // Initialize the Lagrange multiplier corresponding to end-2 voltage angle constraint to zero

} // end function for setting parameter values

void transmissionLine::tpowerangleMessage( double tRho, double Pprevit1, double Pnetavg1, double uprev1, double vprevavg1, double Aprevavg1, double vprev1,  double Pprevit2, double Pnetavg2, double uprev2, double vprevavg2, double Aprevavg2, double vprev2, int contCount ) //const // function tpowerangleMessage begins
{
if ( ptMax == pow( 10.0, -100.0 ) ) {
	Thetat1 = connNodet1Ptr->ThetaavPrevMessageCont( contCount );
	Thetat2 = connNodet2Ptr->ThetaavPrevMessageCont( contCount );
	Pt2 = 0; // get the transmission line end-2 Power iterate
	Pt1 = 0; // get the transmission line end-2 voltage angle iterate
	//*connNodet1Ptr->powerangleMessage( Pt1, &v1, &Thetat1, deviceNature, contCount ); // passes to node object at end 1 the corresponding iterates of power, angle and v
	//*connNodet2Ptr->powerangleMessage( Pt2, &v2, &Thetat2, deviceNature, contCount ); // passes to node object at end 2 the corresponding iterates of power, angle and v
	//*cout << "\nUnder Contingency case, " << contCount << " resistance, reactance, and capacity of line " << translID << " respectively, are " << resT << "\t" << reacT << "\t" << ptMax << endl;
}
else {
	double end1A = reacT * ( Pprevit1 - Pnetavg1 - uprev1 );
	double end1B = ( vprevavg1 + Aprevavg1 - vprev1 );
	double end2C = reacT * ( Pprevit2 - Pnetavg2 - uprev2 );
	double end2D = ( vprevavg2 + Aprevavg2 - vprev2 );
	double Diff; // difference between the bus voltage angles
	if ( getTranslNodeID1() == 1 ) { // if end-1 is the bus-1, or, slack bus, fix the voltage angle of that end to 0
		Thetat1 = 0.0;
		Diff = ( ( pow( reacT, 2.0 ) ) * end2D + end1A - end2C ) / ( 2.0 + pow( reacT, 2.0 ) );
	}
	else {
		if ( getTranslNodeID2() == 1 ) { // if end-2 is the bus-1, or, slack bus, fix the voltage angle of that end to 0
			Thetat1 = ( ( pow( reacT, 2.0 ) ) * end1B - end1A + end2C ) / ( 2.0 + pow( reacT, 2.0 ) );
			Diff = -Thetat1; 
		}
		else { // if none of the ends is the slack bus, consider both the voltage angles as decision variables and calculate them
			Thetat1 = ( ( 2.0 + pow( reacT, 2.0 ) ) * end1B  - end1A + end2C + ( 2.0 * end2D ) ) / ( 4.0 + pow( reacT, 2.0 ) ); // Thetat1 iterate
			Diff = ( ( 2.0 * end1A ) - ( pow( reacT, 2.0 ) ) * end1B - ( 2.0 * end2C ) + ( pow( reacT, 2.0 ) ) * end2D ) / ( 4.0 + pow( reacT, 2.0 ) ); // Magnitude of the difference between the angles at the ends of the transmission line
		}
	}	
	double Limit = reacT * ptMax; // Upper limit of the power flow limit scaled by reactance
	double Obj1 = ( Limit - end1A ) * ( Limit - end1A ) + ( Limit + end2C ) * ( Limit + end2C ) + reacT * reacT * ( Thetat1 - end1B ) * ( Thetat1 - end1B ) + reacT * reacT * ( Thetat1 + Limit - end2D ) * ( Thetat1 + Limit - end2D ); // Objective on assumption that difference between angles is equal to upper limit allowed
	double Obj2 = ( -Limit - end1A ) * ( -Limit - end1A ) + ( -Limit + end2C ) * ( -Limit + end2C ) + reacT * reacT * ( Thetat1 - end1B ) * ( Thetat1 - end1B ) + reacT * reacT * ( Thetat1 - Limit - end2D ) * ( Thetat1 - Limit - end2D ); // Objective on assumption that difference between angles is equal to lower limit allowed
	if ( ( Diff <= Limit ) && ( Diff >= -Limit ) ) { // If the power flow and consequently the angle difference is well within allowed limits
		double Obj3 = ( Diff - end1A ) * ( Diff - end1A ) + ( Diff + end2C ) * ( Diff + end2C ) + reacT * reacT * ( Thetat1 - end1B ) * ( Thetat1 - end1B ) + reacT * reacT * ( Thetat1 + Diff - end2D ) * ( Thetat1 + Diff - end2D ); // Objective on assumption that Difference between angles lies well within the allowed limits
		double Obj = ( Obj1 < Obj2 ? ( Obj1 < Obj3 ? Obj1 : Obj3 ) : ( Obj2 < Obj3 ? Obj2 : Obj3 ) );
		if ( Obj == Obj1 ) { // if Diff == Limit gives the lowest objective
			if ( getTranslNodeID2() == 1 ) { // check if end-2 is slack bus
				Thetat2 = 0.0; // in that case fix the corresponding angle to zero
				Thetat1 = Thetat2 - Limit; // adjust the end-1 angle accordingly
			}
			else { // if end-2 is not the slack bus
				Thetat2 = Thetat1 + Limit; // adjust the end-2 angle accordingly
			}
		}
		else {
			if ( Obj == Obj2 ) { // if Diff == -Limit gives the lowest objective
				if ( getTranslNodeID2() == 1 ) { // check if end-2 is slack bus
					Thetat2 = 0.0; // in that case fix the corresponding angle to zero
					Thetat1 = Thetat2 + Limit; // adjust the end-1 angle accordingly
				}
				else { // if end-2 is not the slack bus
					Thetat2 = Thetat1 - Limit; // adjust the end-2 angle accordingly
				}
			}
			else { // if an intermediate value of Diff gives the lowest objective
				if ( getTranslNodeID2() == 1 ) { // check if end-2 is slack bus
					Thetat2 = 0.0; // in that case fix the corresponding angle to zero
					Thetat1 = Thetat2 - Diff; // adjust the end-1 angle accordingly
				}
				else { // if end-2 is not the slack bus
					Thetat2 = Thetat1 + Diff; // adjust the end-2 angle accordingly
				}
			}
		}
	}
	else { // if value of Diff that minimizes the objective falls outside the range
		double Obj = ( Obj1 < Obj2 ? Obj1 : Obj2 ); // check the objective value at the two limit points
		if ( Obj == Obj1 ) { // if Diff == Limit gives the lowest objective
			if ( getTranslNodeID2() == 1 ) { // check if end-2 is the slack bus
				Thetat2 = 0.0; // in that case set the voltage angle of that end to zero
				Thetat1 = Thetat2 - Limit; // adjust the angle of end-1 accordingly
			}
			else { // if end-2 is not the slack bus
				Thetat2 = Thetat1 + Limit; // adjust the end-2 angle accordingly
			}
		}
		else { // if Diff == -Limit gives the lowest objective
			if ( getTranslNodeID2() == 1 ) { // check if end-2 is the slack bus
				Thetat2 = 0.0; // in that case set the voltage angle of that end to zero
				Thetat1 = Thetat2 + Limit; // adjust the angle of end-1 accordingly
			}
			else { // if end-2 is not the slack bus
				Thetat2 = Thetat1 - Limit; // adjust the end-2 angle accordingly
			}
		}
	} // whichever objective is the minimum, consider that value of angle difference as the optimizer
	Pt2 = ( Thetat1 - Thetat2 ) / reacT; // get the transmission line end-2 Power iterate
	Pt1 = ( Thetat2 - Thetat1 ) / reacT; // get the transmission line end-2 voltage angle iterate
	//*cout << "\nFrom Transmission Line " << translID << " in the scenario " << contCount;
	//*cout << "\nUnder Contingency case, " << contCount << " resistance, reactance, and capacity of line " << translID << " respectively, are " << resT << "\t" << reacT << "\t" << ptMax << endl;
	//*cout << "\nPt1 from transmission line solver: " << Pt1 << "\tThetat1 from transmission line solver: " << Thetat1 << endl;
	//*cout << "\nPt2 from transmission line solver: " << Pt2 << "\tThetat2 from transmission line solver: " << Thetat2 << endl;
//} 
	connNodet1Ptr->powerangleMessage( Pt1, &v1, &Thetat1, deviceNature, contCount ); // passes to node object at end 1 the corresponding iterates of power, angle and v
	connNodet2Ptr->powerangleMessage( Pt2, &v2, &Thetat2, deviceNature, contCount ); // passes to node object at end 2 the corresponding iterates of power, angle and v
}
} // function tpowerangleMessage ends

double transmissionLine::translPower1() //const // function translPower1 begins
{
	return Pt1; // returns the Pt1 iterate
} // function translPower1 ends

double transmissionLine::translPower2() //const // function translPower2 begins
{
	return Pt2; // returns the Pt2 iterate
} // function translPower2 ends

double transmissionLine::calcPtilde1( int contin ) //const // function calcPtilde1 begins
{
	double P_avg1 = connNodet1Ptr->PavMessageCont( contin ); // Gets average power for end-1 from the corresponding node object
	double Ptilde1 = Pt1 - P_avg1; // calculates the difference between power iterate and average
	return Ptilde1; // returns the difference
} // function calcPtilde1 ends

double transmissionLine::calcPavInit1() const // function calcPavInit1 begins
{
	return connNodet1Ptr->devpinitMessage(); // seeks the initial Ptilde from the node at end 1
} // function calcPavInit1 ends

double transmissionLine::calcPtilde2( int contin ) //const // function calcPtilde2 begins
{
	double P_avg2 = connNodet2Ptr->PavMessageCont( contin ); // Gets average power for end-2 from the corresponding node object
	double Ptilde2 = Pt2 - P_avg2; // calculates the difference between power iterate and average
	return Ptilde2; // returns the difference
} // function calcPtilde2 ends

double transmissionLine::calcPavInit2() const // function calcPavInit2 begins
{
	return connNodet2Ptr->devpinitMessage(); // seeks the initial Ptilde from the node at end 2
} // function calcPavInit2 ends


double transmissionLine::getu1( int contin ) const // function getu1 begins
{
	if ( ptMax == pow( 10.0, -100.0 ) ) {
		double u1 = 0.0;
		return u1;
	}
	else {
		double u1 = connNodet1Ptr->uMessage( contin ); // gets the value of the price corresponding to power balance from node
		//cout << "u1: " << u1 << endl;
		return u1; // returns the price
	}
} // function getu1 ends

double transmissionLine::getu2( int contin ) const // function getu2 begins
{
	if ( ptMax == pow( 10.0, -100.0 ) ) {
		double u2 = 0.0;
		return u2;
	}
	else {
		double u2 = connNodet2Ptr->uMessage( contin ); // gets the value of the price corresponding to power balance from node
		//cout << "u2: " << u2 << endl;
		return u2; // returns the price
	}
} // function getu2 ends

double transmissionLine::calcThetatilde1( int contin ) // function calcThetatilde1 begins
{
	if ( ptMax == pow( 10.0, -100.0 ) ) {
		double Theta_tilde1 = 0.0; // claculate the deviation between the voltage angle of the device and the average
		//*cout << "\nTransmission line # " << translID << " contingency scenario# " << contin << " end-1 angle = " << Thetat1 << " , end-1 average angle = " << Theta_avg1 << endl;
		return Theta_tilde1; // return the deviation
	}
	else {
		double Theta_avg1 = connNodet1Ptr->ThetaavMessageCont( contin ); // get the average voltage angle at the particular node
		double Theta_tilde1 = Thetat1 - Theta_avg1; // claculate the deviation between the voltage angle of the device and the average
		//*cout << "\nTransmission line # " << translID << " contingency scenario# " << contin << " end-1 angle = " << Thetat1 << " , end-1 average angle = " << Theta_avg1 << endl;
		return Theta_tilde1; // return the deviation
	}
} // function calcThetatilde1 ends

double transmissionLine::calcThetatilde2( int contin ) // function calcThetatilde2 begins
{
	if ( ptMax == pow( 10.0, -100.0 ) ) {
		double Theta_tilde2 = 0.0; // claculate the deviation between the voltage angle of the device and the average
		//*cout << "\nTransmission line # " << translID << " contingency scenario# " << contin << " end-1 angle = " << Thetat1 << " , end-1 average angle = " << Theta_avg1 << endl;
		return Theta_tilde2; // return the deviation
	}
	else {
		double Theta_avg2 = connNodet2Ptr->ThetaavMessageCont( contin ); // get the average voltage angle at the particular node
		double Theta_tilde2 = Thetat2 - Theta_avg2; // claculate the deviation between the voltage angle of the device and the average
		//*cout << "\nTransmission line # " << translID << " contingency scenario# " << contin << " end-2 angle = " << Thetat2 << " , end-2 average angle = " << Theta_avg2 << endl;
		return Theta_tilde2; // return the deviation
	}
} // function calcThetatilde2 ends

double transmissionLine::calcvtilde1() const // function calcvtilde1 begins
{
	double v_avg1 = connNodet1Ptr->vavMessage(); // get the average of the Lagrange multiplier corresponding to voltage angle balance 
	double v_tilde1 = v1 - v_avg1; // calculate the deviation of the node Lagrange multiplier to the average
	return v_tilde1; // return the deviation
} // function calcvtilde1 ends

double transmissionLine::calcvtilde2() const // function calcvtilde2 begins
{
	double v_avg2 = connNodet2Ptr->vavMessage(); // get the average of the Lagrange multiplier corresponding to voltage angle balance 
	double v_tilde2 = v2 - v_avg2; // calculate the deviation of the node Lagrange multiplier to the average
	return v_tilde2; // return the deviation
} // function calcvtilde2 ends

double transmissionLine::getv1( int contin ) //const // function getv1 begins
{
	if ( ptMax == pow( 10.0, -100.0 ) ) {
		v1 = 0.0; // claculate the deviation between the voltage angle of the device and the average
		//*cout << "\nTransmission line # " << translID << " contingency scenario# " << contin << " end-1 angle = " << Thetat1 << " , end-1 average angle = " << Theta_avg1 << endl;
		return v1; // return the deviation
	}
	else {
		//cout << "v1_initial: " << v1 << endl;
		v1 = v1 + calcThetatilde1( contin ); // Calculate the value of the Lagrange multiplier corresponding to angle constraint
		//cout << "v1_final: " << v1 << endl;
		return v1; // return the voltage angle price
	}
} // function getv1 ends

double transmissionLine::getv2( int contin ) //const // function getv2 begins
{
	if ( ptMax == pow( 10.0, -100.0 ) ) {
		v2 = 0.0; // claculate the deviation between the voltage angle of the device and the average
		//*cout << "\nTransmission line # " << translID << " contingency scenario# " << contin << " end-1 angle = " << Thetat1 << " , end-1 average angle = " << Theta_avg1 << endl;
		return v2; // return the deviation
	}
	else {
		//cout << "v2_initial: " << v2 << endl;
		v2 = v2 + calcThetatilde2( contin ); // Calculate the value of the Lagrange multiplier corresponding to angle constraint
		//cout << "v2_final: " << v2 << endl;
		return v2; // Calculate the value of the Lagrange multiplier corresponding to angle constraint
	}
} // function getv2 ends	


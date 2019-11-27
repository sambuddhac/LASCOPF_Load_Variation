// Supernetwork source code for implementing APMP (Auxiliary Proximal Message Passing) Algorithm for the SCOPF in serial mode
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <cmath>
#include <fstream>
#include <cstring>
using namespace std;
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file
#include "supernetwork.h" // Network class definition
#include "network.h" // Network class definition

superNetwork::superNetwork(int networkID, int choiceSolver, int rhoTuning, int dispInterval, int lastFlag, int nextChoice) // function main begins program execution
{
	netID=networkID;
	solverChoice=choiceSolver; 
	setRhoTuning=rhoTuning; // parameter to select adaptive rho, fixed rho, and type of adaptive rho}
	intervalCount=dispInterval; // count of the dispatch interval to which the particular network instance for the coarse grain belongs
	lastInterval=lastFlag; // Flas to indicate if the network belongs to last interval: 0=not last interval; 1=last interval
	cout << endl << "\n*** NETWORK INITIALIZATION STAGE BEGINS ***\n" << endl << endl;
	Network* network = new Network( netID, 0, 0, 0, solverChoice, intervalCount, lastInterval, nextChoice ); // create network object corresponding to the base case
	numberOfCont = network->retContCount(); // gets the number of contingency scenarios in the variable numberOfCont
	contNetVector.push_back( network ); // push to the vector of network instances 
	for ( int i = 1; i <= numberOfCont; ++i ) {
		int lineOutaged = contNetVector[0]->indexOfLineOut(i); // gets the serial number of transmission line outaged in this scenario 
		Network* network = new Network( netID, i, lineOutaged, 1, solverChoice, intervalCount, lastInterval, nextChoice ); // create the network instances for the contingency scenarios, which includes as many networks as the number of contingency scenarios
		contNetVector.push_back( network ); // push to the vector of network instances
	}

	cout << "\n*** NETWORK INITIALIZATION STAGE ENDS ***\n" << endl;

	numberOfGenerators = contNetVector[0]->getGenNumber(); // get the number of generators in the system
	consLagDim = numberOfCont*numberOfGenerators; // Dimension of the vectors of APP Lagrange Multipliers and Power Generation Consensus
}

superNetwork::~superNetwork() // Destructor
{
	cout << "\nDispatch interval super-network object for dispatch interval " << intervalCount << " destroyed" << endl;
}

double superNetwork::getvirtualNetExecTime(){return virtualNetExecTime;}

void superNetwork::runSimulation(int outerIter, double LambdaOuter[], double powDiffOuter[], double powSelfBel[], double powNextBel[], double powPrevBel[], int innerAPP, GRBEnv* environmentGUROBI) { // runs the distributed SCOPF simulations using ADMM-PMP with CVXGEN custom solver	
	double lambdaAPP[consLagDim]; // Array of APP Lagrange Multipliers for achieving consensus among the values of power generated, as guessed by scenarios
	double powDiff[consLagDim]; // Array of lack of consensus between generation values, as guessed by scenarios
	double alphaAPP[consLagDim]; // APP Parameter/Path-length
	double alphaAPP1[consLagDim]; // Previous value of APP Parameter/Path-length from previous iteration
	double W[consLagDim], Wprev[consLagDim]; // Present and previous values of W for the PID controller for modifying Rho
	double lambdaAdap = 0.1; // Parameter of the Proportional (P) controller for adjusting the APP tuning parameter
	double muAdap = 0.5; // Parameter of the Derivative (D) controller for adjusting the APP tuning parameter
        double xiAdap = 0.0000; // Parameter of the Integral (I) controller for adjusting the APP tuning parameter
        double controllerSum[consLagDim]; // Integral term of the PID controller
	int setInnerAPPTuning=innerAPP; // parameter to select adaptive, fixed, and type of adaptive APP step length
	for ( int i = 0; i < consLagDim; ++i ) {
		lambdaAPP[i] = 0.0; // Initialize lambdaAPP for the first iteration of APP and ADMM-PMP
		powDiff[i] = 0.0; // Initialize powDiff for the first iteration of APP and ADMM-PMP
		alphaAPP[i] = 100.0; // Initialize APP Parameter/Path-length
		controllerSum[i] = 0.0; // Initialize APP Integral term
		W[i] = 0.0;
		Wprev[i] = 0.0;
	}
	iterCountAPP = 1; // Iteration counter for APP coarse grain decomposition algorithm
	finTol = 1000.0; // Initial Guess of the Final tolerance of the APP iteration/Stopping criterion
	if ((solverChoice==1) || (solverChoice==2)) { // APMP Fully distributed, Bi-layer (N-1) SCOPF Simulation 
		string outputAPPFileName;
		if (solverChoice==1)
			outputAPPFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/output/ADMM_PMP_GUROBI/resultAPP-SCOPF"+to_string(intervalCount)+".txt";
		if (solverChoice==2)
			outputAPPFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/output/ADMM_PMP_CVXGEN/resultAPP-SCOPF"+to_string(intervalCount)+".txt";
		ofstream matrixResultAPPOut( outputAPPFileName, ios::out ); // create a new file result.txt to output the results

		// exit program if unable to create file
		if ( !matrixResultAPPOut ) {
			cerr << "File could not be opened" << endl;
			exit( 1 );
		}
		matrixResultAPPOut << endl << "\nInitial Value of the Tolerance to kick-start the APP outer iterations= " << finTol << "\n" << endl << endl;
		matrixResultAPPOut << "APP Iteration Count" << "\t" << "APP Tolerance" << "\n";	
		clock_t start_s = clock(); // begin keeping track of the time
		cout << endl << "\n*** APMP ALGORITHM BASED COARSE+FINE GRAINED BILAYER DECENTRALIZED/DISTRIBUTED SCOPF (SERIAL IMPLEMENTATION) BEGINS ***\n" << endl << endl;
		cout << endl << "\n*** SIMULATION IN PROGRESS; PLEASE DON'T CLOSE ANY WINDOW OR OPEN ANY OUTPUT FILE YET ... ***\n" << endl << endl;

		//*********************************************AUXILIARY PROBLEM PRINCIPLE (APP) COARSE GRAINED DECOMPOSITION COMPONENT******************************************************//	
		//do { // APP Coarse grain iterations start
		largestNetTimeVec.clear();
		double actualNetTime = 0;
		//do {
		for ( iterCountAPP = 1; ((iterCountAPP <= 10) /*&& (finTol>=0.7)*/); ++iterCountAPP ) {
			singleNetTimeVec.clear();
			for ( int netSimCount = 0; netSimCount < (numberOfCont+1); ++netSimCount ) {
				cout << "\nStart of " << iterCountAPP << " -th Innermost APP iteration for " << netSimCount+1 << " -th base/contingency scenario" << endl;
				contNetVector[netSimCount]->runSimulation(outerIter, LambdaOuter, powDiffOuter, setRhoTuning, iterCountAPP, lambdaAPP, powDiff, powSelfBel, powNextBel, powPrevBel, environmentGUROBI); // start simulation
				double singleNetTime = contNetVector[netSimCount]->returnVirtualExecTime();
				actualNetTime += singleNetTime;
				singleNetTimeVec.push_back(singleNetTime);
			}
			double largestNetTime = *max_element(singleNetTimeVec.begin(), singleNetTimeVec.end());
			largestNetTimeVec.push_back(largestNetTime);
			for ( int i = 0; i < numberOfCont; ++i ) {
				for ( int j = 0; j < numberOfGenerators; ++j ) {
					powDiff[i*numberOfGenerators+j]=*(contNetVector[0]->getPowSelf()+j)-*(contNetVector[i+1]->getPowSelf()+j); // what base thinks about itself Vs. what contingency thinks about base
				}
			}
			// Tuning the alphaAPP by a discrete-time PID Controller
			for ( int i = 0; i < numberOfCont; ++i ) {
				for ( int j = 0; j < numberOfGenerators; ++j ) {
					if ( ( iterCountAPP <= 10 ) && (setInnerAPPTuning!=0) ) {
						W[i*numberOfGenerators+j] = abs(( *(contNetVector[i+1]->getPowSelf()+j) / *(contNetVector[0]->getPowSelf()+j) ) - 1); // Definition of W for adaptive path-length
					}
					if ( ( iterCountAPP > 5 ) && ( iterCountAPP <= 10 ) && (setInnerAPPTuning==0) ) {
						W[i*numberOfGenerators+j] = 0.0; // Definition of W for adaptive path-length
						alphaAPP[i*numberOfGenerators+j] = 75.0;
					}
					if ( ( iterCountAPP > 10 ) && ( iterCountAPP <= 15 ) && (setInnerAPPTuning==0) ) {
						W[i*numberOfGenerators+j] = 0.0; // Definition of W for adaptive path-length
						alphaAPP[i*numberOfGenerators+j] = 2.5;
					}
					if ( ( iterCountAPP > 15 ) && ( iterCountAPP <= 20 ) && (setInnerAPPTuning==0) ) {
						W[i*numberOfGenerators+j] = 0.0; // Definition of W for adaptive path-length
						alphaAPP[i*numberOfGenerators+j] = 1.25;
					}
					if ( ( iterCountAPP > 20 ) && (setInnerAPPTuning==0) ) {
						W[i*numberOfGenerators+j] = 0.0; // Definition of W for fixed Rho
						alphaAPP[i*numberOfGenerators+j] = 0.5;
					}
					// Calculation of Adaptive Step-length
					controllerSum[i*numberOfGenerators+j] = controllerSum[i*numberOfGenerators+j] + W[i*numberOfGenerators+j];
					alphaAPP1[i*numberOfGenerators+j] = alphaAPP[i*numberOfGenerators+j]; // Store previous alphaAPP
					alphaAPP[i*numberOfGenerators+j] = ( alphaAPP1[i*numberOfGenerators+j] ) * ( exp( ( lambdaAdap * W[i*numberOfGenerators+j] ) + ( muAdap * ( W[i*numberOfGenerators+j] - Wprev[i*numberOfGenerators+j] ) ) + ( xiAdap * controllerSum[i*numberOfGenerators+j]  ) ) ); // Next iterate value of Rho
					Wprev[i*numberOfGenerators+j] = W[i*numberOfGenerators+j]; // Buffering
				}
			}
			for ( int i = 0; i < numberOfCont; ++i ) {
				for ( int j = 0; j < numberOfGenerators; ++j ) {
					lambdaAPP[i*numberOfGenerators+j] = lambdaAPP[i*numberOfGenerators+j] + alphaAPP[i*numberOfGenerators+j] * (powDiff[i*numberOfGenerators+j]); // what I think about myself Vs. what next door fellow thinks about me
				}
			}
			double tolAPP = 0.0;
			for ( int i = 0; i < consLagDim; ++i ) {
				tolAPP = tolAPP + pow(powDiff[i], 2);
			}
			finTol = sqrt(tolAPP);
			matrixResultAPPOut << iterCountAPP << "\t" << finTol << "\n";
			//++iterCountAPP; // increment the APP iteration counter
		//} while (finTol>=0.5); //Check the termination criterion of the APP iterations
		}
		//****************************************END OF AUXILIARY PROBLEM PRINCIPLE (APP) COARSE GRAINED DECOMPOSITION COMPONENT******************************************************//

		cout << "\n*** SCOPF SIMULATION ENDS ***\n" << endl;
		clock_t stop_s = clock();  // end
		matrixResultAPPOut << "\nExecution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;
		matrixResultAPPOut << "\nVirtual Supernetwork Execution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC  - actualNetTime + accumulate(largestNetTimeVec.begin(), largestNetTimeVec.end(), 0.0)<< endl;
		virtualNetExecTime = static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC  - actualNetTime + accumulate(largestNetTimeVec.begin(), largestNetTimeVec.end(), 0.0);
		//cout << "\nExecution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;
		cout << "\nFinal Value of APP Tolerance " << finTol << endl; 
	}

	if (solverChoice==3) { // Centralized (N-1) SCOPF Simulation
		string outputAPPFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/output/APP_Quasi_Decent_GUROBI/resultAPP-SCOPF"+to_string(intervalCount)+".txt";
		ofstream matrixResultAPPOut( outputAPPFileName, ios::out ); // create a new file result.txt to output the results

		// exit program if unable to create file
		if ( !matrixResultAPPOut ) {
			cerr << "File could not be opened" << endl;
			exit( 1 );
		}
		matrixResultAPPOut << endl << "\nInitial Value of the Tolerance to kick-start the APP outer iterations= " << finTol << "\n" << endl << endl;
		matrixResultAPPOut << "APP Itearion Count" << "\t" << "APP Tolerance" << "\n";	
		clock_t start_s = clock(); // begin keeping track of the time
		cout << endl << "\n*** COARSE GRAINED APP QUASI-DECENTRALIZED (N-1) SCOPF SIMULATION BEGINS ***\n" << endl << endl;
		cout << endl << "\n*** SIMULATION IN PROGRESS; PLEASE WAIT. PLEASE DON'T CLOSE ANY WINDOW OR OPEN ANY OUTPUT FILE ... ***\n" << endl << endl;
		//*********************************************AUXILIARY PROBLEM PRINCIPLE (APP) COARSE GRAINED DECOMPOSITION COMPONENT******************************************************/	
		//do { // APP Coarse grain iterations start
		largestNetTimeVec.clear();
		double actualNetTime = 0;
		for ( iterCountAPP = 1; iterCountAPP <= 2; ++iterCountAPP ) {
			singleNetTimeVec.clear();
			contNetVector[0]->runSimAPPGurobiBase(outerIter, LambdaOuter, powDiffOuter, iterCountAPP, lambdaAPP, powDiff, powSelfBel, powNextBel, powPrevBel, environmentGUROBI); // start simulation
			double singleNetTime = contNetVector[0]->returnVirtualExecTime();
			actualNetTime += singleNetTime;
			singleNetTimeVec.push_back(singleNetTime);
			for ( int netSimCount = 1; netSimCount < (numberOfCont+1); ++netSimCount ) {
				contNetVector[netSimCount]->runSimAPPGurobiCont(outerIter, LambdaOuter, powDiffOuter, iterCountAPP, lambdaAPP, powDiff, environmentGUROBI); // start simulation
				singleNetTime = contNetVector[netSimCount]->returnVirtualExecTime();
				actualNetTime += singleNetTime;
				singleNetTimeVec.push_back(singleNetTime);
			}
			double largestNetTime = *max_element(singleNetTimeVec.begin(), singleNetTimeVec.end());
			largestNetTimeVec.push_back(largestNetTime);
			for ( int i = 0; i < numberOfCont; ++i ) {
				for ( int j = 0; j < numberOfGenerators; ++j ) {
					powDiff[i*numberOfGenerators+j]=*(contNetVector[0]->getPowSelfGUROBI()+j)-*(contNetVector[i+1]->getPowSelfGUROBI()+j); // what base thinks about itself Vs. what contingency thinks about base
				}
			}
			// Tuning the alphaAPP by a discrete-time PID Controller
			for ( int i = 0; i < numberOfCont; ++i ) {
				for ( int j = 0; j < numberOfGenerators; ++j ) {
					if ( ( iterCountAPP <= 100 ) && (setInnerAPPTuning!=0) ) {
						W[i*numberOfGenerators+j] = abs(( *(contNetVector[i+1]->getPowSelf()+j) / *(contNetVector[0]->getPowSelf()+j) ) - 1); // Definition of W for adaptive path-length
					}
					else {
						W[i*numberOfGenerators+j] = 0.0; // Definition of W for fixed Rho
					}
					// Calculation of Adaptive Step-length
					controllerSum[i*numberOfGenerators+j] = controllerSum[i*numberOfGenerators+j] + W[i*numberOfGenerators+j];
					alphaAPP1[i*numberOfGenerators+j] = alphaAPP[i*numberOfGenerators+j]; // Store previous alphaAPP
					alphaAPP[i*numberOfGenerators+j] = ( alphaAPP1[i*numberOfGenerators+j] ) * ( exp( ( lambdaAdap * W[i*numberOfGenerators+j] ) + ( muAdap * ( W[i*numberOfGenerators+j] - Wprev[i*numberOfGenerators+j] ) ) + ( xiAdap * controllerSum[i*numberOfGenerators+j]  ) ) ); // Next iterate value of Rho
					Wprev[i*numberOfGenerators+j] = W[i*numberOfGenerators+j]; // Buffering
				}
			}
			for ( int i = 0; i < numberOfCont; ++i ) {
				for ( int j = 0; j < numberOfGenerators; ++j ) {
					lambdaAPP[i*numberOfGenerators+j] = lambdaAPP[i*numberOfGenerators+j] + alphaAPP[i*numberOfGenerators+j] * (powDiff[i*numberOfGenerators+j]); // what I think about myself Vs. what next door fellow thinks about me
				}
			}
			double tolAPP = 0.0;
			for ( int i = 0; i < consLagDim; ++i ) {
				tolAPP = tolAPP + pow(powDiff[i], 2);
			}
			finTol = sqrt(tolAPP);
			matrixResultAPPOut << iterCountAPP << "\t" << finTol << "\n";
			//++iterCountAPP; // increment the APP iteration counter
		//} while (finTol>=0.05); //Check the termination criterion of the APP iterations
		}
		//****************************************END OF AUXILIARY PROBLEM PRINCIPLE (APP) COARSE GRAINED DECOMPOSITION COMPONENT******************************************************/
		cout << "\n*** SCOPF SIMULATION ENDS ***\n" << endl;
		clock_t stop_s = clock();  // end
		matrixResultAPPOut << "\nExecution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;
		cout << "\nExecution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;
		cout << "\nFinal Value of APP Tolerance " << finTol << endl; 
		cout << "\n*** (N-1) SCOPF SIMULATION ENDS ***\n" << endl;	
	}
	if (solverChoice==4) { // Centralized (N-1) SCOPF Simulation
		contNetVector[0]->runSimulationCentral(outerIter, LambdaOuter, powDiffOuter, powSelfBel, powNextBel, powPrevBel, environmentGUROBI); // start simulation
	}
	else
		cout << "\nInvalid choice of solution method and algorithm." << endl;
} // end main

double *superNetwork::getPowSelf()
{
	if ((solverChoice==1) || (solverChoice==2))
		return (contNetVector[0])->getPowSelf();
	if ((solverChoice==3) || (solverChoice==4))
		return (contNetVector[0])->getPowSelfGUROBI();
} // returns the difference in the values of what I think about myself Vs. what next door fellow thinks about me

double *superNetwork::getPowPrev()
{
	if ((solverChoice==1) || (solverChoice==2))
		return (contNetVector[0])->getPowPrev();
	if ((solverChoice==3) || (solverChoice==4))
		return (contNetVector[0])->getPowPrevGUROBI();
} // returns what I think about previous dispatch interval generators

double *superNetwork::getPowNext()
{
	if ((solverChoice==1) || (solverChoice==2))
		return (contNetVector[0])->getPowNext();
	if ((solverChoice==3) || (solverChoice==4))
		return (contNetVector[0])->getPowNextGUROBI();
} // returns what I think about next door fellow 

int superNetwork::getGenNumber() //Function getGenNumber begins
{
	return numberOfGenerators;
} // end of getGenNumber function

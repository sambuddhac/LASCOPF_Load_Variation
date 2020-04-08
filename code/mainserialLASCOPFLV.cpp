// Main function for implementing APMP Algorithm for the LASCOPF for Load (Demand) Variation Tracking case in serial mode
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <cmath>
#include <fstream>
#include <cstring>
#include "gurobi_c++.h"
#include "omp.h"
using namespace std;
#include "supernetwork.h" // Super Network class definition
#define NUM_THREADS 4
int main() // function main begins program execution
{
	int netID; // Network ID number to indicate the type of the system with specifying the number of buses/nodes
	int solverChoice; // Solver choice among CVXGEN-ADMM-PMP+APP fully distributed, GUROBI-ADMM-PMP+APP fully distributed, or GUROBI APP half distributed
	int dispatchIntervals; // Number of dispatch intervals for the look-ahead simulation
	int nextChoice; // The index which decides whether to include the ramping constraint for the last supernetwork for considering the next imaginary interval or not
	int last = 0; // flag to indicate the last interval; last = 0, for dispatch interval that is not the last one; last = 1, for the last interval
	vector< superNetwork* > futureNetVector; // Vector of future look-ahead dispatch interval network objects
	cout << "\nEnter the number of nodes to initialize the network. (Allowed choices are 2, 3, 5, 14, 30, 48, 57, 118, and 300 Bus IEEE Test Bus Systems as of now. So, please restrict yourself to one of these)\n";
	cin >> netID;
	cout << "\nEnter the choice of the solver for SCOPF of each dispatch interval, 1 for GUROBI-APMP(ADMM/PMP+APP), 2 for CVXGEN-APMP(ADMM/PMP+APP), 3 for GUROBI APP Coarse Grained, 4 for centralized GUROBI SCOPF" << endl;
	cin >> solverChoice;
	cout << "\nEnter the choice pertaining to whether you want to consider the ramping constraint to the next interval, for the last interval: 0 for not considering and 1 for considering" << endl;
	cin >> nextChoice;
	int setRhoTuning; // parameter to select adaptive rho, fixed rho, and type of adaptive rho
	if ((solverChoice==1) || (solverChoice==2)) { // APMP Fully distributed, Bi-layer (N-1) SCOPF Simulation 
		cout << "Enter the tuning mode; Enter 1 for maintaining Rho * primTol = dualTol; 2 for primTol = dualTol; anything else for Adaptive Rho (with mode-1 being implemented for the first 3000 iterations and then Rho is held constant).\n" << endl;
		cin >> setRhoTuning;
	}
	else
		setRhoTuning = 0; // Otherwise, if we aren't using ADMM-PMP, Rho tuning is unnecessary, 0 is a dummy value
	cout << "\nEnter the number of look-ahead dispatch time horizons.\n";
	cin >> dispatchIntervals;

	cout << endl << "\n*** SUPERNETWORK INITIALIZATION STAGE BEGINS ***\n" << endl << endl;
	GRBEnv* environmentGUROBI = new GRBEnv("GUROBILogFile.log"); // GUROBI Environment object for storing the different optimization models
	for ( int i = 0; i < dispatchIntervals; ++i ) {
		if (i==(dispatchIntervals-1))
			last = 1; // set the flag to 1 to indicate the last interval
		superNetwork* supernet = new superNetwork( netID, solverChoice, setRhoTuning, i, last, nextChoice ); // create the network instances for the future dispatch intervals
		futureNetVector.push_back( supernet ); // push to the vector of future network instances
	}

	cout << "\n*** SUPERNETWORK INITIALIZATION STAGE ENDS ***\n" << endl;

	int numberOfGenerators = futureNetVector[0]->getGenNumber(); // get the number of generators in the system
	int iterCountAPP = 1; // Iteration counter for APP coarse grain decomposition algorithm
	int consLagDim = 2*(dispatchIntervals-1)*numberOfGenerators; // Dimension of the vectors of APP Lagrange Multipliers and Power Generation Consensus	
	double lambdaAPP[consLagDim]; // Array of APP Lagrange Multipliers for achieving consensus among the values of power generated, as guessed by different intervals
	double powDiff[consLagDim]; // Array of lack of consensus between generation values, as guessed by different intervals
	double alphaAPP[consLagDim]; // APP Parameter/Path-length
	double alphaAPP1[consLagDim]; // Previous value of APP Parameter/Path-length from previous iteration
	double W[consLagDim], Wprev[consLagDim]; // Present and previous values of W for the PID controller for modifying Rho
	double lambdaAdap = 0.1; // Parameter of the Proportional (P) controller for adjusting the APP tuning parameter
	double muAdap = 0.5; // Parameter of the Derivative (D) controller for adjusting the APP tuning parameter
        double xiAdap = 0.0000; // Parameter of the Integral (I) controller for adjusting the APP tuning parameter
        double controllerSum[consLagDim]; // Integral term of the PID controller
	int setOuterAPPTuning; // parameter to select adaptive, fixed, and type of adaptive APP step length
	int setInnerAPPTuning; // parameter to select adaptive, fixed, and type of adaptive APP step length
	double powerSelfGen[dispatchIntervals*numberOfGenerators]; // what I think about myself
	double powerNextBel[dispatchIntervals*numberOfGenerators]; // what I think about next door fellow
	double powerPrevBel[dispatchIntervals*numberOfGenerators]; // what I think about previous door fellow
	for ( int i = 0; i < consLagDim; ++i ) {
		lambdaAPP[i] = 0.0; // Initialize lambdaAPP for the first iteration of APP and ADMM-PMP
		powDiff[i] = 0.0; // Initialize powDiff for the first iteration of APP and ADMM-PMP
		alphaAPP[i] = 10.0; // Initialize APP Parameter/Path-length
		controllerSum[i] = 0.0; // Initialize APP Integral term
		W[i] = 0.0;
		Wprev[i] = 0.0;
	}
	cout << "\nEnter the choice for tuning the step-length for the external APP iterations for consensus among the different intervals' MW outputs: 0 for constant step-length and anything else for a discrete time PID controller for tuning step-length.\n";
	cin >> setOuterAPPTuning;
	cout << "\nEnter the choice for tuning the step-length for the internal APP iterations for consensus among the different scenarios' MW outputs: 0 for constant step-length and anything else for a discrete time PID controller for tuning step-length.\n" << endl;
	cin >> setInnerAPPTuning;
	// Initializing the self belief, next belief, and previous beliefs about MW generated by a warm start with the respective generation values of last realized dispatch
	for ( int i = 0; i < dispatchIntervals; ++i ) {
		for ( int j = 0; j < numberOfGenerators; ++j ) {
			powerSelfGen[i*numberOfGenerators+j] = *((futureNetVector[0])->getPowPrev()+j); // Use 0.0 if warm start is not desired
			powerNextBel[i*numberOfGenerators+j] = *((futureNetVector[0])->getPowPrev()+j); // Use 0.0 if warm start is not desired
			if (i==0) {
				powerPrevBel[i*numberOfGenerators+j] = *((futureNetVector[0])->getPowPrev()+j); // Actual value of previous interval dispatch for the first interval
			}
			else {
				powerPrevBel[i*numberOfGenerators+j] = *((futureNetVector[0])->getPowPrev()+j); // Use 0.0 if warm start is not desired
			}
		}
	}
	double finTol = 1000.0; // Initial Guess of the Final tolerance of the APP iteration/Stopping criterion
	double finTolDelayed = 1000.0; // Initial Guess of the Final tolerance delayed of the APP iteration/Stopping criterion
	string outputAPPFileName;
	if (solverChoice==1)
		outputAPPFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/output/ADMM_PMP_GUROBI/resultOuterAPP-SCOPF.txt";
	if (solverChoice==2)
		outputAPPFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/output/ADMM_PMP_CVXGEN/resultOuterAPP-SCOPF.txt";
	if (solverChoice==3)
		outputAPPFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/output/APP_Quasi_Decent_GUROBI/resultOuterAPP-SCOPF.txt";
	if (solverChoice==4)
		outputAPPFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/LASCOPF_Load_Variation/output/APP_GUROBI_Centralized_SCOPF/resultOuterAPP-SCOPF.txt";
	ofstream matrixResultAPPOut( outputAPPFileName, ios::out ); // create a new file to output the results
	// exit program if unable to create file
	if ( !matrixResultAPPOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	matrixResultAPPOut << endl << "\n*** APMP ALGORITHM BASED LASCOPF FOR LOAD (DEMAND) VARIATION TRACKING SIMULATION (SERIAL IMPLEMENTATION) SUPERNETWORK LAYER BEGINS ***\n" << endl << endl;
	matrixResultAPPOut << endl << "\n*** SIMULATION IN PROGRESS; PLEASE DON'T CLOSE ANY WINDOW OR OPEN ANY OUTPUT FILE YET ... ***\n" << endl << endl;	
	matrixResultAPPOut << endl << "\nInitial Value of the Tolerance to kick-start the APP outer iterations= " << finTol << "\n" << endl << endl;
	matrixResultAPPOut << "APP Itearion Count" << "\t" << "APP Tolerance" << "\n";	
	clock_t start_s = clock(); // begin keeping track of the time
	cout << endl << "\n*** APMP ALGORITHM BASED LASCOPF FOR LOAD (DEMAND) VARIATION TRACKING SIMULATION (SERIAL IMPLEMENTATION) SUPERNETWORK LAYER BEGINS ***\n" << endl << endl;
	cout << endl << "\n*** SIMULATION IN PROGRESS; PLEASE DON'T CLOSE ANY WINDOW OR OPEN ANY OUTPUT FILE YET ... ***\n" << endl << endl;

//*********************************************AUXILIARY PROBLEM PRINCIPLE (APP) COARSE GRAINED DECOMPOSITION COMPONENT******************************************************//
	//do { // APP Coarse grain iterations start
	vector<double> largestSuperNetTimeVec; // vector largest value of the computational time in a particular outer APP iteration for any supernetwork
	vector<double> singleSuperNetTimeVec; // vector of the computational times in a particular outer APP iteration for all supernetworks
	largestSuperNetTimeVec.clear(); // clear for upcoming iteration
	double actualSuperNetTime = 0; // Initialize the supernetwork computational time
	do { // begin the outer AP iterations for attaining consensus among the MW values of different intervals
	//for ( iterCountAPP = 1; iterCountAPP <= 10; ++iterCountAPP ) { // begin the outer AP iterations for attaining consensus among the MW values of different intervals
		singleSuperNetTimeVec.clear(); // clear for upcoming iteration
		int nthreads;
    	double time_start, time_finish;
    	omp_set_num_threads(NUM_THREADS);
    	time_start = omp_get_wtime();
		#pragma omp parallel 
		{
        	int ID = omp_get_thread_num();
        	int nthrds = omp_get_num_threads();
        	if(ID==0) nthreads = nthrds;
			for ( int netSimCount = ID; netSimCount < dispatchIntervals; netSimCount += nthrds) { // Solve the individual SCOPFs for different dispatch intervals
				cout << "\nStart of " << iterCountAPP << " -th Outermost APP iteration for " << netSimCount+1 << " -th dispatch interval" << endl;
				futureNetVector[netSimCount]->runSimulation(iterCountAPP, lambdaAPP, powDiff, powerSelfGen, powerNextBel, powerPrevBel, setInnerAPPTuning, environmentGUROBI); // start simulation
				double singleSuperNetTime = futureNetVector[netSimCount]->getvirtualNetExecTime(); // get the computational time for each supernetwork under the assumption of nested and complete parallelism of each generator optimization, within each coarse grain optimization in the supernetworks
				#pragma omp critical
				{
					actualSuperNetTime += singleSuperNetTime; // Actual time
					singleSuperNetTimeVec.push_back(singleSuperNetTime); // Vector of all independent supernet solve times
				}
			}
		}
		double largestSuperNetTime = *max_element(singleSuperNetTimeVec.begin(), singleSuperNetTimeVec.end()); // get the laziest solve-time for this iteration
		largestSuperNetTimeVec.push_back(largestSuperNetTime); // vector of all te laziest supernet calculations over all iterations
		// Calculate the power generation opinions and disagreements between the different dispatch interval coarse grains
		for ( int i = 0; i < ( dispatchIntervals - 1 ); ++i ) {
			for ( int j = 0; j < numberOfGenerators; ++j ) {
				powDiff[2*i*numberOfGenerators+j]=*(futureNetVector[i]->getPowSelf()+j)-*(futureNetVector[i+1]->getPowPrev()+j); // what I think about myself Vs. what next door fellow thinks about me
				powerSelfGen[i*numberOfGenerators+j]=*(futureNetVector[i]->getPowSelf()+j); // what I think about myself
				powerNextBel[i*numberOfGenerators+j]=*(futureNetVector[i]->getPowNext()+j); // what I think about next door fellow
				powerPrevBel[i*numberOfGenerators+j]=*(futureNetVector[i]->getPowPrev()+j); // what I think about previous interval
				powDiff[(2*i+1)*numberOfGenerators+j]=*(futureNetVector[i]->getPowNext()+j)-*(futureNetVector[i+1]->getPowSelf()+j); // what I think about next door fellow Vs. what next door fellow thinks about himself
			}
		}
		for ( int j = 0; j < numberOfGenerators; ++j ) {
			powerSelfGen[( dispatchIntervals - 1 )*numberOfGenerators+j]=*(futureNetVector[( dispatchIntervals - 1 )]->getPowSelf()+j); // what I think about myself
			powerNextBel[( dispatchIntervals - 1 )*numberOfGenerators+j]=*(futureNetVector[( dispatchIntervals - 1 )]->getPowNext()+j); // what I think about next door fellow
			powerPrevBel[( dispatchIntervals - 1 )*numberOfGenerators+j]=*(futureNetVector[( dispatchIntervals - 1 )]->getPowPrev()+j); // what I think about previous interval
		}
		// Tuning the alphaAPP by a discrete-time PID Controller
		for ( int i = 0; i < ( dispatchIntervals - 1 ); ++i ) {
			for ( int j = 0; j < numberOfGenerators; ++j ) {
				if ( ( iterCountAPP <= 1 ) && (setOuterAPPTuning!=0) ) {
					W[2*i*numberOfGenerators+j] = abs(( *(futureNetVector[i+1]->getPowPrev()+j) / *(futureNetVector[i]->getPowSelf()+j) ) - 1); // Definition of W for adaptive path-length
					W[(2*i+1)*numberOfGenerators+j] = abs(( *(futureNetVector[i+1]->getPowSelf()+j) / *(futureNetVector[i]->getPowNext()+j) ) - 1); // Definition of W for adaptive path-length
				}
				if ( ( iterCountAPP > 5 ) && ( iterCountAPP <= 10 ) && (setOuterAPPTuning==0) ) {
					W[2*i*numberOfGenerators+j] = 0.0; // Definition of W for adaptive path-length
					W[(2*i+1)*numberOfGenerators+j] = 0.0; // Definition of W for adaptive path-length
					alphaAPP[2*i*numberOfGenerators+j] = 5.0;
					alphaAPP[(2*i+1)*numberOfGenerators+j] = 5.0;
				}
				if ( ( iterCountAPP > 10 ) && ( iterCountAPP <= 15 ) && (setOuterAPPTuning==0) ) {
					W[2*i*numberOfGenerators+j] = 0.0; // Definition of W for adaptive path-length
					W[(2*i+1)*numberOfGenerators+j] = 0.0; // Definition of W for adaptive path-length
					alphaAPP[2*i*numberOfGenerators+j] = 2.5;
					alphaAPP[(2*i+1)*numberOfGenerators+j] = 2.5;
				}
				if ( ( iterCountAPP > 15 ) && ( iterCountAPP <= 20 ) && (setOuterAPPTuning==0) ) {
					W[2*i*numberOfGenerators+j] = 0.0; // Definition of W for adaptive path-length
					W[(2*i+1)*numberOfGenerators+j] = 0.0; // Definition of W for adaptive path-length
					alphaAPP[2*i*numberOfGenerators+j] = 1.25;
					alphaAPP[(2*i+1)*numberOfGenerators+j] = 1.25;
				}
				if ( ( iterCountAPP > 20 ) && (setOuterAPPTuning==0) ) {
					W[2*i*numberOfGenerators+j] = 0.0; // Definition of W for fixed Rho
					W[(2*i+1)*numberOfGenerators+j] = 0.0; // Definition of W for fixed Rho
					alphaAPP[2*i*numberOfGenerators+j] = 0.5;
					alphaAPP[(2*i+1)*numberOfGenerators+j] = 0.5;
				}
				// Calculation of Adaptive Step-length
				controllerSum[2*i*numberOfGenerators+j] = controllerSum[2*i*numberOfGenerators+j] + W[2*i*numberOfGenerators+j];
				alphaAPP1[2*i*numberOfGenerators+j] = alphaAPP[2*i*numberOfGenerators+j]; // Store previous alphaAPP
				alphaAPP[2*i*numberOfGenerators+j] = ( alphaAPP1[2*i*numberOfGenerators+j] ) * ( exp( ( lambdaAdap * W[2*i*numberOfGenerators+j] ) + ( muAdap * ( W[2*i*numberOfGenerators+j] - Wprev[2*i*numberOfGenerators+j] ) ) + ( xiAdap * controllerSum[2*i*numberOfGenerators+j]  ) ) ); // Next iterate value of Rho
				Wprev[2*i*numberOfGenerators+j] = W[2*i*numberOfGenerators+j]; // Buffering
				controllerSum[(2*i+1)*numberOfGenerators+j] = controllerSum[(2*i+1)*numberOfGenerators+j] + W[(2*i+1)*numberOfGenerators+j];
				alphaAPP1[(2*i+1)*numberOfGenerators+j] = alphaAPP[(2*i+1)*numberOfGenerators+j]; // Store previous alphaAPP
				alphaAPP[(2*i+1)*numberOfGenerators+j] = ( alphaAPP1[(2*i+1)*numberOfGenerators+j] ) * ( exp( ( lambdaAdap * W[(2*i+1)*numberOfGenerators+j] ) + ( muAdap * ( W[(2*i+1)*numberOfGenerators+j] - Wprev[(2*i+1)*numberOfGenerators+j] ) ) + ( xiAdap * controllerSum[(2*i+1)*numberOfGenerators+j]  ) ) ); // Next iterate value of Rho
				Wprev[(2*i+1)*numberOfGenerators+j] = W[(2*i+1)*numberOfGenerators+j]; // Buffering
			}
		}
		// Update power disagreement Lagrange Multipliers
		for ( int i = 0; i < ( dispatchIntervals - 1 ); ++i ) {
			for ( int j = 0; j < numberOfGenerators; ++j ) {
				lambdaAPP[2*i*numberOfGenerators+j] = lambdaAPP[2*i*numberOfGenerators+j] + alphaAPP[2*i*numberOfGenerators+j] * (powDiff[2*i*numberOfGenerators+j]); // what I think about myself Vs. what next door fellow thinks about me
				lambdaAPP[(2*i+1)*numberOfGenerators+j] = lambdaAPP[(2*i+1)*numberOfGenerators+j] + alphaAPP[(2*i+1)*numberOfGenerators+j] * (powDiff[(2*i+1)*numberOfGenerators+j]); // what I think about next door fellow Vs. what next door fellow thinks about himself
			}
		}
		//++iterCountAPP; // increment the APP iteration counter
		double tolAPP = 0.0;
		double tolAPPDelayed = 0.0; // APP tolerance, excluding the first (dummy) interval
		matrixResultAPPOut << (iterCountAPP-1) << "\t";
		for ( int i = 0; i < consLagDim; ++i ) {
			tolAPP = tolAPP + pow(powDiff[i], 2);
			if (i>=2*numberOfGenerators)
				tolAPPDelayed = tolAPPDelayed + pow(powDiff[i], 2);
			matrixResultAPPOut << powDiff[i] << "\t";
		}
		matrixResultAPPOut << "\n";
		finTol = sqrt(tolAPP);
		finTolDelayed = sqrt(tolAPPDelayed);
		matrixResultAPPOut << (iterCountAPP-1) << "\t" << finTol << "\t" << finTolDelayed <<"\n";
		++iterCountAPP;
		cout << "\nFinal Value of Outer APP Tolerance " << finTol << "\nAnd Final Value of Outer APP Delayed Tolerance " << finTolDelayed << endl;
	} while (finTolDelayed>=0.7); //Check the termination criterion of the APP iterations
	//}
//****************************************END OF AUXILIARY PROBLEM PRINCIPLE (APP) COARSE GRAINED DECOMPOSITION COMPONENT******************************************************//

	cout << "\n*** LASCOPF FOR LOAD (DEMAND) VARIATION TRACKING SIMULATION SUPERNETWORK LAYER ENDS ***\n" << endl;
	clock_t stop_s = clock();  // end
	matrixResultAPPOut << "\nExecution Outermost layer time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;
	matrixResultAPPOut << "\nVirtual Outermost layer Execution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC  - actualSuperNetTime + accumulate(largestSuperNetTimeVec.begin(), largestSuperNetTimeVec.end(), 0.0)<< endl;
	cout << "\nVirtual Outermost layer Execution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC  - actualSuperNetTime + accumulate(largestSuperNetTimeVec.begin(), largestSuperNetTimeVec.end(), 0.0);
	cout << "\nExecution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;
	delete environmentGUROBI; // Free the memory of the GUROBI environment object
	for ( int i = 0; i < dispatchIntervals; ++i ) {
		delete futureNetVector[i]; // Free the memory of the vector of future network instances
	}

	return 0; // indicates successful program termination

} // end main

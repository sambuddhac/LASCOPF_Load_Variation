// Member functions for class Generator.
#include <iostream>
#include <iomanip>
// include Generator class definition from generator.h
#include "generator.h"
// Include definition of Node class 
#include "node.h"
// Include Generator solver class defintion
#include "gensolverFirstBase.h" // definition of Gensolver class for base case scenario
#include "gensolverIntermediateBase.h" // definition of Gensolver class for base case scenario
#include "gensolverLastBase.h" // definition of Gensolver class for base case scenario
#include "gensolverCont.h" // definition of Gensolver class for contingency scenario
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file

using namespace std;

Generator::Generator( int idOfGen, int interval, int lastFlag, int contScenarioCount, int baseCont, Node *nodeConng, GensolverFirstBase &paramOfGenFirstBase, GensolverInterBase &paramOfGenInterBase, GensolverLastBase &paramOfGenLastBase, GensolverCont &paramOfGenCont, int countOfContingency, int genTotal ) // constructor begins
	: genID( idOfGen ),
	  numberOfGenerators(genTotal),
	  dispatchInterval(interval),
	  flagLast(lastFlag),
	  scenarioContCount( contScenarioCount ),
	  baseContScenario( baseCont ),
	  connNodegPtr( nodeConng ),
	  genSolverFirstBase( paramOfGenFirstBase ),
	  genSolverInterBase( paramOfGenInterBase ),
	  genSolverLastBase( paramOfGenLastBase ),
	  genSolverCont( paramOfGenCont ),
	  contCountGen( countOfContingency )
{
	//cout << "\nInitializing the parameters of the generator with ID: " << genID << endl;
	connNodegPtr->setgConn(idOfGen); // increments the generation connection variable to node
	PgenPrev=genSolverFirstBase.getPgPrev();
	setGenData(); // calls setGenData member function to set the parameter values

} // constructor ends

Generator::~Generator() // destructor
{
	//cout << "\nThe generator object having ID " << genID << " have been destroyed.\n";

} // end of destructor

int Generator::getGenID() // function getGenID begins
{
	return genID; // returns the ID of the generator object
} // end of getGenID function

int Generator::getGenNodeID() // function getGenNodeID begins
{
	return connNodegPtr->getNodeID(); // returns the ID number of the node to which the generator object is connected
} // end of getGenNodeID function

void Generator::setGenData() // start setGenData function
{
	Pg = 0.0; // Initialize power iterate
	Thetag = 0.0; // Initialize angle iterate
	v = 0.0; // Initialize the Lagrange multiplier corresponding voltage angle constraint to zero
	
} // end of setGenData function

void Generator::gpowerangleMessage(int outerAPPIt, int  APPItCount, double gsRho, double Pgenprev, double Pgenavg, double Powerprice, double Angpriceavg, double Angavg, double Angprice, double PgenPrevAPP, double PgenAPP, double PgenAPPInner, double PgenNextAPP, double AAPPExternal, double BAPPExternal, double DAPPExternal, double LambAPP1External, double LambAPP2External, double LambAPP3External, double LambAPP4External, double BAPP[], double LambAPP1[]) //const // function gpowerangleMessage begins
{
	double BAPPNew[contCountGen];
	double LambdaAPPNew[contCountGen];
	for (int counterCont = 0; counterCont < contCountGen; ++counterCont) {
		BAPPNew[counterCont]=0; 
		LambdaAPPNew[counterCont]=0;
	}
	if ( baseContScenario == 0 ) { // Use the solver for first dispatch interval
		if ((dispatchInterval==0) && (flagLast==0)) {
			for (int counterCont = 0; counterCont < contCountGen; ++counterCont) {
				BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
				LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
			}
			genSolverFirstBase.mainsolve( outerAPPIt, APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPP, PgenAPPInner, PgenNextAPP, BAPPExternal, DAPPExternal, LambAPP1External, LambAPP2External, BAPPNew, LambdaAPPNew ); // calls the Generator optimization solver
			Pg = genSolverFirstBase.getPSol(); // get the Generator Power iterate
			PgenNext = genSolverFirstBase.getPNextSol();
			PgenPrev = genSolverFirstBase.getPgPrev();
			Thetag = *(genSolverFirstBase.getThetaPtr());
			//Thetag = genSolver.getThetaPtr(); // get the Generator voltage angle iterate
		}
		if ((dispatchInterval!=0) && (flagLast==0)) {
			for (int counterCont = 0; counterCont < contCountGen; ++counterCont) {
				BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
				LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
			}
			genSolverInterBase.mainsolve( outerAPPIt, APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPP, PgenAPPInner, PgenNextAPP, PgenPrevAPP, AAPPExternal, BAPPExternal, DAPPExternal, LambAPP1External, LambAPP2External, LambAPP3External, LambAPP4External, BAPPNew, LambdaAPPNew ); // calls the Generator optimization solver
			Pg = genSolverInterBase.getPSol(); // get the Generator Power iterate
			PgenNext = genSolverInterBase.getPNextSol();
			PgenPrev = genSolverInterBase.getPPrevSol();
			Thetag = *(genSolverInterBase.getThetaPtr());
			//Thetag = genSolver.getThetaPtr(); // get the Generator voltage angle iterate
		}
		if ((dispatchInterval!=0) && (flagLast==1)) {
			for (int counterCont = 0; counterCont < contCountGen; ++counterCont) {
				BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
				LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
			}
			genSolverLastBase.mainsolve( outerAPPIt, APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPP, PgenAPPInner, PgenPrevAPP, AAPPExternal, BAPPExternal, LambAPP3External, LambAPP4External,  BAPPNew, LambdaAPPNew ); // calls the Generator optimization solver
			Pg = genSolverLastBase.getPSol(); // get the Generator Power iterate
			PgenNext = genSolverLastBase.getPNextSol();
			PgenPrev = genSolverLastBase.getPPrevSol();
			Thetag = *(genSolverLastBase.getThetaPtr());
			//Thetag = genSolver.getThetaPtr(); // get the Generator voltage angle iterate
		}
	}
	else { // Use the solver for first dispatch interval
		genSolverCont.mainsolve( APPItCount, gsRho, Pgenprev, Pgenavg, Powerprice, Angpriceavg, Angavg, Angprice, PgenAPPInner, -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)], LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)] ); // calls the Generator optimization solver
		Pg = genSolverCont.getPSol(); // get the Generator Power iterate
		//Thetag = genSolver.getThetaPtr(); // get the Generator voltage angle iterate
		Thetag = *(genSolverCont.getThetaPtr());
		//*cout << "\nThiterate from generator: " << *(Thetag+i) << endl;
	}

	connNodegPtr->powerangleMessage( Pg, v, Thetag ); // passes to node object the corresponding iterates of power, angle, v, and number of scenarios
} // function gpowerangleMessage ends

void Generator::gpowerangleMessageGUROBI(int outerAPPIt, int  APPItCount, double gsRho, double Pprevit, double Pnetavg, double uprev, double vprevavg, double Aprevavg, double vprev, double PgenPrevAPP, double PgenAPP, double PgenAPPInner, double PgenNextAPP, double AAPPExternal, double BAPPExternal, double DAPPExternal, double LambAPP1External, double LambAPP2External, double LambAPP3External, double LambAPP4External, double BAPP[], double LambAPP1[], GRBEnv* environmentGUROBI) //const // function gpowerangleMessage begins
{
	// CREATION OF THE MIP SOLVER INSTANCE //
        int dimRow; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        int dimCol; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
	if ( ( baseContScenario == 0 ) && (dispatchInterval==0) ){ 
        	dimRow = 6; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        	dimCol = 3; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
	}
	if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==0) ){ 
        	dimRow = 6; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        	dimCol = 4; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
	}
	if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==1) ){ 
        	dimRow = 6; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        	dimCol = 3; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
	}
	if (( baseContScenario != 0 )){ 
        	dimRow = 2; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        	dimCol = 2; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
	}
	// Instantiate GUROBI Problem model
	GRBModel *modelGenQP = new GRBModel(*environmentGUROBI);
    	modelGenQP->set(GRB_StringAttr_ModelName, "assignment");
	modelGenQP->set(GRB_IntParam_OutputFlag, 0);
	GRBVar decvar[dimCol+1];
	double z; // variable to store the objective value

	// SPECIFICATION OF PROBLEM PARAMETERS //
	// Dummy Decision Variable //
	decvar[0] = modelGenQP->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
	//Decision Variable Definitions, Bounds, and Objective Function Co-efficients//
	int colCount = 1;
	if ( ( baseContScenario == 0 ) && (dispatchInterval==0) ){ 
		//Columns corresponding to Power Generation continuous variables for different generators//
		decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		++colCount;
		//Columns corresponding to Power Generation continuous variables for different generators for next interval//
		decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		++colCount;
		//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
		decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
	}
	if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==0) ){ 
		//Columns corresponding to Power Generation continuous variables for different generators//
		decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		++colCount;
		//Columns corresponding to Power Generation continuous variables for different generators for next interval//
		decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		++colCount;
		//Columns corresponding to Power Generation continuous variables for different generators for previous interval//
		decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		++colCount;
		//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
		decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
	}
	if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==1) ){ 
		//Columns corresponding to Power Generation continuous variables for different generators//
		decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		++colCount;
		//Columns corresponding to Power Generation continuous variables for different generators for previous interval//
		decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		++colCount;
		//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
		decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
	}
	if (( baseContScenario != 0 )){ 
		//Columns corresponding to Power Generation continuous variables for different generators//
		decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		++colCount;
		//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
		decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
	}

	//Setting Objective//
	GRBQuadExpr obj = 0.0;
	// Objective Contribution from Dummy Decision Variable //
	obj += 0*(decvar[0]);
	colCount = 1;
	double BAPPNew[contCountGen];
	double LambdaAPPNew[contCountGen];
	double BAPPSum = 0;
	double LambdaAPPSum = 0;
	for (int counterCont = 0; counterCont < contCountGen; ++counterCont) {
		BAPPNew[counterCont]=0; 
		LambdaAPPNew[counterCont]=0;
	}
	if ( baseContScenario == 0 ) { // Use the solver for first dispatch interval
		for (int counterCont = 0; counterCont < contCountGen; ++counterCont) {
			BAPPNew[counterCont]=BAPP[counterCont*numberOfGenerators+(genID-1)]; 
			LambdaAPPNew[counterCont]=LambAPP1[counterCont*numberOfGenerators+(genID-1)];
			BAPPSum += BAPPNew[counterCont];
			LambdaAPPSum += LambdaAPPNew[counterCont];
		}
	}
	else { // Use the solver for first dispatch interval
		BAPPSum = -BAPP[(scenarioContCount-1)*numberOfGenerators+(genID-1)]; 
		LambdaAPPSum = LambAPP1[(scenarioContCount-1)*numberOfGenerators+(genID-1)];
	}
	//Columns corresponding to Power Generation continuous variables for different generators//
	if ( ( baseContScenario == 0 ) && (dispatchInterval==0) ){ 
		obj += (genSolverFirstBase.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverFirstBase.getLinCoeff())*(decvar[colCount])+(genSolverFirstBase.getConstCoeff())+((genSolverFirstBase.getBeta())/2)*((decvar[colCount])-PgenAPP)*((decvar[colCount])-PgenAPP)+(genSolverFirstBase.getIntGamma())*((decvar[colCount])*BAPPSum)+(decvar[colCount])*LambdaAPPSum+(genSolverFirstBase.getGamma())*((decvar[colCount])*BAPPExternal)+(decvar[colCount])*LambAPP1External+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
		++colCount;
		obj += ((genSolverFirstBase.getBeta())/2)*((decvar[colCount])-PgenNextAPP)*((decvar[colCount])-PgenNextAPP)+(genSolverFirstBase.getGamma())*((decvar[colCount])*DAPPExternal)+(decvar[colCount])*LambAPP2External;
	}
	if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==0) ){ 
		obj += (genSolverInterBase.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverInterBase.getLinCoeff())*(decvar[colCount])+(genSolverInterBase.getConstCoeff())+((genSolverInterBase.getBeta())/2)*((decvar[colCount])-PgenAPP)*((decvar[colCount])-PgenAPP)+(genSolverInterBase.getIntGamma())*((decvar[colCount])*BAPPSum)+(decvar[colCount])*LambdaAPPSum+(genSolverInterBase.getGamma())*((decvar[colCount])*BAPPExternal)+(decvar[colCount])*(LambAPP1External-LambAPP4External)+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
		++colCount;
		obj += ((genSolverInterBase.getBeta())/2)*((decvar[colCount])-PgenNextAPP)*((decvar[colCount])-PgenNextAPP)+(genSolverInterBase.getGamma())*((decvar[colCount])*DAPPExternal)+(decvar[colCount])*LambAPP2External;
		++colCount;
		obj += ((genSolverInterBase.getBeta())/2)*((decvar[colCount])-PgenPrevAPP)*((decvar[colCount])-PgenPrevAPP)+(genSolverInterBase.getGamma())*((decvar[colCount])*AAPPExternal)-(decvar[colCount])*LambAPP3External;
	}
	if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==1) ){ 
		obj += (genSolverLastBase.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverLastBase.getLinCoeff())*(decvar[colCount])+(genSolverLastBase.getConstCoeff())+((genSolverLastBase.getBeta())/2)*((decvar[colCount])-PgenAPP)*((decvar[colCount])-PgenAPP)+(genSolverLastBase.getIntGamma())*((decvar[colCount])*BAPPSum)+(decvar[colCount])*LambdaAPPSum+(genSolverLastBase.getGamma())*((decvar[colCount])*BAPPExternal)+(decvar[colCount])*(-LambAPP4External)+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
		++colCount;
		obj += ((genSolverLastBase.getBeta())/2)*((decvar[colCount])-PgenPrevAPP)*((decvar[colCount])-PgenPrevAPP)+(genSolverLastBase.getGamma())*((decvar[colCount])*AAPPExternal)-(decvar[colCount])*LambAPP3External;
	}
	if ( baseContScenario != 0 )
		obj += (genSolverCont.getQuadCoeff())*(decvar[colCount])*(decvar[colCount])+(genSolverCont.getLinCoeff())*(decvar[colCount])+(genSolverCont.getConstCoeff())+((genSolverCont.getBeta())/2)*((decvar[colCount])-PgenAPP)*((decvar[colCount])-PgenAPP)+(genSolverCont.getGamma())*((decvar[colCount])*BAPPSum)-(decvar[colCount])*LambdaAPPSum+(gsRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
	++colCount;
	//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
	obj += (gsRho/2)*(decvar[colCount]-vprevavg-Aprevavg+vprev )*(decvar[colCount]-vprevavg-Aprevavg+vprev );

	modelGenQP->setObjective(obj, GRB_MINIMIZE);
	//Row Definitions: Specification of b<=Ax<=b//
	GRBLinExpr lhs[dimRow+1];
	//Row Definitions and Bounds Corresponding to Constraints/
	//Non-Zero entries of A matrix (Constraint/Coefficient matrix entries)//
	// Dummy Constraint //
	lhs[0] = 0*(decvar[0]);
	modelGenQP->addConstr(lhs[0], GRB_EQUAL, 0);
	int rCount = 1; // Initialize the row count
	if ( ( baseContScenario == 0 ) && (dispatchInterval==0) ){ 
		// Coefficients corresponding to lower generation limits
		lhs[rCount] = 0;
		lhs[rCount] += decvar[rCount];
		modelGenQP->addConstr(lhs[rCount] >= (genSolverFirstBase.getPMin()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to upper generation limits
		lhs[rCount] = 0;
		lhs[rCount] += decvar[rCount - 1];
		modelGenQP->addConstr(lhs[rCount] <= (genSolverFirstBase.getPMax()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to lower ramp limits for next interval
		lhs[rCount] = 0;
		lhs[rCount] += (decvar[rCount-1]-decvar[rCount-2]);
		modelGenQP->addConstr(lhs[rCount] >= (genSolverFirstBase.getRMin()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to upper ramp limits for next interval
		lhs[rCount] = 0;
		lhs[rCount] += (decvar[rCount-2]-decvar[rCount-3]);
		modelGenQP->addConstr(lhs[rCount] <= (genSolverFirstBase.getRMax()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to lower ramp limits for previous interval
		lhs[rCount] = 0;
		lhs[rCount] += (decvar[rCount-4]-genSolverFirstBase.getPgPrev());
		modelGenQP->addConstr(lhs[rCount] >= (genSolverFirstBase.getRMin()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to upper ramp limits for previous interval
		lhs[rCount] = 0;
		lhs[rCount] += (decvar[rCount-5]-genSolverFirstBase.getPgPrev());
		modelGenQP->addConstr(lhs[rCount] <= (genSolverFirstBase.getRMax()));
	}
	if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==0) ){ 
		// Coefficients corresponding to lower generation limits
		lhs[rCount] = 0;
		lhs[rCount] += decvar[rCount];
		modelGenQP->addConstr(lhs[rCount] >= (genSolverInterBase.getPMin()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to upper generation limits
		lhs[rCount] = 0;
		lhs[rCount] += decvar[rCount - 1];
		modelGenQP->addConstr(lhs[rCount] <= (genSolverInterBase.getPMax()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to lower ramp limits for next interval
		lhs[rCount] = 0;
		lhs[rCount] += (decvar[rCount-1]-decvar[rCount-2]);
		modelGenQP->addConstr(lhs[rCount] >= (genSolverInterBase.getRMin()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to upper ramp limits for next interval
		lhs[rCount] = 0;
		lhs[rCount] += (decvar[rCount-2]-decvar[rCount-3]);
		modelGenQP->addConstr(lhs[rCount] <= (genSolverInterBase.getRMax()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to lower ramp limits for previous interval
		lhs[rCount] = 0;
		lhs[rCount] += (decvar[rCount-4]-decvar[rCount-2]);
		modelGenQP->addConstr(lhs[rCount] >= (genSolverInterBase.getRMin()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to upper ramp limits for previous interval
		lhs[rCount] = 0;
		lhs[rCount] += (decvar[rCount-5]-decvar[rCount-3]);
		modelGenQP->addConstr(lhs[rCount] <= (genSolverInterBase.getRMax()));
	}
	if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==1) ){ 
		// Coefficients corresponding to lower generation limits
		lhs[rCount] = 0;
		lhs[rCount] += decvar[rCount];
		modelGenQP->addConstr(lhs[rCount] >= (genSolverLastBase.getPMin()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to upper generation limits
		lhs[rCount] = 0;
		lhs[rCount] += decvar[rCount - 1];
		modelGenQP->addConstr(lhs[rCount] <= (genSolverLastBase.getPMax()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to lower ramp limits for next interval
		lhs[rCount] = 0;
		lhs[rCount] += (PgenAPP-decvar[rCount-2]);
		modelGenQP->addConstr(lhs[rCount] >= (genSolverLastBase.getRMin()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to upper ramp limits for next interval
		lhs[rCount] = 0;
		lhs[rCount] += (PgenAPP-decvar[rCount-3]);
		modelGenQP->addConstr(lhs[rCount] <= (genSolverLastBase.getRMax()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to lower ramp limits for previous interval
		lhs[rCount] = 0;
		lhs[rCount] += (decvar[rCount-4]-decvar[rCount-3]);
		modelGenQP->addConstr(lhs[rCount] >= (genSolverLastBase.getRMin()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to upper ramp limits for previous interval
		lhs[rCount] = 0;
		lhs[rCount] += (decvar[rCount-5]-decvar[rCount-4]);
		modelGenQP->addConstr(lhs[rCount] <= (genSolverLastBase.getRMax()));
	}
	if ( baseContScenario != 0 ){ 
		// Coefficients corresponding to lower generation limits
		lhs[rCount] = 0;
		lhs[rCount] += decvar[rCount];
		modelGenQP->addConstr(lhs[rCount] >= (genSolverFirstBase.getPMin()));
		++rCount; // Increment the row count to point to the next generator object
		// Coefficients corresponding to upper generation limits
		lhs[rCount] = 0;
		lhs[rCount] += decvar[rCount - 1];
		modelGenQP->addConstr(lhs[rCount] <= (genSolverFirstBase.getPMax()));
	}
	// RUN THE OPTIMIZATION SIMULATION ALGORITHM //
	modelGenQP->optimize(); // Solves the optimization problem
	int stat = modelGenQP->get(GRB_IntAttr_Status); // Outputs the solution status of the problem 

	// DISPLAY THE SOLUTION DETAILS //
	if (stat == GRB_INFEASIBLE){
		cout << "\nThe solution to the problem is INFEASIBLE." << endl;
		delete modelGenQP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_INF_OR_UNBD) {
		cout << "\nNO FEASIBLE or BOUNDED solution to the problem exists." << endl;
		delete modelGenQP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_UNBOUNDED) {
		cout << "\nThe solution to the problem is UNBOUNDED." << endl;
		delete modelGenQP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_OPTIMAL) {
		//Get the Optimal Objective Value results//
		z = modelGenQP->get(GRB_DoubleAttr_ObjVal);
		// writing results of different variables
		vector<double> x; // Vector for storing decision variable output 
		x.push_back(0); // Initialize the decision Variable vector
		objOpt = 0;
		//Power Generation
		int arrayInd = 1;
		if ( ( baseContScenario == 0 ) && (dispatchInterval==0) ){ 
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
			objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
			++arrayInd;
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			PgenNext = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for next interval	
			// Internal node voltage phase angle variables
			++arrayInd;
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate	
		}
		if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==0) ){ 
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
			objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
			++arrayInd;
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			PgenNext = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for next interval	
			// Internal node voltage phase angle variables
			++arrayInd;
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			PgenPrev = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for previous interval	
			// Internal node voltage phase angle variables
			++arrayInd;
			// Internal node voltage phase angle variables
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate	
		}
		if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==1) ){ 
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
			objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
			++arrayInd;
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			PgenPrev = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate for previous interval	
			// Internal node voltage phase angle variables
			++arrayInd;
			// Internal node voltage phase angle variables
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate	
		}
		if ( baseContScenario != 0 ){ 
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
			objOpt += (genSolverFirstBase.getQuadCoeff())*(Pg)*(Pg)+(genSolverFirstBase.getLinCoeff())*(Pg)+(genSolverFirstBase.getConstCoeff());
			++arrayInd;
			// Internal node voltage phase angle variables
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate	
		}
		connNodegPtr->powerangleMessage( Pg, v, Thetag ); // passes to node object the corresponding iterates of power, angle and v
		delete modelGenQP; 
	}
} // function gpowerangleMessage ends

double Generator::genPower() //const // function genPower begins
{
	return Pg; // returns the Pg iterate
} // function genPower ends

double Generator::genPowerPrev() //const // function genPower begins
{
	if (dispatchInterval==0)
		return genSolverFirstBase.getPgPrev(); // returns the Pg iterate
	else 
		return PgenPrev;
} // function genPower ends

double Generator::genPowerNext() //const // function genPower begins
{
	if (flagLast==1)
		return Pg; // returns the Pg iterate
	else 
		return PgenNext;
} // function genPower ends

double Generator::objectiveGen() // function objectiveGen begins
{
	if ( ( baseContScenario == 0 ) && (dispatchInterval==0) ) {
		return genSolverFirstBase.getObj(); //returns the evaluated objective
	}
	if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==0) ) {
		return genSolverInterBase.getObj(); //returns the evaluated objective
	}
	if ( ( baseContScenario == 0 ) && (dispatchInterval!=0) && (flagLast==1) ) {
		return genSolverLastBase.getObj(); //returns the evaluated objective
	}
	if ( baseContScenario != 0 ) {
		return genSolverCont.getObj(); //returns the evaluated objective
	}
} // function objectiveGen ends

double Generator::objectiveGenGUROBI() // Objective from GUROBI ADMM 
{
	return objOpt; //returns the evaluated objective
} // function objectiveGen ends

double Generator::calcPtilde() //const // function calcPtilde begins
{
	double P_avg = connNodegPtr->PavMessage(); // Gets average power from the corresponding node object
	double Ptilde = Pg - P_avg; // calculates the difference between power iterate and average
	return Ptilde; // returns the difference
} // function calcPtilde ends

double Generator::calcPavInit() const // function calcPavInit begins
{
	return connNodegPtr->devpinitMessage(); // seeks the initial Ptilde from the node
} // function calcPavInit ends

double Generator::getu() const // function getu begins
{
	double u = connNodegPtr->uMessage(); // gets the value of the price corresponding to power balance from node
	//cout << "u: " << u << endl;
	return u; // returns the price
} // function getu ends

double Generator::calcThetatilde() //const // function calcThetatilde begins
{
	//cout << "Thetag: " << Thetag << endl;
	double Theta_avg = connNodegPtr->ThetaavMessage(); // get the average voltage angle at the particular node
	//cout << "Theta_avg: " << Theta_avg << endl;
	double Theta_tilde = Thetag - Theta_avg; // claculate the deviation between the voltage angle of the device and the average
	return Theta_tilde; // return the deviation
} // function calcThetatilde ends

double Generator::calcvtilde() const // function calcvtilde begins
{
	double v_avg = connNodegPtr->vavMessage(); // get the average of the Lagrange multiplier corresponding to voltage angle balance
	//cout << "v_avg: " << v_avg << endl;
	double v_tilde = v - v_avg; // calculate the deviation of the node Lagrange multiplier to the average
	return v_tilde; // return the deviation
} // function calcvtilde ends

double Generator::getv() // function getv begins
{
	//cout << "v_initial: " << v << endl;
	v = v + calcThetatilde(); // Calculate the value of the Lagrange multiplier corresponding to angle constraint
	//cout << "v_final: " << v << endl;
	return v; // Calculate the value of the Lagrange multiplier corresponding to angle constraint
} // function getv ends		
double Generator::getPMax(){return genSolverFirstBase.getPMax();}
double Generator::getPMin(){return genSolverFirstBase.getPMin();}
double Generator::getQuadCoeff(){return genSolverFirstBase.getQuadCoeff();}
double Generator::getLinCoeff(){return genSolverFirstBase.getLinCoeff();}
double Generator::getConstCoeff(){return genSolverFirstBase.getConstCoeff();}
double Generator::getPgenPrev(){return genSolverFirstBase.getPgPrev();}
double Generator::getPgenNext(){return genSolverFirstBase.getPNextSol();}
double Generator::getRMax(){return genSolverFirstBase.getRMax();}
double Generator::getRMin(){return genSolverFirstBase.getRMin();}

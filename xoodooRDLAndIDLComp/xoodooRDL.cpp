#include<iostream>
#include<vector>
#include<cmath>
#include<omp.h>
#include"xoodoo.h"
#include"differential.h"
#include"linear.h"
#include"differential-linear.h"
#include"miniZincReader.h"
#include"internalDifferential.h"
#include"internalLinear.h"
#include"internalDifferential-linear.h"

using namespace std;

int main() {
	int rot = 1;
	XoodooBase xoodooBase;
	MiniZincReader reader(xoodooBase);
	RotationalDifferentialBase rdBase(xoodooBase, rot);
	InternalDifferentialBase idBase(xoodooBase, rdBase);
	LinearBase linBase(xoodooBase, rdBase);
	InternalLinearBase iLinBase(xoodooBase, idBase, linBase);
	RotationalDifferentialLinearBase rdlBase(xoodooBase, linBase, rdBase, 2);
	InternalDifferentialLinearBase idlBase(xoodooBase, iLinBase, idBase, rdlBase);



	// enter the input internal difference of EM here
	vector<int> iDiff13DArr = { 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	// enter the linear mask on the output internal difference of EM here
	vector<int> ilin03DArr = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 };

	// enter the input rotational difference of EM here
	vector<int> diff13DArr = { 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	// enter the linear mask on the output rotational difference of EM here
	vector<int> linear43DArr = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	StateTypeBySheet stateDiff1 = reader.ReadStateFrom3DArray(diff13DArr);
	StateTypeBySheet stateLin4 = reader.ReadStateFrom3DArray(linear43DArr);

	HalfStateTypeBySheet halfStateDiff1 = reader.ReadHalfStateFrom3DArray(iDiff13DArr);

	HalfStateTypeBySheet halfStateLin0 = reader.ReadHalfStateFrom3DArray(ilin03DArr);

	// double idProb = idBase.ComputeDifferentialProbabilityManually(halfStateDiff0, halfStateDiff1, 0, 1, 6);
	// cout << "ID prob " << idProb << " log2 " << log2(idProb) << endl;

	// double experimentalIDProb = idBase.ExperimentalComputation(halfStateDiff0, halfStateDiff1, 0, 1, 6);
	// cout << experimentalIDProb << " log2 " << log2(experimentalIDProb) << endl;

	// double ilinCor = iLinBase.ComputeApproximationCorrelationSquareAutomatically(halfStateLin0, halfStateLin1, 5, 5, 6);
	// cout << "iLin cor " << ilinCor << " log2 " << log2(abs(ilinCor)) << endl;




	int startr = 2; // the starting number of rounds for EM
	int endr = 5; // the ending number of rounds for EM

	double finalCor = 0.0;
	int linRot = 0;



	// finalCor = rdlBase.ComputeThreeRoundRDLCor(stateDiff1, stateLin4, startr, endr, 6); // if you are considering RDL analysis, call this to compute the RDL correlation of the 3-round EM; otherwise, comment it out; if the program stuck here, try to adjust the order of forward and backward division calls in the function ComputeThreeRoundRDLCor in differential-linear.h
	finalCor = idlBase.ComputeThreeRoundIDLCor(halfStateDiff1, halfStateLin0, startr, endr, 6); // if you are considering IDL analysis, call this to compute the IDL correlation of the 3-round EM; otherwise, comment it out; if the program stuck here, try to adjust the order of forward and backward division calls in the function ComputeThreeRoundIDLCor in internalDifferential-linear.h

	// finalCor = rdlBase.ComputeFourRoundRDLCor(stateDiff1, stateLin4, startr, endr, 6); // if you would like to compute the RDL correlation of the 4-round EM, call this to compute the RDL correlation of the 4-round EM (remember to change the startr and endr correspondingly); otherwise, comment it out; if the program stuck here, try to adjust the order of forward and backward division calls in the function ComputeFourRoundRDLCor in differential-linear.h

	// finalCor = idlBase.ComputeFourRoundIDLCor(halfStateDiff1, halfStateLin0, startr, endr, 6); // if you would like to compute the IDL correlation of the 4-round EM, call this to compute the IDL correlation of the 4-round EM (remember to change the startr and endr correspondingly); otherwise, comment it out; if the program stuck here, try to adjust the order of forward and backward division calls in the function ComputeFourRoundIDLCor in internalDifferential-linear.h

	if (finalCor != 0.0)
	{
		cout << "final cor log2 " << log2(abs(finalCor))  << endl;
	}





	// experimental computation
	double expCor = 0.0;



	int threadsNum = omp_get_max_threads();
	cout << "threads num " << threadsNum << endl;

#pragma omp parallel for num_threads(threadsNum)
	for (int i = 0; i < threadsNum; i++)
	{
		double tmpCor = idlBase.ExperimentalComputation(halfStateDiff1, halfStateLin0, startr, endr, 6); // If you are considering IDL analysis, call this to evaluate the RDL correlation of EM experimentally; otherwise, comment it out
		// double tmpCor = rdlBase.ExperimentalComputation(stateDiff1, stateLin4, startr, endr, 6); // If you are considering RDL analysis, call this to evaluate the RDL correlation of EM experimentally; otherwise, comment it out
		cout << omp_get_thread_num() << " " << i << " " << tmpCor << endl;
#pragma omp critical 
		{
			expCor += tmpCor;
		}
	}
	cout << "Experimental cor log2 " << log2(abs(expCor / threadsNum)) << endl;
}




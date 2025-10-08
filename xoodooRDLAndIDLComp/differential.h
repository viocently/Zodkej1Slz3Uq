#pragma once
#include<iostream>
#include<array>
#include<vector>
#include<random>
#include"xoodoo.h"
#include"gurobi_c++.h"


using namespace std;

using RDLaneMILPType = array<GRBVar, 32>;
using RDSheetMILPType = array<RDLaneMILPType, 3>;
using RDPlaneMILPType = array<RDLaneMILPType, 4>;
using RDStateMILPTypeBySheet = array<RDSheetMILPType, 4>;

RDLaneMILPType RDRotateLeft(RDLaneMILPType laneVars, int shift) {
	RDLaneMILPType newLaneVars;
	for (int i = 0; i < 32; i++)
		newLaneVars[i] = laneVars[(i - shift + 32) % 32];
	return newLaneVars;
}


class RotationalDifferentialBase
{
private:
	const XoodooBase& xoodooBase;

	int rot;
	int milpThreads;
	array<LaneType, 12> rotationalDifferencePerRoundConstants;
	array<array<int, 8>, 8> ddt;
	array<array<double, 8>, 8> ddtProb;
	array<array<double, 8>, 8> ddtProbLog2;
	vector<vector<int>> outputDifferencesPerInputDifferenceForChi;

	vector<array<int, 7>> milpChiInequalities;

	bool SwitchToNextPossibleOutputDifferences(vector<int>& activeColumnOutputDiffIndexes, const vector<int>& activeColumnInputDiffs) const {
		for (int i = 0; i < activeColumnInputDiffs.size(); i++)
		{
			const int& activeColumnInputDiff = activeColumnInputDiffs[i];
			const vector<int> & possibleOutputDiffs = outputDifferencesPerInputDifferenceForChi[activeColumnInputDiff];
			int& activeColumnOutputDiffIndex = activeColumnOutputDiffIndexes[i];

			if (activeColumnOutputDiffIndex < possibleOutputDiffs.size() - 1)
			{
				activeColumnOutputDiffIndex++;
				return true;
			}
			else
			{
				activeColumnOutputDiffIndex = 0;
			}

		}

		return false;
	}

	// Generate the rotational difference for constant addition
	void GenerateRotationalDifferenceForConstantAddition() {
		for (int r = 0; r < 12; r++)
		{
			LaneType roundConstant = xoodooBase.GetRoundConstant(r);
			rotationalDifferencePerRoundConstants[r] = RotateLeft(roundConstant, rot) ^ roundConstant;
		}
	}

	// Generate DDT for the Chi operation in the round function of Xoodoo
	void GenerateDDTForChi() {
		for (int a = 0; a < 8; a++)
			for (int b = 0; b < 8; b++)
				ddt[a][b] = 0;

		for (int a = 0; a < 8; a++)
		{
			for (uint8_t x0 = 0; x0 < 8; x0++)
			{
				uint8_t x1 = x0 ^ a;
				uint8_t y0 = xoodooBase.ChiColumnByOperation(x0);
				uint8_t y1 = xoodooBase.ChiColumnByOperation(x1);
				uint8_t b = y0 ^ y1;

				ddt[a][b]++;
			}

		}

		for (int a = 0; a < 8; a++)
			for (int b = 0; b < 8; b++)
			{
				ddtProb[a][b] = ddt[a][b] / 8.0;
				ddtProbLog2[a][b] = log2(ddtProb[a][b]);
			}
	}

	// Generate possible output differences for each input difference of Chi operation
	void GenerateOutputDifferencesPerInputDifferenceForChi() {
		outputDifferencesPerInputDifferenceForChi.resize(8);
		for (int a = 0; a < 8; a++)
		{
			for (int b = 0; b < 8; b++)
				if (ddt[a][b] > 0)
					outputDifferencesPerInputDifferenceForChi[a].emplace_back(b);
		}
	}

	

public:
	RotationalDifferentialBase(const XoodooBase& xoodooBase, int rot, int milpThreads = 2) : xoodooBase(xoodooBase), rot(rot), milpThreads(milpThreads) {
		GenerateRotationalDifferenceForConstantAddition();
		GenerateDDTForChi();
		GenerateOutputDifferencesPerInputDifferenceForChi();
		milpChiInequalities = { {0, -1, 0, 0, 1, 1, 1},{0, 0, 0, 1, 1, 1, -1},{0, 0, 1, 0, 1, -1, 1},{0, 1, -1, 0, 0, 1, 1},{0, 1, 1, 0, -1, 0, 1},{0, 1, 1, 1, 0, 0, -1},{3, -1, -1, 0, -1, -1, 1},{3, -1, 0, -1, -1, 1, -1},{3, 0, -1, -1, 1, -1, -1},{0, -1, 1, 1, 1, 0, 0},{0, 0, 0, -1, 1, 1, 1},{0, 1, -1, 1, 0, 1, 0},{3, -1, -1, 1, -1, -1, 0},{3, -1, 1, -1, -1, 0, -1},{3, 1, -1, -1, 0, -1, -1}
		};
	}

	// Given the input state difference, compute the output state difference after propagating through the round linear layer of xoodoo
	void PropagateThroughLinearLayerWithProb(StateTypeBySheet& inputStateDiff, double& prob, const LaneType rcDiff) const {
		xoodooBase.LinearLayer(inputStateDiff, rcDiff);
	}

	void PropagateThroughLinearLayer(StateTypeBySheet& inputStateDiff, const LaneType rcDiff) const {
		xoodooBase.LinearLayer(inputStateDiff, rcDiff);
	}

	array<LaneType, 12> getRotationalDifferencePerRoundConstants() const {
		return rotationalDifferencePerRoundConstants;
	}

	vector<vector<int>> getOutputDiffsPerInputDiffForChi() const {
		return outputDifferencesPerInputDifferenceForChi;
	}


	array<array<double, 8>, 8> getDDTProb() const {
		return ddtProb;
	}

	int getRotationNumber() const {
		return rot;
	}

	int getActiveColumnsNum(const StateTypeBySheet& state) const
	{
		int activeColumnsNum = 0;
		for (int x = 0; x < 4; x++)
		{
			for (int z = 0; z < 32; z++)
			{
				if (xoodooBase.GetColumn(state, x, z) != 0)
					activeColumnsNum++;
			}
		}
		return activeColumnsNum;
	}

	// Given the input state difference of Chi at round startr and the output difference of Chi at round endr, compute the probability of the differential manually
    double ComputeDifferentialProbabilityManually(const StateTypeBySheet& inputStateDiff, const StateTypeBySheet& outputStateDiff, int startr, int endr, int totalRounds) const {
        int r = endr - startr;

       
		if (r < 0 || startr < 0 || endr < 0 || startr > totalRounds || endr > totalRounds)
		{
			throw MyException(__func__ + string(": Parameters error"));
		}
		else if (r == 0)
		{
			double difProb = 1.0;

			int activeSboxesNum = 0;
			for (int x = 0; x < 4; x++)
			{
				for (int z = 0; z < 32; z++)
				{
					uint8_t inputColumnDiff = xoodooBase.GetColumn(inputStateDiff, x, z);
					uint8_t outputColumnDiff = xoodooBase.GetColumn(outputStateDiff, x, z);
					difProb *= ddtProb[inputColumnDiff][outputColumnDiff];
					
					if (inputColumnDiff != 0)
						activeSboxesNum++;
					
				}
			}
			cout << activeSboxesNum << endl;
			return difProb;
		}
        else
        {
			double difProb = 0.0;

			StateTypeBySheet nextInputStateDiffBase = { 0 };
            vector<pair<int, int>> activeColumns;
            vector<int> activeColumnOutputDiffIndexes;
            vector<int> activeColumnInputDiffs;

			int rotationalDifferenceIndexNextRoundConstant = xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + 1);

            for(int x = 0; x < 4; x++)
            {
                for (int z = 0; z < 32; z++)
                {
                    uint8_t inputColumnDiff = xoodooBase.GetColumn(inputStateDiff, x, z);
                    if (inputColumnDiff == 0)
                        xoodooBase.SetColumn(nextInputStateDiffBase, x, z, 0);
                    else
                    {
                        activeColumns.emplace_back(x, z);
                        activeColumnOutputDiffIndexes.emplace_back(0);
                        activeColumnInputDiffs.emplace_back(inputColumnDiff);
                    }
                }
            }

            int activeColumnsNum = activeColumns.size();

			bool isAllActiveColumnsProcessed = false;
            while (!isAllActiveColumnsProcessed)
            {
                StateTypeBySheet nextInputStateDiff = nextInputStateDiffBase;

                double probCurNonLinearLayer = 1.0;

                for (int activeColumnIndex = 0; activeColumnIndex < activeColumnsNum; activeColumnIndex++)
                {
                    const pair<int,int> & activeColumnPos = activeColumns[activeColumnIndex];
                    const int & activeColumnOutputDiffIndex = activeColumnOutputDiffIndexes[activeColumnIndex];
                    const int & activeColumnInputDiff = activeColumnInputDiffs[activeColumnIndex];
                    const int & activeColumnOutputDiff = outputDifferencesPerInputDifferenceForChi[activeColumnInputDiff][activeColumnOutputDiffIndex];

                    xoodooBase.SetColumn(nextInputStateDiff, activeColumnPos.first, activeColumnPos.second, activeColumnOutputDiff);

                    probCurNonLinearLayer *= ddtProb[activeColumnInputDiff][activeColumnOutputDiff];
                }

                PropagateThroughLinearLayerWithProb(nextInputStateDiff, probCurNonLinearLayer, rotationalDifferencePerRoundConstants[rotationalDifferenceIndexNextRoundConstant]);
				double nextDifProb = ComputeDifferentialProbabilityManually(nextInputStateDiff, outputStateDiff, startr + 1, endr, totalRounds);

				
				if (nextDifProb > 0.0)
				{
					difProb += probCurNonLinearLayer * nextDifProb;
				}

                // switch to next possible output differences for the active columns
                if (SwitchToNextPossibleOutputDifferences(activeColumnOutputDiffIndexes, activeColumnInputDiffs))
                    isAllActiveColumnsProcessed = false;
                else
                    isAllActiveColumnsProcessed = true;
            }

            return difProb;
        }

		

    }

	double ExperimentalComputation(const StateTypeBySheet& stateInputDiff, const StateTypeBySheet& stateOutputDiff, int startr, int endr, int totalRounds) const
	{
		random_device rd;
		mt19937 gen(rd());
		uniform_int_distribution<uint32_t> distribution(0, numeric_limits<uint32_t>::max());

		int totalTestNum = 1 << 29;
		int zeroCount = 0;
		int rot = getRotationNumber();
		cout << "Rot is set to " << rot << endl;
		for (int testIndex = 0; testIndex < totalTestNum; testIndex++)
		{
			StateTypeBySheet state0;
			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					state0[x][y] = distribution(gen);

			StateTypeBySheet state1;
			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
				{
					state1[x][y] = state0[x][y] ^ stateInputDiff[x][y];
					state1[x][y] = RotateRight(state1[x][y], rot);
				}

			xoodooBase.ReducedRoundForRotDiffVerification(state0, startr, endr, totalRounds);
			xoodooBase.ReducedRoundForRotDiffVerification(state1, startr, endr, totalRounds);

			StateTypeBySheet outputStateDiff;
			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					outputStateDiff[x][y] = state0[x][y] ^ RotateLeft(state1[x][y], rot);

			bool isEqual = true;
			for (int x = 0; x < 4; x++)
			{
				for (int y = 0; y < 3; y++)
					if (outputStateDiff[x][y] != stateOutputDiff[x][y])
					{
						isEqual = false;
						break;
					}

				if (!isEqual)
					break;
			}

			if (isEqual)
			{
				// cout << "Find equal." << endl;
				zeroCount++;
			}
		}



		double zeroProb = zeroCount / (double)totalTestNum;
		// cout << zeroCount << endl;
		// cout << totalTestNum << endl;
		// cout << zeroProb << endl;

		return zeroProb;
	}

	// Given the input state differences and output state differences of Chi at each round, compute the probability of this trail
	double ComputeDifferentialTrailProbability(const vector<StateTypeBySheet>& inputStateDiffs, const vector<StateTypeBySheet>& outputStateDiffs, int startr, int endr, int totalRounds) const
	{
		int r = endr - startr;

		double trailProb = 1.0;
		if (r < 0 || startr < 0 || endr < 0 || startr > totalRounds || endr > totalRounds || inputStateDiffs.size() != endr - startr + 1 || outputStateDiffs.size() != endr - startr + 1)
		{
			throw MyException(__func__ + string(": Parameters error"));
		}
		else
		{
			for (int i = 0; i < r; i++)
			{
				StateTypeBySheet curInputStateDiff = inputStateDiffs[i];
				StateTypeBySheet curOutputStateDiff = outputStateDiffs[i];

				double probCurNonLinearLayer = ComputeDifferentialProbabilityManually(curInputStateDiff, curOutputStateDiff, startr + i, startr + i, totalRounds);
				trailProb *= probCurNonLinearLayer;
			}

			trailProb *= ComputeDifferentialProbabilityManually(inputStateDiffs[r], outputStateDiffs[r], startr + r, endr, totalRounds);

			return trailProb;
		}
	}


	// Given the input state differences and output state differences of Chi at each round, verify the validity and the probability of this trail
	bool VerifyDifferentialTrail(const vector<StateTypeBySheet>& inputStateDiffs, const vector<StateTypeBySheet> & outputStateDiffs, int startr, int endr, int totalRounds, double toVerifyTrailProb) const
	{
		int r = endr - startr;

		if (r < 0 || startr < 0 || endr < 0 || startr > totalRounds || endr > totalRounds || inputStateDiffs.size() != endr - startr + 1 || outputStateDiffs.size() != endr - startr + 1)
		{
			throw MyException(__func__ + string(": Parameters error"));
		}
		else
		{
			for (int i = 0; i < r; i++)
			{
				StateTypeBySheet curInputStateDiff = inputStateDiffs[i];
				StateTypeBySheet nextInputStateDiff = inputStateDiffs[i + 1]; 
				StateTypeBySheet curOutputStateDiff = outputStateDiffs[i];

				int curRotationalDifferenceIndexNextRoundConstant = xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + i + 1);

				PropagateThroughLinearLayer(curOutputStateDiff, rotationalDifferencePerRoundConstants[curRotationalDifferenceIndexNextRoundConstant]);
				if(curOutputStateDiff != nextInputStateDiff)
					return false;
			}

			double trailProb = ComputeDifferentialTrailProbability(inputStateDiffs, outputStateDiffs, startr, endr, totalRounds);

			if (abs(log2(trailProb) - log2(toVerifyTrailProb)) < 1e-2)
				return true;
			else
				return false;
		}
	}

// =====================================================================================
// ================= Functions below are for the milp models of RD propagation through xoodoo ==================
// =====================================================================================
	// Model the xor operation on three inputs
	void MilpXor3(GRBModel& model, const GRBVar& a0, const GRBVar& a1, const GRBVar& a2, const GRBVar& b) const
	{
		model.addConstr(b - a2 - a1 - a0 >= -2);
		model.addConstr(-b + a2 - a1 - a0 >= -2);
		model.addConstr(-b - a2 + a1 - a0 >= -2);
		model.addConstr(b + a2 + a1 - a0 >= 0);
		model.addConstr(-b - a2 - a1 + a0 >= -2);
		model.addConstr(b + a2 - a1 + a0 >= 0);
		model.addConstr(b - a2 + a1 + a0 >= 0);
		model.addConstr(-b + a2 + a1 + a0 >= 0);
	}


	RDStateMILPTypeBySheet MilpAddStateVars(GRBModel& model) const
	{
		RDStateMILPTypeBySheet stateDiffVars;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					stateDiffVars[x][y][z] = model.addVar(0, 1, 0, GRB_BINARY);
		return stateDiffVars;
	}

	StateTypeBySheet MilpReadStateVars(RDStateMILPTypeBySheet& stateDiffVars) const
	{
		StateTypeBySheet stateDiff;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					xoodooBase.SetBit(stateDiff, x, y, z, round(stateDiffVars[x][y][z].get(GRB_DoubleAttr_Xn)));
		return stateDiff;
	}

	RDPlaneMILPType MilpAddPlaneVars(GRBModel& model) const
	{
		RDPlaneMILPType planeDiffVars;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 32; z++)
				planeDiffVars[x][z] = model.addVar(0, 1, 0, GRB_BINARY);
		return planeDiffVars;
	}


	void MilpPropagateThroughRhoEast(RDStateMILPTypeBySheet& stateDiffVars) const
	{
		RDStateMILPTypeBySheet tmpStateDiffVars = stateDiffVars;
		for (int x = 0; x < 4; x++)
			stateDiffVars[x][1] = RDRotateLeft(tmpStateDiffVars[x][1], 1);
		
		for (int x = 0; x < 4; x++)
			stateDiffVars[x][2] = RDRotateLeft(tmpStateDiffVars[(x + 2) % 4][2], 8);
	}

	void MilpPropagateThroughTheta(GRBModel& model, RDStateMILPTypeBySheet& stateDiffVars) const
	{
		RDPlaneMILPType parityDiffVars = MilpAddPlaneVars(model);
		for(int x = 0; x < 4; x++)
			for(int z = 0; z < 32; z++)
				MilpXor3(model, stateDiffVars[x][0][z], stateDiffVars[x][1][z], stateDiffVars[x][2][z], parityDiffVars[x][z]);

		RDStateMILPTypeBySheet newStateDiffVars = MilpAddStateVars(model);

		for (int x = 0; x < 4; x++)
		{
			RDLaneMILPType one = RDRotateLeft(parityDiffVars[(x + 3) % 4], 5);
			RDLaneMILPType two = RDRotateLeft(parityDiffVars[(x + 3) % 4], 14);
			for(int y = 0; y < 3; y++)
				for(int z = 0; z < 32; z++)
					MilpXor3(model, stateDiffVars[x][y][z], one[z], two[z], newStateDiffVars[x][y][z]);
		}

		stateDiffVars = newStateDiffVars;
	}

	void MilpPropagateThroughRhoWest(RDStateMILPTypeBySheet& stateDiffVars) const
	{
		RDStateMILPTypeBySheet tmpStateDiffVars = stateDiffVars;
		for (int x = 0; x < 4; x++)
			stateDiffVars[x][1] = tmpStateDiffVars[(x + 3) % 4][1];

		for (int x = 0; x < 4; x++)
			stateDiffVars[x][2] = RDRotateLeft(tmpStateDiffVars[x][2], 11);
	}

	void MilpPropagateThroughLota(GRBModel& model, RDStateMILPTypeBySheet& stateDiffVars, const LaneType rcDiff) const
	{
		for(int z = 0; z < 32; z++)
			if ((rcDiff >> z) & 1)
			{
				GRBVar newBitVar = model.addVar(0, 1, 0, GRB_BINARY);
				model.addConstr(newBitVar == 1 - stateDiffVars[0][0][z]);
				stateDiffVars[0][0][z] = newBitVar;
			}
	}

	void MilpPropagateThroughChiColumn(GRBModel& model, GRBVar& v0, GRBVar& v1, GRBVar& v2) const
	{
		GRBVar u0 = model.addVar(0, 1, 0, GRB_BINARY);
		GRBVar u1 = model.addVar(0, 1, 0, GRB_BINARY);
		GRBVar u2 = model.addVar(0, 1, 0, GRB_BINARY);

		for (auto& ieq : milpChiInequalities)
		{
			GRBLinExpr ieqExp = ieq[0];
			ieqExp += v0 * ieq[1] + v1 * ieq[2] + v2 * ieq[3];
			ieqExp += u0 * ieq[4] + u1 * ieq[5] + u2 * ieq[6];
			model.addConstr(ieqExp >= 0);
		}

		v0 = u0;
		v1 = u1;
		v2 = u2;
	}

	void MilpPropagateThroughChi(GRBModel& model, RDStateMILPTypeBySheet& stateDiffVars) const
	{
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 32; z++)
				MilpPropagateThroughChiColumn(model, stateDiffVars[x][0][z], stateDiffVars[x][1][z], stateDiffVars[x][2][z]);
	}

	


	// Model the rotational difference propagation through the linear layer of Xoodoo
	void MilpPropagateThroughLinearLayer(GRBModel& model, RDStateMILPTypeBySheet& stateDiffVars, const LaneType rcDiff) const
	{
		MilpPropagateThroughRhoEast(stateDiffVars);
		MilpPropagateThroughTheta(model, stateDiffVars);
		MilpPropagateThroughRhoWest(stateDiffVars);
		MilpPropagateThroughLota(model, stateDiffVars, rcDiff);
	}

	void MilpPropagateThroughNonLinearLayer(GRBModel& model, RDStateMILPTypeBySheet& stateDiffVars) const
	{
		MilpPropagateThroughChi(model, stateDiffVars);
	}


	// Given the input state difference of Chi at round startr and the output difference of Chi at round endr, compute the probability of the differential automatically (using MILP model)
	double ComputeDifferentialProbabilityAutomatically(const StateTypeBySheet& inputStateDiff, const StateTypeBySheet& outputStateDiff, int startr, int endr, int totalRounds) const
	{
		int r = endr - startr;

		
		if (r < 0 || startr < 0 || endr < 0 || startr > totalRounds || endr > totalRounds)
		{
			throw MyException(__func__ + string(": Parameters error"));
		}
		else if (r == 0)
		{
			double difProb = 1.0;
			for (int x = 0; x < 4; x++)
			{
				for (int z = 0; z < 32; z++)
				{
					uint8_t inputColumnDiff = xoodooBase.GetColumn(inputStateDiff, x, z);
					uint8_t outputColumnDiff = xoodooBase.GetColumn(outputStateDiff, x, z);
					difProb *= ddtProb[inputColumnDiff][outputColumnDiff];
				}
			}
			return difProb;
		}
		else
		{

			GRBEnv env = GRBEnv();

			env.set(GRB_IntParam_LogToConsole, 0);
			env.set(GRB_IntParam_PoolSearchMode, 2);
			env.set(GRB_IntParam_PoolSolutions, 200000000);
			env.set(GRB_IntParam_Threads, milpThreads);

			GRBModel model = GRBModel(env);

			RDStateMILPTypeBySheet stateDiffVars = MilpAddStateVars(model);

			for(int x = 0; x < 4; x++)
				for(int y  = 0; y < 3; y++)
					for(int z = 0; z < 32; z++)
						model.addConstr(stateDiffVars[x][y][z] == xoodooBase.GetBit(inputStateDiff, x, y, z));


			vector<RDStateMILPTypeBySheet> inputStateDiffVarsPerRoundChi;
			vector<RDStateMILPTypeBySheet> outputStateDiffVarsPerRoundChi;

			inputStateDiffVarsPerRoundChi.emplace_back(stateDiffVars);
			MilpPropagateThroughNonLinearLayer(model, stateDiffVars);
			outputStateDiffVarsPerRoundChi.emplace_back(stateDiffVars);

			int midR = (startr + endr + 1) / 2;
			GRBLinExpr objExp = 0;

			for (int i = startr+1; i <= endr; i++)
			{
				int rotationalDifferenceIndexNextRoundConstant = xoodooBase.ComputeRoundConstantIndex(totalRounds, i);
				MilpPropagateThroughLinearLayer(model, stateDiffVars, rotationalDifferencePerRoundConstants[rotationalDifferenceIndexNextRoundConstant]);

				if(i == midR)
					for(int x =0; x < 4; x++)
						for(int y = 0; y < 3; y++)
							for(int z = 0; z < 32; z++)
								objExp += stateDiffVars[x][y][z];

				inputStateDiffVarsPerRoundChi.emplace_back(stateDiffVars);
				MilpPropagateThroughNonLinearLayer(model, stateDiffVars);
				outputStateDiffVarsPerRoundChi.emplace_back(stateDiffVars);
			}

			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					for (int z = 0; z < 32; z++)
						model.addConstr(stateDiffVars[x][y][z] == xoodooBase.GetBit(outputStateDiff, x, y, z));

			model.setObjective(objExp, GRB_MINIMIZE);

			model.optimize();

			int solCount = model.get(GRB_IntAttr_SolCount);
			cout << "Get " << solCount << " differential trails." << endl;

			if (solCount == 0)
				return 0.0;
			else
			{
				double difProb = 0.0;
				for (int solIndex = 0; solIndex < solCount; solIndex++)
				{
					// cout << "Trail " << solIndex << endl;

					model.set(GRB_IntParam_SolutionNumber, solIndex);
					vector<StateTypeBySheet> inputStateDiffs, outputStateDiffs;
					for (int i = 0; i <= r; i++)
					{
						StateTypeBySheet inputStateDiff = MilpReadStateVars(inputStateDiffVarsPerRoundChi[i]);
						StateTypeBySheet outputStateDiff = MilpReadStateVars(outputStateDiffVarsPerRoundChi[i]);

						inputStateDiffs.emplace_back(inputStateDiff);
						outputStateDiffs.emplace_back(outputStateDiff);
						// xoodooBase.OutputState(inputStateDiff);
						// xoodooBase.OutputState(outputStateDiff);
					}

					
					difProb += ComputeDifferentialTrailProbability(inputStateDiffs, outputStateDiffs, startr, endr, totalRounds);
				}

				return difProb;
			}
			
			
		}
	}
};

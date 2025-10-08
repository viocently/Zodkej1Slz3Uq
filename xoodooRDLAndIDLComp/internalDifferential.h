#pragma once
#include<iostream>
#include<array>
#include<vector>
#include"xoodoo.h"
#include"gurobi_c++.h"
#include"differential.h"

using namespace std;

using IDLaneMILPType = array<GRBVar, 16>;
using IDSheetMILPType = array<IDLaneMILPType, 3>;
using IDPlaneMILPType = array<IDLaneMILPType, 4>;
using IDStateMILPTypeBySheet = array<IDSheetMILPType, 4>;

IDLaneMILPType IDRotateLeft(IDLaneMILPType laneVars, int shift) {
	IDLaneMILPType newLaneVars;
	for (int i = 0; i < 16; i++)
		newLaneVars[i] = laneVars[(i - shift + 16) % 16];
	return newLaneVars;
}




class InternalDifferentialBase
{
private:
	const XoodooBase& xoodooBase;
	const RotationalDifferentialBase& rdBase;

	int offset;
	int milpThreads;
	array<HalfLaneType, 12> internalDifferencePerRoundConstants;
	array<array<double, 8>, 8> ddtProb;
	vector<vector<int>> outputDifferencesPerInputDifferenceForChi;
	vector<array<int, 7>> milpChiInequalities;

	bool SwitchToNextPossibleOutputDifferences(vector<int>& activeColumnOutputDiffIndexes, const vector<int>& activeColumnInputDiffs) const {
		for (int i = 0; i < activeColumnInputDiffs.size(); i++)
		{
			const int& activeColumnInputDiff = activeColumnInputDiffs[i];
			const vector<int>& possibleOutputDiffs = outputDifferencesPerInputDifferenceForChi[activeColumnInputDiff];
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

	void GenerateInternalDifferenceForConstantAddition() {
		for (int r = 0; r < 12; r++)
		{
			LaneType roundConstant = xoodooBase.GetRoundConstant(r);
			internalDifferencePerRoundConstants[r] = RotateLeft((roundConstant & 0x0000ffff) ^ ((roundConstant >> 16) & 0x0000ffff), offset);
		}
	}


public:
	InternalDifferentialBase(const XoodooBase& xoodooBase, const RotationalDifferentialBase & rdBase, int offset = 0, int milpThreads = 2) : xoodooBase(xoodooBase), rdBase(rdBase), offset(offset), milpThreads(milpThreads)
	{
		GenerateInternalDifferenceForConstantAddition();
		ddtProb = rdBase.getDDTProb();
		outputDifferencesPerInputDifferenceForChi = rdBase.getOutputDiffsPerInputDiffForChi();

		milpChiInequalities = { {0, -1, 0, 0, 1, 1, 1},{0, 0, 0, 1, 1, 1, -1},{0, 0, 1, 0, 1, -1, 1},{0, 1, -1, 0, 0, 1, 1},{0, 1, 1, 0, -1, 0, 1},{0, 1, 1, 1, 0, 0, -1},{3, -1, -1, 0, -1, -1, 1},{3, -1, 0, -1, -1, 1, -1},{3, 0, -1, -1, 1, -1, -1},{0, -1, 1, 1, 1, 0, 0},{0, 0, 0, -1, 1, 1, 1},{0, 1, -1, 1, 0, 1, 0},{3, -1, -1, 1, -1, -1, 0},{3, -1, 1, -1, -1, 0, -1},{3, 1, -1, -1, 0, -1, -1}
		};
	}

	array<array<double, 8>, 8> getDDTProb() const
	{
		return ddtProb;
	}

	int getActiveColumnsNum(const HalfStateTypeBySheet& state) const
	{
		int activeColumnsNum = 0;
		for (int x = 0; x < 4; x++)
		{
			for (int z = 0; z < 16; z++)
			{
				if (xoodooBase.GetColumn(state, x, z) != 0)
					activeColumnsNum++;
			}
		}
		return activeColumnsNum;
	}


	double ExperimentalComputation(const HalfStateTypeBySheet& stateInputDiff, const HalfStateTypeBySheet& stateOutputDiff, int startr, int endr, int totalRounds) const
	{
		random_device rd;
		mt19937 gen(rd());
		uniform_int_distribution<uint16_t> distribution(0, numeric_limits<uint16_t>::max());

		int totalTestNum = 1 << 25;
		int zeroCount = 0;
		for (int testIndex = 0; testIndex < totalTestNum; testIndex++)
		{
			StateTypeBySheet state0;
			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
				{
					state0[x][y] = distribution(gen);
					state0[x][y] ^= ((state0[x][y] ^ stateInputDiff[x][y]) << 16);
				}

			

			xoodooBase.ReducedRoundForRotDiffVerification(state0, startr, endr, totalRounds);

			HalfStateTypeBySheet outputStateDiff;
			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					outputStateDiff[x][y] = state0[x][y] ^ (state0[x][y] >> 16);

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

	void PropagateThroughRhoEast(HalfStateTypeBySheet& stateDiff) const
	{
		HalfStateTypeBySheet tmp(stateDiff);
		for (int x = 0; x < 4; x++)
			stateDiff[x][1] = RotateLeft(tmp[x][1], 1);

		for (int x = 0; x < 4; x++)
			stateDiff[x][2] = RotateLeft(tmp[(x + 2) % 4][2], 8);
	}

	void ProagateThroughTheta(HalfStateTypeBySheet& stateDiff) const
	{
		HalfPlaneType parity = { 0 };
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				parity[x] ^= stateDiff[x][y];

		for (int x = 0; x < 4; x++) {
			HalfLaneType one = RotateLeft(parity[(x + 3) % 4], 5);
			HalfLaneType two = RotateLeft(parity[(x + 3) % 4], 14);
			for (int y = 0; y < 3; y++)
				stateDiff[x][y] ^= one ^ two;
		}
	}

	void PropagateThroughRhoWest(HalfStateTypeBySheet& stateDiff) const
	{
		HalfStateTypeBySheet tmp(stateDiff);
		for (int x = 0; x < 4; x++)
			stateDiff[x][1] = tmp[(x + 3) % 4][1];

		for (int x = 0; x < 4; x++)
			stateDiff[x][2] = RotateLeft(tmp[x][2], 11);
	}

	void PropagateThroughLota(HalfStateTypeBySheet& stateDiff, const HalfLaneType rcDiff) const
	{
		stateDiff[0][0] ^= rcDiff;
	}

	void PropagateThroughLinearLayer(HalfStateTypeBySheet& stateDiff, const HalfLaneType rcDiff) const
	{
		PropagateThroughRhoEast(stateDiff);
		ProagateThroughTheta(stateDiff);
		PropagateThroughRhoWest(stateDiff);
		PropagateThroughLota(stateDiff, rcDiff);
	}

	array<HalfLaneType, 12> getInternalDifferencePerRoundConstants() const {
		return internalDifferencePerRoundConstants;
	}

	int getOffsetNumber() const {
		return offset;
	}

	int getActiveColumnNum(const HalfStateTypeBySheet & stateDiff) const {
		int activeColumnsNum = 0;
		for (int x = 0; x < 4; x++)
		{
			for (int z = 0; z < 16; z++)
			{
				if (xoodooBase.GetColumn(stateDiff, x, z) != 0)
					activeColumnsNum++;
			}
		}
		return activeColumnsNum;
	}

	double ComputeDifferentialProbabilityManually(const HalfStateTypeBySheet& inputStateDiff, const HalfStateTypeBySheet& outputStateDiff, int startr, int endr, int totalRounds) const {
		int r = endr - startr;


		if (r < 0 || startr < 0 || endr < 0 || startr > totalRounds || endr > totalRounds)
		{
			throw MyException(__func__ + string(": Parameters error"));
		}
		else if (r == 0)
		{
			double difProb = 1.0;

			// int activeSboxesNum = 0;
			for (int x = 0; x < 4; x++)
			{
				for (int z = 0; z < 16; z++)
				{
					uint8_t inputColumnDiff = xoodooBase.GetColumn(inputStateDiff, x, z);
					uint8_t outputColumnDiff = xoodooBase.GetColumn(outputStateDiff, x, z);
					difProb *= ddtProb[inputColumnDiff][outputColumnDiff];
					/*
					if (inputColumnDiff != 0)
						activeSboxesNum++;
					*/
				}
			}
			// cout << activeSboxesNum << endl;
			return difProb;
		}
		else
		{
			double difProb = 0.0;

			HalfStateTypeBySheet nextInputStateDiffBase = { 0 };
			vector<pair<int, int>> activeColumns;
			vector<int> activeColumnOutputDiffIndexes;
			vector<int> activeColumnInputDiffs;

			int rotationalDifferenceIndexNextRoundConstant = xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + 1);

			for (int x = 0; x < 4; x++)
			{
				for (int z = 0; z < 16; z++)
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
				HalfStateTypeBySheet nextInputStateDiff = nextInputStateDiffBase;

				double probCurNonLinearLayer = 1.0;

				for (int activeColumnIndex = 0; activeColumnIndex < activeColumnsNum; activeColumnIndex++)
				{
					const pair<int, int>& activeColumnPos = activeColumns[activeColumnIndex];
					const int& activeColumnOutputDiffIndex = activeColumnOutputDiffIndexes[activeColumnIndex];
					const int& activeColumnInputDiff = activeColumnInputDiffs[activeColumnIndex];
					const int& activeColumnOutputDiff = outputDifferencesPerInputDifferenceForChi[activeColumnInputDiff][activeColumnOutputDiffIndex];

					xoodooBase.SetColumn(nextInputStateDiff, activeColumnPos.first, activeColumnPos.second, activeColumnOutputDiff);

					probCurNonLinearLayer *= ddtProb[activeColumnInputDiff][activeColumnOutputDiff];
				}

				PropagateThroughLinearLayer(nextInputStateDiff, internalDifferencePerRoundConstants[rotationalDifferenceIndexNextRoundConstant]);
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

	double ComputeDifferentialTrailProbability(const vector<HalfStateTypeBySheet>& inputStateDiffs, const vector<HalfStateTypeBySheet>& outputStateDiffs, int startr, int endr, int totalRounds) const
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
				HalfStateTypeBySheet curInputStateDiff = inputStateDiffs[i];
				HalfStateTypeBySheet curOutputStateDiff = outputStateDiffs[i];

				double probCurNonLinearLayer = ComputeDifferentialProbabilityManually(curInputStateDiff, curOutputStateDiff, startr + i, startr + i, totalRounds);
				trailProb *= probCurNonLinearLayer;
			}

			trailProb *= ComputeDifferentialProbabilityManually(inputStateDiffs[r], outputStateDiffs[r], startr + r, endr, totalRounds);

			return trailProb;
		}
	}

	bool VerifyDifferentialTrail(const vector<HalfStateTypeBySheet>& inputStateDiffs, const vector<HalfStateTypeBySheet>& outputStateDiffs, int startr, int endr, int totalRounds, double toVerifyTrailProb) const
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
				HalfStateTypeBySheet curInputStateDiff = inputStateDiffs[i];
				HalfStateTypeBySheet nextInputStateDiff = inputStateDiffs[i + 1];
				HalfStateTypeBySheet curOutputStateDiff = outputStateDiffs[i];

				int curRotationalDifferenceIndexNextRoundConstant = xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + i + 1);

				PropagateThroughLinearLayer(curOutputStateDiff, internalDifferencePerRoundConstants[curRotationalDifferenceIndexNextRoundConstant]);
				if (curOutputStateDiff != nextInputStateDiff)
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
// ================= Functions below are for the milp models of ID propagation through xoodoo ==================
// =====================================================================================
	void MilpXor3(GRBModel& model, const GRBVar& a0, const GRBVar& a1, const GRBVar& a2, const GRBVar& b) const
	{
		rdBase.MilpXor3(model, a0, a1, a2, b);
	}

	IDStateMILPTypeBySheet MilpAddStateVars(GRBModel& model) const
	{
		IDStateMILPTypeBySheet stateDiffVars;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					stateDiffVars[x][y][z] = model.addVar(0, 1, 0, GRB_BINARY);
		return stateDiffVars;
	}

	HalfStateTypeBySheet MilpReadStateVars(IDStateMILPTypeBySheet& stateDiffVars) const
	{
		HalfStateTypeBySheet stateDiff;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					xoodooBase.SetBit(stateDiff, x, y, z, round(stateDiffVars[x][y][z].get(GRB_DoubleAttr_Xn)));
		return stateDiff;
	}

	IDPlaneMILPType MilpAddPlaneVars(GRBModel& model) const
	{
		IDPlaneMILPType planeDiffVars;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
				planeDiffVars[x][z] = model.addVar(0, 1, 0, GRB_BINARY);
		return planeDiffVars;
	}

	void MilpPropagateThroughRhoEast(IDStateMILPTypeBySheet& stateDiffVars) const
	{
		IDStateMILPTypeBySheet tmpStateDiffVars = stateDiffVars;
		for (int x = 0; x < 4; x++)
			stateDiffVars[x][1] = IDRotateLeft(tmpStateDiffVars[x][1], 1);

		for (int x = 0; x < 4; x++)
			stateDiffVars[x][2] = IDRotateLeft(tmpStateDiffVars[(x + 2) % 4][2], 8);
	}

	void MilpPropagateThroughTheta(GRBModel& model, IDStateMILPTypeBySheet& stateDiffVars) const
	{

		
		IDPlaneMILPType parityDiffVars = MilpAddPlaneVars(model);
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
				MilpXor3(model, stateDiffVars[x][0][z], stateDiffVars[x][1][z], stateDiffVars[x][2][z], parityDiffVars[x][z]);

		IDStateMILPTypeBySheet newStateDiffVars = MilpAddStateVars(model);

		for (int x = 0; x < 4; x++)
		{
			IDLaneMILPType one = IDRotateLeft(parityDiffVars[(x + 3) % 4], 5);
			IDLaneMILPType two = IDRotateLeft(parityDiffVars[(x + 3) % 4], 14);
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					MilpXor3(model, stateDiffVars[x][y][z], one[z], two[z], newStateDiffVars[x][y][z]);
		}

		stateDiffVars = newStateDiffVars;
	}

	void MilpPropagateThroughRhoWest(IDStateMILPTypeBySheet& stateDiffVars) const
	{
		IDStateMILPTypeBySheet tmpStateDiffVars = stateDiffVars;
		for (int x = 0; x < 4; x++)
			stateDiffVars[x][1] = tmpStateDiffVars[(x + 3) % 4][1];

		for (int x = 0; x < 4; x++)
			stateDiffVars[x][2] = IDRotateLeft(tmpStateDiffVars[x][2], 11);
	}

	void MilpPropagateThroughLota(GRBModel& model, IDStateMILPTypeBySheet& stateDiffVars, const HalfLaneType rcDiff) const
	{
		for (int z = 0; z < 16; z++)
			if ((rcDiff >> z) & 1)
			{
				GRBVar newBitVar = model.addVar(0, 1, 0, GRB_BINARY);
				model.addConstr(newBitVar == 1 - stateDiffVars[0][0][z]);
				stateDiffVars[0][0][z] = newBitVar;
			}
	}

	void MilpPropagateThroughChiColumn(GRBModel& model, GRBVar& v0, GRBVar& v1, GRBVar& v2) const
	{
		rdBase.MilpPropagateThroughChiColumn(model, v0, v1, v2);
	}

	void MilpPropagateThroughChi(GRBModel& model, IDStateMILPTypeBySheet& stateDiffVars) const
	{
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
				MilpPropagateThroughChiColumn(model, stateDiffVars[x][0][z], stateDiffVars[x][1][z], stateDiffVars[x][2][z]);
	}

	void MilpPropagateThroughLinearLayer(GRBModel& model, IDStateMILPTypeBySheet& stateDiffVars, const HalfLaneType rcDiff) const
	{
		MilpPropagateThroughRhoEast(stateDiffVars);
		MilpPropagateThroughTheta(model, stateDiffVars);
		MilpPropagateThroughRhoWest(stateDiffVars);
		MilpPropagateThroughLota(model, stateDiffVars, rcDiff);
	}

	void MilpPropagateThroughNonLinearLayer(GRBModel& model, IDStateMILPTypeBySheet& stateDiffVars) const
	{
		MilpPropagateThroughChi(model, stateDiffVars);
	}

	double ComputeDifferentialProbabilityAutomatically(const HalfStateTypeBySheet& inputStateDiff, const HalfStateTypeBySheet& outputStateDiff, int startr, int endr, int totalRounds) const
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
				for (int z = 0; z < 16; z++)
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

			IDStateMILPTypeBySheet stateDiffVars = MilpAddStateVars(model);

			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					for (int z = 0; z < 16; z++)
						model.addConstr(stateDiffVars[x][y][z] == xoodooBase.GetBit(inputStateDiff, x, y, z));


			vector<IDStateMILPTypeBySheet> inputStateDiffVarsPerRoundChi;
			vector<IDStateMILPTypeBySheet> outputStateDiffVarsPerRoundChi;

			inputStateDiffVarsPerRoundChi.emplace_back(stateDiffVars);
			MilpPropagateThroughNonLinearLayer(model, stateDiffVars);
			outputStateDiffVarsPerRoundChi.emplace_back(stateDiffVars);

			int midR = (startr + endr + 1) / 2;
			GRBLinExpr objExp = 0;

			for (int i = startr + 1; i <= endr; i++)
			{
				int rotationalDifferenceIndexNextRoundConstant = xoodooBase.ComputeRoundConstantIndex(totalRounds, i);
				MilpPropagateThroughLinearLayer(model, stateDiffVars, internalDifferencePerRoundConstants[rotationalDifferenceIndexNextRoundConstant]);

				if (i == midR)
					for (int x = 0; x < 4; x++)
						for (int y = 0; y < 3; y++)
							for (int z = 0; z < 16; z++)
								objExp += stateDiffVars[x][y][z];

				inputStateDiffVarsPerRoundChi.emplace_back(stateDiffVars);
				MilpPropagateThroughNonLinearLayer(model, stateDiffVars);
				outputStateDiffVarsPerRoundChi.emplace_back(stateDiffVars);
			}

			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					for (int z = 0; z < 16; z++)
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
					vector<HalfStateTypeBySheet> inputStateDiffs, outputStateDiffs;
					for (int i = 0; i <= r; i++)
					{
						HalfStateTypeBySheet inputStateDiff = MilpReadStateVars(inputStateDiffVarsPerRoundChi[i]);
						HalfStateTypeBySheet outputStateDiff = MilpReadStateVars(outputStateDiffVarsPerRoundChi[i]);
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

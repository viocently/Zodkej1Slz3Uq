#pragma once
#include<iostream>
#include<array>
#include<vector>
#include<random>
#include"ascon.h"
#include"gurobi_c++.h"


using namespace std;



using RDRowMILPType = array<GRBVar, 64>;
using RDStateMILPType = array<RDRowMILPType, 5>;

RDRowMILPType RDRotateLeft(const RDRowMILPType& value, int shift) {
	RDRowMILPType result;
	for (int i = 0; i < 64; i++) {
		result[i] = value[(i - shift + 64) % 64];
	}
	return result;
}

RDRowMILPType RDRotateRight(const RDRowMILPType& value, int shift) {
	return RDRotateLeft(value, 64 - shift);
}

class RotationalDifferentialBase
{
private:
	const AsconBase& asconBase;

	int rot;
	int milpThreads;
	array<RowType, 12> rotationalDifferencePerRoundConstant;

	array<array<int, 32>, 32> ddt;
	array<array<double, 32>, 32> ddtProb;
	array<array<double, 32>, 32> ddtProbLog2;
	vector<vector<int>> outputDifferencesPerInputDifferenceForSbox;

	vector<array<int, 11>> milpSboxInequalities;

	bool SwitchToNextPossibleOutputDifferences(vector<int>& activeColumnOutputDiffIndexes, const vector<int>& activeColumnInputDiffs) const {
		for (int i = 0; i < activeColumnInputDiffs.size(); i++)
		{
			const int& activeColumnInputDiff = activeColumnInputDiffs[i];
			const vector<int>& possibleOutputDiffs = outputDifferencesPerInputDifferenceForSbox[activeColumnInputDiff];
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

	void GenerateRotationalDifferenceForConstantAddition()
	{
		for (int i = 0; i < 12; i++)
		{
			RowType rc = asconBase.GetRoundConstant(i);
			RowType rcRot = RotateLeft(rc, rot);
			rotationalDifferencePerRoundConstant[i] = rc ^ rcRot;
		}
	}

	void GenerateDDTForSbox()
	{
		for (int a = 0; a < 32; a++)
			for (int b = 0; b < 32; b++)
			{
				ddt[a][b] = 0;
			}

		for (int a = 0; a < 32; a++)
			for (int x = 0; x < 32; x++)
			{
				int y = asconBase.Sbox(x);
				int xPrime = x ^ a;
				int yPrime = asconBase.Sbox(xPrime);
				int b = y ^ yPrime;
				ddt[a][b]++;
			}

		for (int a = 0; a < 32; a++)
			for (int b = 0; b < 32; b++)
			{
				ddtProb[a][b] = (double)ddt[a][b] / 32.0;
				ddtProbLog2[a][b] = log2(ddtProb[a][b]);
			}
	}

	void GenerateOutputDifferencesPerInputDifferenceForSbox()
	{
		outputDifferencesPerInputDifferenceForSbox.resize(32);
		for (int a = 0; a < 32; a++)
		{
			for (int b = 0; b < 32; b++)
			{
				if (ddt[a][b] > 0)
					outputDifferencesPerInputDifferenceForSbox[a].push_back(b);
			}
		}
	}

public:
	RotationalDifferentialBase(const AsconBase& asconBase, int rot, int milpThreads = 2) : asconBase(asconBase), rot(rot), milpThreads(milpThreads)
	{
		GenerateRotationalDifferenceForConstantAddition();
		GenerateDDTForSbox();
		GenerateOutputDifferencesPerInputDifferenceForSbox();

		milpSboxInequalities = { {0, 5, 9, 3, 9, 4, -6, 0, -1, 2, -2}, {0, -2, 4, 0, 3, -2, 3, 0, 2, 1, 1}, {0, 8, -5, -5, 7, 13, 13, 4, 4, -2, -1}, {7, -2, -3, 1, 3, -1, -2, -1, 0, 0, -1}, {0, 1, -2, 5, 5, 3, 3, -2, -2, 0, 3}, {0, 4, 5, 3, -1, 1, -2, -1, 6, 2, -2}, {2, 1, 1, -1, -1, 0, 0, 0, -1, 0, 0}, {0, -2, -3, -1, -1, -1, 5, 8, 8, 7, 5}, {8, -1, -2, -2, -2, 1, 0, -2, 0, -1, 0}, {9, -3, -2, 2, -1, -3, -2, 5, 5, -2, -1}, {26, -2, -3, -2, -4, -5, -3, -7, -7, 5, 1}, {4, -2, 4, 2, 3, -3, 1, 0, -1, -2, 1}, {1, -1, 1, -2, 2, 2, 0, 1, 0, 0, 2}, {14, -1, 5, -5, -1, -2, 6, 2, -3, -4, -4}, {6, 1, -1, 1, -2, 3, -3, 3, -3, 1, 1}, {13, -2, 2, -2, -4, -1, -3, -1, 1, -4, 4}, {2, 1, -1, -1, 1, -1, 0, 0, 0, 0, 1}, {5, -3, 4, 2, 5, -5, 7, -1, -1, 3, -2}, {2, 1, 2, -1, 2, 0, -1, 0, 1, -1, -1}, {6, -3, -2, -1, -1, 3, -1, 2, 0, 3, -1}, {0, 2, 1, 2, -2, 1, 1, 1, 1, 0, 0}, {23, -6, -4, -3, -2, -3, -1, -3, -3, -3, -1}, {3, 1, -1, -1, -1, -1, 0, 0, 0, 1, 0}, {0, -1, 2, -1, 0, 0, 2, 0, 1, 2, 2}, {6, -2, -2, -1, 1, 2, 0, -1, 0, 0, -2}, {6, 1, 3, -3, -1, 0, -2, 0, -1, 2, -2}, {2, -1, -1, 1, -1, 1, 0, 1, 0, 0, 0}, {2, 0, 1, 0, -1, 0, 1, 0, 1, -1, -1}, {4, 4, 2, 1, -2, 3, -3, -3, 5, -1, 1}, {0, -2, 3, 2, 1, -4, 1, 4, 3, 5, 5}, {8, 1, -2, 2, -2, -1, 2, -2, -1, -2, 0}, {12, -2, 2, 4, 1, 1, -4, 1, -3, -3, -4}, {0, 1, -5, 1, 3, 2, 4, 3, 2, 0, 4}, {4, -1, -2, 1, 2, -1, -1, 0, 0, 0, -1}, {2, 1, -1, 0, -1, -1, 1, 1, 0, 1, 0}, {2, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0}, {3, 1, -1, -1, -1, 1, 0, 0, 0, -1, 0}, {3, 1, -1, 0, -1, 0, -1, -1, 1, 1, 0}, {5, -2, -2, -1, 1, -1, -1, 1, 1, 0, 0}, {0, -2, 6, 5, 5, -1, 6, -1, -1, -2, 6}, {0, -1, 3, 3, 2, 0, 3, 0, -1, 3, -1}, {0, 1, 1, 1, -1, 0, 0, 0, 1, 0, 0}, {4, 1, -1, 1, -1, -1, -1, 1, 0, -1, 0}, {0, -4, 7, -1, 3, -1, 4, -1, 3, 3, 3}, {0, -2, 4, 1, 6, -1, -2, -1, 5, 5, 6}, {0, -2, 4, 6, 5, 0, -1, -2, -1, 6, 6}, {0, -1, 0, 1, 2, 3, 3, 2, -1, -1, 3}, {0, 2, 3, 2, -1, 2, -1, 0, 3, -1, 1}, {3, -1, 0, 0, 1, -1, 1, -1, -1, 0, 0} };

	}

	void PropagateThroughLinearLayer(StateType& inputStateDiff, const RowType rcDiff) const
	{
		asconBase.LinearLayer(inputStateDiff, rcDiff);
	}

	array<RowType, 12> GetRotationalDifferencePerRoundConstant() const
	{
		return rotationalDifferencePerRoundConstant;
	}


	vector<vector<int>> GetOutputDiffsPerInputDiffForSbox() const
	{
		return outputDifferencesPerInputDifferenceForSbox;
	}

	array<array<double, 32>, 32> GetDDTProb() const
	{
		return ddtProb;
	}

	int GetRotationNumber() const
	{
		return rot;
	}

	int GetActiveColumnsNum(const StateType& stateDiff) const
	{
		int activeCols = 0;
		for (int col = 0; col < 64; col++)
		{
			if (asconBase.GetColumn(stateDiff, col) != 0)
				activeCols++;
		}
		return activeCols;
	}

	double ComputeDifferentialProbabilityManually(const StateType& inputStateDiff, const StateType& outputStateDiff, int startr, int endr) const
	{
		int r = endr - startr;

		if (r == 0)
		{
			double difProb = 1.0;
			RowType constantDiff = rotationalDifferencePerRoundConstant[startr + 1];

			StateType sboxesInputDiff = inputStateDiff;

			StateType sboxesOutputDiff = outputStateDiff;

			asconBase.LinearLayerInverse(sboxesOutputDiff, constantDiff);



			for (int col = 0; col < 64; col++)
			{
				uint8_t inputColumnDiff = asconBase.GetColumn(sboxesInputDiff, col);
				uint8_t outputColumnDiff = asconBase.GetColumn(sboxesOutputDiff, col);

				difProb *= ddtProb[inputColumnDiff][outputColumnDiff];
			}

			return difProb;
		}
		else
		{
			double difProb = 0.0;

			StateType nextInputStateDiffBase = { 0 };
			vector<int> activeColumns;
			vector<int> activeColumnOutputDiffIndexes;
			vector<int> activeColumnInputDiffs;

			for (int col = 0; col < 64; col++)
			{
				uint8_t inputColumnDiff = asconBase.GetColumn(inputStateDiff, col);
				if (inputColumnDiff == 0)
					asconBase.SetColumn(nextInputStateDiffBase, col, 0);
				else
				{
					activeColumns.emplace_back(col);
					activeColumnOutputDiffIndexes.emplace_back(0);
					activeColumnInputDiffs.emplace_back(inputColumnDiff);
				}
			}

			int activeColumnsNum = activeColumns.size();

			bool isAllActiveColumnsProcessed = false;
			while (!isAllActiveColumnsProcessed)
			{
				StateType nextInputStateDiff = nextInputStateDiffBase;

				double probCurNonLinearLayer = 1.0;
				for (int activeColumnIndex = 0; activeColumnIndex < activeColumnsNum; activeColumnIndex++)
				{
					const int& activeColumnPos = activeColumns[activeColumnIndex];
					const int& activeColumnOutputDiffIndex = activeColumnOutputDiffIndexes[activeColumnIndex];
					const int& activeColumnInputDiff = activeColumnInputDiffs[activeColumnIndex];
					const int& activeColumnOutputDiff = outputDifferencesPerInputDifferenceForSbox[activeColumnInputDiff][activeColumnOutputDiffIndex];

					asconBase.SetColumn(nextInputStateDiff, activeColumnPos, activeColumnOutputDiff);

					probCurNonLinearLayer *= ddtProb[activeColumnInputDiff][activeColumnOutputDiff];
				}

				PropagateThroughLinearLayer(nextInputStateDiff, rotationalDifferencePerRoundConstant[startr + 1]);

				double nextDifProb = ComputeDifferentialProbabilityManually(nextInputStateDiff, outputStateDiff, startr + 1, endr);


				if (nextDifProb > 0)
					difProb += probCurNonLinearLayer * nextDifProb;

				// switch to next possible output differences for the active columns
				if (SwitchToNextPossibleOutputDifferences(activeColumnOutputDiffIndexes, activeColumnInputDiffs))
					isAllActiveColumnsProcessed = false;
				else
					isAllActiveColumnsProcessed = true;

			}


			return difProb;
		}

	}

	double ExperimentalComputation(const StateType& stateInputDiff, const StateType& stateOutputDiff, int startr, int endr) const
	{
		random_device rd;
		mt19937 gen(rd());
		uniform_int_distribution<uint64_t> distribution(0, numeric_limits<uint64_t>::max());


		int totalTestNum = 1 << 20;
		int zeroCount = 0;
		int rot = GetRotationNumber();
		cout << "Rot is set to " << rot << endl;
		for (int testIndex = 0; testIndex < totalTestNum; testIndex++)
		{
			StateType state0;
			for (int row = 0; row < 5; row++)
				state0[row] = distribution(gen);

			StateType state1;
			for (int row = 0; row < 5; row++)
			{
				state1[row] = state0[row] ^ stateInputDiff[row];
				state1[row] = RotateRight(state1[row], rot);
			}

			asconBase.ReducedRoundsForVerification(state0, startr, endr);
			asconBase.ReducedRoundsForVerification(state1, startr, endr);

			StateType outputStateDiff;
			for (int row = 0; row < 5; row++)
				outputStateDiff[row] = state0[row] ^ RotateLeft(state1[row], rot);

			bool isEqual = true;
			for (int row = 0; row < 5; row++)
			{
				if (outputStateDiff[row] != stateOutputDiff[row])
				{
					isEqual = false;
					break;
				}
			}

			if (isEqual)
				zeroCount++;
		}

		double zeroProb = (double)zeroCount / (double)totalTestNum;

		return zeroProb;
	}

	double ComputeDifferentialTrailProbability(const vector<StateType>& stateDiffs, int startr, int endr) const
	{
		int r = endr - startr;


		double trailProb = 1.0;
		for (int i = 0; i < r + 1; i++)
		{
			StateType curInputStateDiff = stateDiffs[i];
			StateType curOutputStateDiff = stateDiffs[i + 1];

			double curDifProb = ComputeDifferentialProbabilityManually(curInputStateDiff, curOutputStateDiff, startr + i, startr + i);

			trailProb *= curDifProb;
		}

		return trailProb;
	}

	//=====================================================================================
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


	RDStateMILPType MilpAddStateVars(GRBModel& model) const
	{
		RDStateMILPType stateVars;
		for (int row = 0; row < 5; row++)
		{
			for (int col = 0; col < 64; col++)
			{
				stateVars[row][col] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}
		return stateVars;
	}

	StateType MilpReadStateVars(const RDStateMILPType& stateVars) const
	{
		StateType state;
		for (int row = 0; row < 5; row++)
		{
			for (int col = 0; col < 64; col++)
			{
				asconBase.SetBit(state, row, col, round(stateVars[row][col].get(GRB_DoubleAttr_Xn)));
			}
		}
		return state;
	}

	void MilpPropagateThroughRotateRows(GRBModel& model, RDStateMILPType& stateVars) const
	{
		for (int row = 0; row < 5; row++)
		{
			int rot0 = asconBase.GetRotationNum(row, 0);
			int rot1 = asconBase.GetRotationNum(row, 1);
			auto tmp0 = RDRotateRight(stateVars[row], rot0);
			auto tmp1 = RDRotateRight(stateVars[row], rot1);
			for (int col = 0; col < 64; col++)
			{
				GRBVar newVar = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				MilpXor3(model, stateVars[row][col], tmp0[col], tmp1[col], newVar);
				stateVars[row][col] = newVar;
			}
		}
	}

	void MilpPropagateThroughAddRoundConstant(GRBModel& model, RDStateMILPType& stateVars, const RowType rcDiff) const
	{
		for (int col = 0; col < 64; col++)
		{
			if ((rcDiff >> col) & 1ull)
			{
				GRBVar newBitVar = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				model.addConstr(newBitVar == 1 - stateVars[2][col]);
				stateVars[2][col] = newBitVar;
			}
		}
	}

	void MilpPropagateThroughSbox(GRBModel& model, GRBVar& v0, GRBVar& v1, GRBVar& v2, GRBVar& v3, GRBVar& v4) const
	{
		GRBVar u0 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		GRBVar u1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		GRBVar u2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		GRBVar u3 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		GRBVar u4 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

		for (auto& ieq : milpSboxInequalities)
		{
			GRBLinExpr ieqExp = ieq[0];
			ieqExp += v0 * ieq[1] + v1 * ieq[2] + v2 * ieq[3] + v3 * ieq[4] + v4 * ieq[5];
			ieqExp += u0 * ieq[6] + u1 * ieq[7] + u2 * ieq[8] + u3 * ieq[9] + u4 * ieq[10];
			model.addConstr(ieqExp >= 0);
		}

		v0 = u0;
		v1 = u1;
		v2 = u2;
		v3 = u3;
		v4 = u4;
	}

	void MilpPropagateThroughSboxes(GRBModel& model, RDStateMILPType& stateVars) const
	{
		for (int col = 0; col < 64; col++)
		{
			MilpPropagateThroughSbox(model, stateVars[4][col], stateVars[3][col], stateVars[2][col], stateVars[1][col], stateVars[0][col]);
		}
	}

	void MilpPropagateThroughLinearLayer(GRBModel& model, RDStateMILPType& stateVars, const RowType rcDiff) const
	{
		MilpPropagateThroughRotateRows(model, stateVars);
		MilpPropagateThroughAddRoundConstant(model, stateVars, rcDiff);
	}

	void MilpPropagateThroughNonLinearLayer(GRBModel& model, RDStateMILPType& stateVars) const
	{
		MilpPropagateThroughSboxes(model, stateVars);
	}

	double ComputeDifferentialProbabilityAutomatically(const StateType& inputStateDiff, const StateType& outputStateDiff, int startr, int endr) const
	{
		int r = endr - startr;

		if (r == 0)
			return ComputeDifferentialProbabilityManually(inputStateDiff, outputStateDiff, startr, endr);
		else
		{
			GRBEnv env = GRBEnv();

			env.set(GRB_IntParam_LogToConsole, 0);
			env.set(GRB_IntParam_PoolSearchMode, 2);
			env.set(GRB_IntParam_PoolSolutions, 200000000);
			env.set(GRB_IntParam_Threads, milpThreads);

			GRBModel model = GRBModel(env);

			RDStateMILPType stateDiffVars = MilpAddStateVars(model);

			for (int row = 0; row < 5; row++)
				for (int col = 0; col < 64; col++)
					model.addConstr(stateDiffVars[row][col] == asconBase.GetBit(inputStateDiff, row, col));

			vector<RDStateMILPType> stateDiffVarsPerRound;

			stateDiffVarsPerRound.emplace_back(stateDiffVars);

			GRBLinExpr objExp = 0;

			int midR = (startr + endr) / 2;

			for (int i = startr; i < endr + 1; i++)
			{
				MilpPropagateThroughNonLinearLayer(model, stateDiffVars);
				MilpPropagateThroughLinearLayer(model, stateDiffVars, rotationalDifferencePerRoundConstant[i + 1]);

				if (i == midR)
					for (int row = 0; row < 5; row++)
						for (int col = 0; col < 64; col++)
							objExp += stateDiffVars[row][col];

				stateDiffVarsPerRound.emplace_back(stateDiffVars);
			}



			for (int row = 0; row < 5; row++)
				for (int col = 0; col < 64; col++)
					model.addConstr(stateDiffVars[row][col] == asconBase.GetBit(outputStateDiff, row, col));

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
					model.set(GRB_IntParam_SolutionNumber, solIndex);
					vector<StateType> stateDiffs;

					for (int i = 0; i <= r + 1; i++)
					{
						StateType stateDiff = MilpReadStateVars(stateDiffVarsPerRound[i]);
						stateDiffs.emplace_back(stateDiff);
					}

					difProb += ComputeDifferentialTrailProbability(stateDiffs, startr, endr);
				}

				return difProb;

			}
		}


	}
};

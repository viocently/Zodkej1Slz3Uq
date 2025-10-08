#pragma once
#include<iostream>
#include<array>
#include<vector>
#include"ascon.h"
#include"differential.h"
#include"gurobi_c++.h"

using namespace std;

using LINRowMILPType = array<GRBVar, 64>;
using LINStateMILPType = array<LINRowMILPType, 5>;


LINRowMILPType LINRotateLeft(const LINRowMILPType& value, int shift) {
	LINRowMILPType result;
	for (int i = 0; i < 64; i++) {
		result[i] = value[(i - shift + 64) % 64];
	}
	return result;
}

LINRowMILPType LINRotateRight(const LINRowMILPType& value, int shift) {
	return LINRotateLeft(value, 64 - shift);
}

int HammingWeight(uint64_t n) {
	int count = 0;
	while (n) {
		count += n & 1;
		n >>= 1;
	}
	return count;
}


class LinearBase
{
private:
	const AsconBase& asconBase;
	const RotationalDifferentialBase& rdBase;

	int rot;
	int milpThreads;
	array<RowType, 12> rotationalDifferencePerRoundConstant;
	array<array<int, 32>, 32> lat;
	array<array<double, 32>, 32> latCor2;
	array<array<double, 32>, 32> latCor2Log2;

	vector<vector<int>> inputMasksPerOutputMaskForSbox;

	vector<array<int, 11>> milpSboxInequalities;

	bool SwitchToNextPossibleInputMasks(vector<int>& activeColumnInputMaskIndexes, const vector<int>& activeColumnOutputMasks) const {
		for (int i = 0; i < activeColumnOutputMasks.size(); i++)
		{
			const int& activeColumnOutputMask = activeColumnOutputMasks[i];
			const vector<int>& possibleInputMasks = inputMasksPerOutputMaskForSbox[activeColumnOutputMask];
			int& activeColumnInputMaskIndex = activeColumnInputMaskIndexes[i];

			if (activeColumnInputMaskIndex < possibleInputMasks.size() - 1)
			{
				activeColumnInputMaskIndex++;
				return true;
			}
			else
			{
				activeColumnInputMaskIndex = 0;
			}

		}

		return false;
	}

	void GenerateLATForSbox()
	{
		for (int a = 0; a < 32; a++)
			for (int b = 0; b < 32; b++)
				lat[a][b] = 0;

		for (int a = 0; a < 32; a++)
			for (int b = 0; b < 32; b++)
			{
				for (uint8_t x = 0; x < 32; x++)
				{
					uint8_t y = asconBase.Sbox(x);
					uint8_t v0 = x & a;
					uint8_t v1 = y & b;
					if ((HammingWeight(v0) % 2) == (HammingWeight(v1) % 2))
						lat[a][b]++;
				}
			}

		for (int a = 0; a < 32; a++)
			for (int b = 0; b < 32; b++)
			{
				lat[a][b] = 2 * lat[a][b] - 32;
				latCor2[a][b] = lat[a][b] / 32.0;
				latCor2[a][b] = latCor2[a][b] * latCor2[a][b];
				latCor2Log2[a][b] = log2(latCor2[a][b]);
			}
	}

	void GenerateInputMasksPerOutputMaskForSbox()
	{
		inputMasksPerOutputMaskForSbox.resize(32);
		for (int a = 0; a < 32; a++)
		{
			for (int b = 0; b < 32; b++)
			{
				if (lat[a][b] != 0)
					inputMasksPerOutputMaskForSbox[b].push_back(a);
			}
		}
	}

public:
	LinearBase(const AsconBase& asconBase, const RotationalDifferentialBase& rdBase, int milpThreads = 2) : asconBase(asconBase), rdBase(rdBase), milpThreads(milpThreads)
	{
		this->rot = rdBase.GetRotationNumber();
		GenerateLATForSbox();
		GenerateInputMasksPerOutputMaskForSbox();

		milpSboxInequalities = {{0, -5, 4, 6, -3, -5, 4, 3, 3, 10, 10}, {0, 1, 1, 0, 0, 3, -2, 3, 2, 4, -2}, {0, 1, 0, 0, 0, 1, 1, 2, 1, -2, 1}, {4, -1, 0, 0, 0, -1, -1, 2, 1, -2, -1}, {0, -2, 2, -3, 4, -3, 7, 2, 1, 7, 4}, {0, 1, 0, 3, 0, 2, 1, 3, -3, 1, 1}, {2, 0, 0, -1, 0, 0, 0, -1, -1, 1, 1}, {0, 2, 2, 3, 3, 2, -1, 0, 0, -1, -1}, {4, 0, 0, -1, -1, 0, -1, 0, 0, -1, -1}, {4, -2, -2, 1, 0, 2, 0, -1, -1, 2, 0}, {8, -2, -4, -1, 0, -2, -2, 3, 4, 1, -1}, {2, 1, -1, 2, -2, 0, 1, -2, 1, 2, 3}, {0, 4, 2, 1, 0, -1, -1, 2, 4, 3, -2}, {3, -1, 0, -1, 0, 2, 1, 2, -1, -2, 1}, {5, 2, -2, -1, 0, -2, 1, -1, -1, 2, 0}, {0, 2, 5, 4, 3, 2, -3, -1, -1, 2, 3}, {6, 1, 0, -1, 0, -1, -1, 2, -1, -2, -2}, {3, -1, 0, 1, 1, -1, 0, 0, 0, -1, -1}, {0, -3, 3, -2, 0, 1, 3, 2, 5, 4, 1}, {2, 3, -3, -4, 3, 1, 6, -1, 2, 3, 4}, {8, -1, 2, -2, -2, -1, -1, -1, -1, -1, 0}, {8, 1, -3, 2, -2, -1, -2, -1, -2, 1, 2}, {0, 2, 2, 1, 0, 2, 0, -1, -1, 2, 0}, {7, 0, -2, 1, 2, -2, -2, -2, 2, -1, 1}, {3, -1, 1, 0, 0, -1, 0, -1, -1, 1, 0}, {2, -1, 0, 1, 0, 2, -1, 2, -1, 1, -1}, {5, -2, -2, 1, 1, 0, -1, -1, 2, 1, -1}, {8, -1, -2, -2, 2, -1, -1, -1, -1, -1, 2}, {9, -2, -2, -2, -2, 1, 1, -1, 1, -2, 1}, {0, 2, 0, 2, 0, 1, -1, 2, -1, 0, 2}, {3, -1, 1, 1, -1, 1, 1, -1, 0, -1, 0}, {4, -2, 2, -2, 2, 1, 1, -1, 1, -1, 2}, {4, 1, 0, 0, -1, -1, 1, -1, 1, -1, -1}, {0, -5, -5, 4, 4, 5, 10, 3, 2, 3, 10}, {0, 2, 0, 0, 1, 2, 2, 1, 2, -1, -1}, {0, 4, -2, 3, -2, 4, 0, 1, 1, 1, 3}, {2, 0, 0, -1, 0, -2, 2, 1, -1, 2, 1}, {3, 1, 1, -1, -1, 1, 0, -1, -1, 0, 0}, {4, 0, 1, 0, -1, -1, -1, -1, 1, -1, 1}, {0, 2, 3, 2, 3, 3, -1, -1, 1, -1, 1}, {0, 3, 1, -1, -1, 0, 3, 2, 3, -1, 3}, {2, 0, 0, 1, 1, 0, -1, 0, 0, -1, -1}, {3, -1, 0, 0, 1, -1, 1, 0, 1, -1, -1}, {3, -1, 0, 1, 0, -1, 0, 1, 0, -1, -1}, {3, -1, 1, 1, 1, -1, 0, -1, -1, 0, 0}, {3, 0, -1, 1, -1, 1, -1, -1, 0, 0, 1}, {4, -1, -1, 1, -1, 1, 1, -1, 0, 0, -1}, {4, -1, 0, -1, 0, -1, -1, 1, 0, -1, 0}, {4, -1, 0, -1, 1, 1, 1, 0, -1, -1, -1}, {4, 1, 0, -1, 1, -1, 1, 0, -1, -1, -1}, {4, 1, 0, 1, -1, -1, 1, -1, 0, -1, -1}, {0, -4, 2, 3, 0, -4, 2, 3, 3, 8, 6}, {0, -4, 5, 4, -3, -1, 2, 2, 3, 5, 6}, {0, -2, -3, 3, 3, -2, 5, 1, 6, 4, 5}, {0, -1, -1, -1, -1, 0, 3, 4, 4, 3, 3}, {0, -1, -1, 3, 0, -4, 5, 3, 4, 6, 2}, {0, 1, -1, -1, 0, 2, 1, 1, 2, 1, 0}, {0, 1, 1, 1, 1, -2, 1, 0, 2, 1, 1}, {0, 2, 0, 2, 1, 2, 1, 1, 0, -1, -1}, {0, 4, -1, 0, -1, 4, -1, 3, 4, -1, 4}, {3, -1, -1, -1, 0, 1, -1, 1, 1, 1, -1}, {3, 1, -1, -1, 1, 1, 0, -1, -1, 0, 1}, {3, 1, -1, 0, 0, -1, 0, -1, -1, 1, 0}, {3, 1, -1, 1, 1, -1, 1, -1, -1, 0, 1}, {3, 1, 1, -2, 1, -1, 1, -1, -1, 1, 2}, {3, 1, 1, 1, -2, -1, 2, -1, -1, 1, 2}, {4, -1, 0, 0, -1, 1, 1, -1, 1, -1, -1}, {4, 1, -1, -1, -1, -1, 1, -1, -1, 1, 1}, {5, -1, -1, 1, -2, -1, 0, -1, -1, 1, 2}, {5, 1, 0, -1, -1, 1, 0, -1, -1, -1, -1}, {7, -1, 0, -1, -1, -1, 0, -1, -1, -1, -1}};

	}

	array<array<double, 32>, 32> GetLATCor2() const
	{
		return latCor2;
	}

	array<array<double, 32>, 32> GetLATCor2Log2() const
	{
		return latCor2Log2;
	}

	vector<vector<int>> GetInputMasksPerOutputMaskForSbox() const
	{
		return inputMasksPerOutputMaskForSbox;
	}

	int GetRotationNumber() const
	{
		return rot;
	}

	int GetActiveColumnsNum(const StateType& stateMask) const
	{
		int activeCols = 0;
		for (int col = 0; col < 64; col++)
		{
			if (asconBase.GetColumn(stateMask, col) != 0)
				activeCols++;
		}
		return activeCols;
	}

	void PropagateThroughInvRotateRows(StateType& inputStateMask) const
	{
		asconBase.RotateRowsTranspose(inputStateMask);
	}

	void PropagateThroughAddRoundConstant(const StateType& inputStateMask, double& cor2, const RowType rcDiff) const
	{
		if (HammingWeight(rcDiff & inputStateMask[2]) % 2 == 1)
			cor2 = -cor2;
	}

	void PropagateThroughInvLinearLayer(StateType& inputStateMask) const
	{
		PropagateThroughInvRotateRows(inputStateMask);
	}

	void PropagateThroughInvLinearLayerWithCor2(StateType& inputStateMask, double& cor2, const RowType rcDiff) const
	{
		PropagateThroughAddRoundConstant(inputStateMask, cor2, rcDiff);
		PropagateThroughInvRotateRows(inputStateMask);
	}

	double ComputeApproximationCorrelationSquareManually(const StateType& inputStateMask, const StateType& outputStateMask, int startr, int endr) const
	{
		int r = endr - startr;


		if (r == 0)
		{
			double approxCor2 = 1.0;

			StateType sboxesInputMask = inputStateMask;
			// asconBase.OutputState(cout, sboxesInputMask);
			StateType sboxesOutputMask = outputStateMask;
			// asconBase.OutputState(cout, sboxesOutputMask);
			PropagateThroughInvLinearLayerWithCor2(sboxesOutputMask, approxCor2, rotationalDifferencePerRoundConstant[startr + 1]);

			for (int col = 0; col < 64; col++)
			{
				uint8_t inputColumnMask = asconBase.GetColumn(sboxesInputMask, col);
				uint8_t outputColumnMask = asconBase.GetColumn(sboxesOutputMask, col);
				approxCor2 *= latCor2[inputColumnMask][outputColumnMask];
			}

			return approxCor2;
		}
		else
		{
			double approxCor2 = 0.0;
			StateType prevOutputStateMaskBase = { 0 };
			vector<int> activeColumns;
			vector<int> activeColumnInputMaskIndexes;
			vector<int> activeColumnOutputMasks;
			StateType sboxesInputMask = inputStateMask;
			StateType sboxesOutputMask = outputStateMask;
			PropagateThroughInvLinearLayer(sboxesOutputMask);

			double sign = 1.0;
			PropagateThroughAddRoundConstant(outputStateMask, sign, rotationalDifferencePerRoundConstant[endr + 1]);



			for (int col = 0; col < 64; col++)
			{
				uint8_t outputColumnMask = asconBase.GetColumn(sboxesOutputMask, col);
				if (outputColumnMask == 0)
					asconBase.SetColumn(prevOutputStateMaskBase, col, 0);
				else
				{
					activeColumns.emplace_back(col);
					activeColumnInputMaskIndexes.emplace_back(0);
					activeColumnOutputMasks.emplace_back(outputColumnMask);
				}
			}

			int activeColumnsNum = activeColumns.size();

			bool isAllActiveColumnsProcessed = false;
			while (!isAllActiveColumnsProcessed)
			{
				StateType prevOutputStateMask = prevOutputStateMaskBase;

				double cor2CurNonLinearLayer = sign;

				for (int activeColumnIndex = 0; activeColumnIndex < activeColumnsNum; activeColumnIndex++)
				{
					const int& activeColumnPos = activeColumns[activeColumnIndex];
					const int& activeColumnInputMaskIndex = activeColumnInputMaskIndexes[activeColumnIndex];
					const int& activeColumnOutputMask = activeColumnOutputMasks[activeColumnIndex];
					const int& activeColumnInputMask = inputMasksPerOutputMaskForSbox[activeColumnOutputMask][activeColumnInputMaskIndex];

					asconBase.SetColumn(prevOutputStateMask, activeColumnPos, activeColumnInputMask);

					cor2CurNonLinearLayer *= latCor2[activeColumnInputMask][activeColumnOutputMask];
				}

				double prevLinCor2 = ComputeApproximationCorrelationSquareManually(sboxesInputMask, prevOutputStateMask, startr, endr - 1);

				if (prevLinCor2 != 0.0)
					approxCor2 += prevLinCor2 * cor2CurNonLinearLayer;

				if (SwitchToNextPossibleInputMasks(activeColumnInputMaskIndexes, activeColumnOutputMasks))
					isAllActiveColumnsProcessed = false;
				else
					isAllActiveColumnsProcessed = true;

			}

			return approxCor2;
		}
	}

	double ComputeLinearTrailCorrelationSquare(vector<StateType>& stateMasks, int startr, int endr) const
	{
		int r = endr - startr;


		double trailCor2 = 1.0;
		for (int i = 0; i < r + 1; i++)
		{
			StateType curInputStateMask = stateMasks[i];
			StateType curOutputStateDMask = stateMasks[i + 1];

			double curLinCor2 = ComputeApproximationCorrelationSquareManually(curInputStateMask, curOutputStateDMask, startr + i, startr + i);

			trailCor2 *= curLinCor2;
		}

		return trailCor2;
	}



	// =====================================================================================
	// ================= Functions below are for the milp models of Linear propagation through xoodoo ==================
	// =====================================================================================

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

	LINStateMILPType MilpAddStateVars(GRBModel& model) const
	{
		LINStateMILPType stateVars;
		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
				stateVars[row][col] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		return stateVars;
	}

	StateType MilpReadStateVars(const LINStateMILPType& stateVars) const
	{
		StateType state;
		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
				asconBase.SetBit(state, row, col, round(stateVars[row][col].get(GRB_DoubleAttr_Xn)));
		return state;
	}

	void MilpPropagateThroughInvRotateRows(GRBModel & model, LINStateMILPType& stateVars) const
	{
		for (int row = 0; row < 5; row++)
		{
			int rot0 = asconBase.GetRotationNum(row, 0);
			int rot1 = asconBase.GetRotationNum(row, 1);
			auto tmp0 = LINRotateLeft(stateVars[row], rot0);
			auto tmp1 = LINRotateLeft(stateVars[row], rot1);
			for (int col = 0; col < 64; col++)
			{
				GRBVar newVar = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				MilpXor3(model, stateVars[row][col], tmp0[col], tmp1[col], newVar);
				stateVars[row][col] = newVar;
			}
		}
	}

	void MilpPropagateThroughInvSbox(GRBModel& model, GRBVar& u0, GRBVar& u1, GRBVar& u2, GRBVar& u3, GRBVar & u4) const
	{
		GRBVar v0 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		GRBVar v1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		GRBVar v2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		GRBVar v3 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		GRBVar v4 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

		for (auto& ieq : milpSboxInequalities)
		{
			GRBLinExpr ieqExp = ieq[0];
			ieqExp += v0 * ieq[1] + v1 * ieq[2] + v2 * ieq[3] + v3 * ieq[4] + v4 * ieq[5];
			ieqExp += u0 * ieq[6] + u1 * ieq[7] + u2 * ieq[8] + u3 * ieq[9] + u4 * ieq[10];
			model.addConstr(ieqExp >= 0);
		}

		u0 = v0;
		u1 = v1;
		u2 = v2;
		u3 = v3;
		u4 = v4;
	}

	void MilpPropagateThroughInvSboxes(GRBModel& model, RDStateMILPType& stateVars) const
	{
		for (int col = 0; col < 64; col++)
		{
			MilpPropagateThroughInvSbox(model, stateVars[4][col], stateVars[3][col], stateVars[2][col], stateVars[1][col], stateVars[0][col]);
		}
	}

	void MilpPropagateThroughInvLinearLayer(GRBModel& model, RDStateMILPType& stateVars) const
	{
		MilpPropagateThroughInvRotateRows(model, stateVars);
	}

	void MilpPropagateThroughInvNonLinearLayer(GRBModel& model, RDStateMILPType& stateVars) const
	{
		MilpPropagateThroughInvSboxes(model, stateVars);
	}

	double ComputeApproximationCorrelationSquareAutomatically(const StateType& inputStateMask, const StateType& outputStateMask, int startr, int endr) const
	{
		int r = endr - startr;

		if (r == 0)
			return ComputeApproximationCorrelationSquareManually(inputStateMask, outputStateMask, startr, endr);
		else
		{
			GRBEnv env = GRBEnv();

			env.set(GRB_IntParam_LogToConsole, 0);
			env.set(GRB_IntParam_PoolSearchMode, 2);
			env.set(GRB_IntParam_PoolSolutions, 200000000);
			env.set(GRB_IntParam_Threads, milpThreads);

			GRBModel model = GRBModel(env);

			LINStateMILPType stateMaskVars = MilpAddStateVars(model);

			for (int row = 0; row < 5; row++)
				for (int col = 0; col < 64; col++)
					model.addConstr(stateMaskVars[row][col] == asconBase.GetBit(outputStateMask, row, col));

			vector<RDStateMILPType> stateMaskVarsPerRound;

			stateMaskVarsPerRound.emplace_back(stateMaskVars);

			GRBLinExpr objExp = 0;

			int midR = (startr + endr) / 2;

			for (int i = endr; i >= startr; i--)
			{

				MilpPropagateThroughInvLinearLayer(model, stateMaskVars);
				MilpPropagateThroughInvNonLinearLayer(model, stateMaskVars);

				if (i == midR)
					for (int row = 0; row < 5; row++)
						for (int col = 0; col < 64; col++)
							objExp += stateMaskVars[row][col];

				stateMaskVarsPerRound.insert(stateMaskVarsPerRound.begin(), stateMaskVars);
			}



			for (int row = 0; row < 5; row++)
				for (int col = 0; col < 64; col++)
					model.addConstr(stateMaskVars[row][col] == asconBase.GetBit(inputStateMask, row, col));

			model.setObjective(objExp, GRB_MINIMIZE);

			model.optimize();

			int solCount = model.get(GRB_IntAttr_SolCount);
			cout << "Get " << solCount << " linear trails." << endl;

			if (solCount == 0)
				return 0.0;
			else
			{
				double linCor2 = 0.0;
				for (int solIndex = 0; solIndex < solCount; solIndex++)
				{
					model.set(GRB_IntParam_SolutionNumber, solIndex);
					vector<StateType> stateMasks;

					for (int i = 0; i <= r + 1; i++)
					{
						StateType stateMask = MilpReadStateVars(stateMaskVarsPerRound[i]);
						stateMasks.emplace_back(stateMask);
					}

					linCor2 += ComputeLinearTrailCorrelationSquare(stateMasks, startr, endr);
				}

				return linCor2;

			}
		}


	}
};

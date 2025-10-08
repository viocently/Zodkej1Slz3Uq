#pragma once
#include<iostream>
#include<array>
#include<vector>
#include"xoodoo.h"
#include"differential.h"
#include"gurobi_c++.h"

using namespace std;

using LINLaneMILPType = array<GRBVar, 32>;
using LINSheetMILPType = array<LINLaneMILPType, 3>;
using LINPlaneMILPType = array<LINLaneMILPType, 4>;
using LINStateMILPTypeBySheet = array<LINSheetMILPType, 4>;

LINLaneMILPType LINRotateLeft(LINLaneMILPType laneVars, int shift) {
	LINLaneMILPType newLaneVars;
	for (int i = 0; i < 32; i++)
		newLaneVars[i] = laneVars[(i - shift + 32) % 32];
	return newLaneVars;
}

// Compute hamming weight of a 32-bit unsigned integer
int HammingWeight(uint32_t n)  {
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
	const XoodooBase& xoodooBase;
	const RotationalDifferentialBase& rdBase;

	int rot;
	int milpThreads;
	array<LaneType, 12> rotationalDifferencePerRoundConstants;
	array<array<int, 8>, 8> lat;
	array<array<double, 8>, 8> latCor2;
	array<array<double, 8>, 8> latCor2Log2;
	vector<vector<int>> inputMasksPerOutputMaskForChi;

	vector<array<int, 7>> milpChiInequalities;



	bool SwitchToNextPossibleInputMasks(vector<int>& activeColumnInputMaskIndexes, const vector<int>& activeColumnOutputMasks) const {
		for (int i = 0; i < activeColumnOutputMasks.size(); i++)
		{
			const int& activeColumnOutputMask = activeColumnOutputMasks[i];
			const vector<int>& possibleInputMasks = inputMasksPerOutputMaskForChi[activeColumnOutputMask];
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


	// Generate the square of LAT for the Chi operation in the round function of Xoodoo
	void GenerateLATForChi() {
		for (int a = 0; a < 8; a++)
			for (int b = 0; b < 8; b++)
				lat[a][b] = 0;

		for (int a = 0; a < 8; a++)
		{
			for (int b = 0; b < 8; b++)
			{
				for (uint8_t x = 0; x < 8; x++)
				{
					uint8_t y = xoodooBase.ChiColumnByOperation(x);
					uint8_t v0 = x & a;
					uint8_t v1 = y & b;
					if ((HammingWeight(v0) + HammingWeight(v1)) % 2 == 0)
						lat[a][b]++;
				}
			}
		}

		for (int a = 0; a < 8; a++)
			for (int b = 0; b < 8; b++)
			{
				lat[a][b] = 2 * lat[a][b] - 8;
				latCor2[a][b] = lat[a][b] / 8.0;
				latCor2[a][b] = latCor2[a][b] * latCor2[a][b];
				latCor2Log2[a][b] = log2(latCor2[a][b]);
			}
	}

	// Generate possible input masks for each output masks of Chi operation
	void GenerateInputMasksPerOutputMaskForChi() {
		inputMasksPerOutputMaskForChi.resize(8);
		for (int a = 0; a < 8; a++)
		{
			for (int b = 0; b < 8; b++)
				if (lat[a][b] != 0)
					inputMasksPerOutputMaskForChi[b].emplace_back(a);
		}
	}

	



public:
	LinearBase(const XoodooBase& xoodooBase, const RotationalDifferentialBase& rdBase, int milpThreads = 2) : xoodooBase(xoodooBase), rdBase(rdBase), milpThreads(milpThreads) {
		rot = rdBase.getRotationNumber();
		rotationalDifferencePerRoundConstants = rdBase.getRotationalDifferencePerRoundConstants();

		GenerateLATForChi();
		// print LAT in a 8 x 8 table with row index being the input linear mask and the column index being the output linear mask

		for (int i = 0; i < 8; i++)
			cout << "\t" << i;
		cout << endl;
		for (int i = 0; i < 8; i++)
		{
			cout << i << "\t";
			for (int j = 0; j < 8; j++)
				cout << lat[i][j] << "\t";
			cout << endl;
		}


		GenerateInputMasksPerOutputMaskForChi();
		milpChiInequalities = { {0, -1, 0, 0, 1, 1, 1},{0, 0, 0, 1, 1, 1, -1},{0, 0, 1, 0, 1, -1, 1},{0, 1, -1, 0, 0, 1, 1},{0, 1, 1, 0, -1, 0, 1},{0, 1, 1, 1, 0, 0, -1},{3, -1, -1, 0, -1, -1, 1},{3, -1, 0, -1, -1, 1, -1},{3, 0, -1, -1, 1, -1, -1},{0, -1, 1, 1, 1, 0, 0},{0, 0, 0, -1, 1, 1, 1},{0, 1, -1, 1, 0, 1, 0},{3, -1, -1, 1, -1, -1, 0},{3, -1, 1, -1, -1, 0, -1},{3, 1, -1, -1, 0, -1, -1}
		};
	}

	array<array<double, 8>, 8> getLATCor2() const {
		return latCor2;
	}

	array<array<double, 8>, 8> getLATCor2Log2() const {
		return latCor2Log2;
	}

	vector<vector<int>> getInputMasksPerOutputMaskForChi() const {
		return inputMasksPerOutputMaskForChi;
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

	void PropagateThroughRevLinearLayer(StateTypeBySheet& inputStateMask, const LaneType rcDiff) const {
		xoodooBase.RhoWestTranspose(inputStateMask);
		xoodooBase.ThetaTranspose(inputStateMask);
		xoodooBase.RhoEastTranspose(inputStateMask);
	}

	void PropagateThroughLota(const StateTypeBySheet& inputStateMask, double& cor2, const LaneType rcDiff) const {
		LaneType lane00 = inputStateMask[0][0];
		// cout << "Diff " << hex << rcDiff << endl;
		// cout << "lane00 " << hex << lane00 << endl;
		if (HammingWeight(rcDiff & lane00) % 2 == 1)
		{
			// cout << "flipped" << endl;
			cor2 = -cor2;
		}
	}


	void PropagateThroughRevLinearLayerWithCor2(StateTypeBySheet& inputStateMask, double& cor2, const LaneType rcDiff) const
	{
		PropagateThroughLota(inputStateMask, cor2, rcDiff);
		PropagateThroughRevLinearLayer(inputStateMask, rcDiff);
	}

	
	double ComputeApproximationCorrelationSquareManually(const StateTypeBySheet& inputStateMask, const StateTypeBySheet& outputStateMask, int startr, int endr, int totalRounds) const {
		int r = endr - startr;

		if (r < 0 || startr < 0 || endr < 0 || startr > totalRounds || endr > totalRounds)
		{
			throw MyException(__func__ + string(": Parameters error"));
		}
		else if (r == 0)
		{
			double approxCor2 = 1.0;
			for(int x = 0; x < 4; x++)
				for (int z = 0; z < 32; z++)
				{
					uint8_t inputColumnMask = xoodooBase.GetColumn(inputStateMask, x, z);
					uint8_t outputColumnMask = xoodooBase.GetColumn(outputStateMask, x, z);
					approxCor2 *= latCor2[inputColumnMask][outputColumnMask];
				}

			return approxCor2;
		}
		else
		{
			double approxCor2 = 0.0;

			StateTypeBySheet prevOutputStateMaskBase = { 0 };
			vector<pair<int, int>> activeColumns;
			vector<int> activeColumnInputMaskIndexes;
			vector<int> activeColumnOutputMasks;

			int rotationalDifferenceIndexPrevRoundConstant = xoodooBase.ComputeRoundConstantIndex(totalRounds, endr);

			for(int x = 0; x < 4;x++)
				for (int z = 0; z < 32; z++)
				{
					uint8_t outputColumnMask = xoodooBase.GetColumn(outputStateMask, x, z);
					if (outputColumnMask == 0)
						xoodooBase.SetColumn(prevOutputStateMaskBase, x, z, 0);
					else
					{
						activeColumns.emplace_back(x, z);
						activeColumnInputMaskIndexes.emplace_back(0);
						activeColumnOutputMasks.emplace_back(outputColumnMask);
					}
				}

			int activeColumnsNum = activeColumns.size();

			bool isAllActiveColumnsProcessed = false;
			while (!isAllActiveColumnsProcessed)
			{
				StateTypeBySheet prevOutputStateMask = prevOutputStateMaskBase;

				double cor2CurNonLinearLayer = 1.0;

				for (int activeColumnIndex = 0; activeColumnIndex < activeColumnsNum; activeColumnIndex++)
				{
					const pair<int, int>& activeColumnPos = activeColumns[activeColumnIndex];
					const int & activeColumnInputMaskIndex = activeColumnInputMaskIndexes[activeColumnIndex];
					const int & activeColumnOutputMask = activeColumnOutputMasks[activeColumnIndex];
					const int & activeColumnInputMask = inputMasksPerOutputMaskForChi[activeColumnOutputMask][activeColumnInputMaskIndex];

					xoodooBase.SetColumn(prevOutputStateMask, activeColumnPos.first, activeColumnPos.second, activeColumnInputMask);

					cor2CurNonLinearLayer *= latCor2[activeColumnInputMask][activeColumnOutputMask];
				}

				PropagateThroughRevLinearLayerWithCor2(prevOutputStateMask, cor2CurNonLinearLayer, rotationalDifferencePerRoundConstants[rotationalDifferenceIndexPrevRoundConstant]);
				double prevLinCor2 = ComputeApproximationCorrelationSquareManually(inputStateMask, prevOutputStateMask, startr, endr - 1, totalRounds);

				if (prevLinCor2 != 0.0)
					approxCor2 += prevLinCor2 * cor2CurNonLinearLayer;

				// switch to next possible input masks for the active columns
				if (SwitchToNextPossibleInputMasks(activeColumnInputMaskIndexes, activeColumnOutputMasks))
					isAllActiveColumnsProcessed = false;
				else
					isAllActiveColumnsProcessed = true;
				
			}

			return approxCor2;
		}
	}


	double ComputeLinearTrailCorrelation2(const vector<StateTypeBySheet>& inputStateMasks, const vector<StateTypeBySheet>& outputStateMasks, int startr, int endr, int totalRounds) const
	{
		int r = endr - startr;

		double trailCor2 = 1.0;
		if (r < 0 || startr < 0 || endr < 0 || startr > totalRounds || endr > totalRounds || inputStateMasks.size() != endr - startr + 1 || outputStateMasks.size() != endr - startr + 1)
		{
			throw MyException(__func__ + string(": Parameters error"));
		}
		else
		{
			for (int i = 0; i < r; i++)
			{
				StateTypeBySheet curInputStateMask = inputStateMasks[i];
				StateTypeBySheet curOutputStateMask = outputStateMasks[i];

				double cor2CurNonLinearLayer = ComputeApproximationCorrelationSquareManually(curInputStateMask, curOutputStateMask, startr + i, startr + i, totalRounds);

				if (i > 0)
				{
					PropagateThroughLota(curInputStateMask, cor2CurNonLinearLayer, rotationalDifferencePerRoundConstants[xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + i)]);
				}
				

				trailCor2 *= cor2CurNonLinearLayer;
			}

			double cor2LastNonLinearLayer = ComputeApproximationCorrelationSquareManually(inputStateMasks[r], outputStateMasks[r], startr + r, endr, totalRounds);

			if(r > 0)
				PropagateThroughLota(inputStateMasks[r], cor2LastNonLinearLayer, rotationalDifferencePerRoundConstants[xoodooBase.ComputeRoundConstantIndex(totalRounds, endr)]);

			trailCor2 *= cor2LastNonLinearLayer;
			


			return trailCor2;
		}
	}


	bool VerifyLinearTrail(const vector<StateTypeBySheet>& inputStateMasks, const vector<StateTypeBySheet>& outputStateMasks, int startr, int endr, int totalRounds, double toVerifyTrailCor2) const
	{
		int r = endr - startr;

		if (r < 0 || startr < 0 || endr < 0 || startr > totalRounds || endr > totalRounds || inputStateMasks.size() != endr - startr + 1 || outputStateMasks.size() != endr - startr + 1)
		{
			throw MyException(__func__ + string(": Parameters error"));
		}
		else
		{
			for (int i = 0; i < r; i++)
			{
				StateTypeBySheet curInputStateMask = inputStateMasks[i];
				StateTypeBySheet curOutputStateMask = inputStateMasks[i];
				StateTypeBySheet nextInputStateMask = inputStateMasks[i + 1];

				int curRotationalDifferenceIndexNextRoundConstant = xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + i + 1);

				PropagateThroughRevLinearLayer(nextInputStateMask, rotationalDifferencePerRoundConstants[curRotationalDifferenceIndexNextRoundConstant]);

				if (curOutputStateMask != nextInputStateMask)
					return false;
			}

			double trailCor2 = ComputeLinearTrailCorrelation2(inputStateMasks, outputStateMasks, startr, endr, totalRounds);

			if(abs(trailCor2 - toVerifyTrailCor2) < 1e-2)
				return true;
			else
				return false;
		}
	}



// =====================================================================================
// ================= Functions below are for the milp models of Linear propagation through xoodoo ==================
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


	LINStateMILPTypeBySheet MilpAddStateVars(GRBModel& model) const
	{
		LINStateMILPTypeBySheet stateMaskVars;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					stateMaskVars[x][y][z] = model.addVar(0, 1, 0, GRB_BINARY);
		return stateMaskVars;
	}

	StateTypeBySheet MilpReadStateVars(LINStateMILPTypeBySheet& stateMaskVars) const
	{
		StateTypeBySheet stateMask;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					xoodooBase.SetBit(stateMask, x, y, z, round(stateMaskVars[x][y][z].get(GRB_DoubleAttr_Xn)));
		return stateMask;
	}

	LINPlaneMILPType MilpAddPlaneVars(GRBModel& model) const
	{
		LINPlaneMILPType planeMaskVars;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 32; z++)
				planeMaskVars[x][z] = model.addVar(0, 1, 0, GRB_BINARY);
		return planeMaskVars;
	}

	void MilpPropagateThroughRevRhoEast(LINStateMILPTypeBySheet& stateMaskVars) const
	{
		LINStateMILPTypeBySheet tmpStateMaskVars = stateMaskVars;
		for(int x = 0; x < 4; x++)
			stateMaskVars[x][1] = LINRotateLeft(tmpStateMaskVars[x][1], -1);

		for (int x = 0; x < 4; x++)
			stateMaskVars[x][2] = LINRotateLeft(tmpStateMaskVars[(x + 2) % 4][2], -8);
	}

	void MilpPropagateThroughRevTheta(GRBModel& model, LINStateMILPTypeBySheet& stateMaskVars) const
	{
		LINPlaneMILPType parityMaskVars = MilpAddPlaneVars(model);
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 32; z++)
				MilpXor3(model, stateMaskVars[x][0][z], stateMaskVars[x][1][z], stateMaskVars[x][2][z], parityMaskVars[x][z]);

		LINStateMILPTypeBySheet newStateMaskVars = MilpAddStateVars(model);

		for (int x = 0; x < 4; x++)
		{
			LINLaneMILPType one = LINRotateLeft(parityMaskVars[(x + 1) % 4], -5);
			LINLaneMILPType two = LINRotateLeft(parityMaskVars[(x + 1) % 4], -14);
			for(int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					MilpXor3(model, stateMaskVars[x][y][z], one[z], two[z], newStateMaskVars[x][y][z]);
		}

		stateMaskVars = newStateMaskVars;
	}

	void MilpPropagateThroughRevRhoWest(LINStateMILPTypeBySheet& stateMaskVars) const
	{
		LINStateMILPTypeBySheet tmpStateMaskVars = stateMaskVars;
		for (int x = 0; x < 4; x++)
			stateMaskVars[x][1] = tmpStateMaskVars[(x + 1) % 4][1];

		for (int x = 0; x < 4; x++)
			stateMaskVars[x][2] = LINRotateLeft(tmpStateMaskVars[x][2], -11);
	}

	void MilpPropagateThroughRevChiColumn(GRBModel& model, GRBVar& u0, GRBVar& u1, GRBVar& u2) const
	{
		GRBVar v0 = model.addVar(0, 1, 0, GRB_BINARY);
		GRBVar v1 = model.addVar(0, 1, 0, GRB_BINARY);
		GRBVar v2 = model.addVar(0, 1, 0, GRB_BINARY);

		for (auto& ieq : milpChiInequalities)
		{
			GRBLinExpr ieqExp = ieq[0];
			ieqExp += v0 * ieq[1] + v1 * ieq[2] + v2 * ieq[3];
			ieqExp += u0 * ieq[4] + u1 * ieq[5] + u2 * ieq[6];
			model.addConstr(ieqExp >= 0);
		}

		u0 = v0;
		u1 = v1;
		u2 = v2;
	}

	void MilpPropagateThroughRevChi(GRBModel& model, LINStateMILPTypeBySheet& stateMaskVars) const
	{
		for (int x = 0; x < 4; x++)
			for(int z = 0; z < 32; z++)
				MilpPropagateThroughRevChiColumn(model, stateMaskVars[x][0][z], stateMaskVars[x][1][z], stateMaskVars[x][2][z]);
	}
	
	void MilpPropagateThroughRevLinearLayer(GRBModel& model, LINStateMILPTypeBySheet& stateMaskVars) const
	{
		MilpPropagateThroughRevRhoWest(stateMaskVars);
		MilpPropagateThroughRevTheta(model, stateMaskVars);
		MilpPropagateThroughRevRhoEast(stateMaskVars);
	}

	void MilpPropagateThroughRevNonLinearLayer(GRBModel& model, LINStateMILPTypeBySheet & stateMaskVars) const
	{
		MilpPropagateThroughRevChi(model, stateMaskVars);
	}

	double ComputeApproximationCorrelationSquareAutomatically(const StateTypeBySheet& inputStateMask, const StateTypeBySheet& outputStateMask, int startr, int endr, int totalRounds) const
	{
		int r = endr - startr;

		if (r < 0 || startr < 0 || endr < 0 || startr > totalRounds || endr > totalRounds)
		{
			throw MyException(__func__ + string(": Parameters error"));
		}
		else if (r == 0)
		{
			double approxCor2 = 1.0;
			for (int x = 0; x < 4; x++)
				for (int z = 0; z < 32; z++)
				{
					uint8_t inputColumnMask = xoodooBase.GetColumn(inputStateMask, x, z);
					uint8_t outputColumnMask = xoodooBase.GetColumn(outputStateMask, x, z);
					approxCor2 *= latCor2[inputColumnMask][outputColumnMask];
				}

			return approxCor2;
		}
		else
		{

			GRBEnv env = GRBEnv();
			env.set(GRB_IntParam_LogToConsole, 0);
			env.set(GRB_IntParam_PoolSearchMode, 2);
			env.set(GRB_IntParam_PoolSolutions, 200000000);
			env.set(GRB_IntParam_Threads, milpThreads);

			GRBModel model = GRBModel(env);

			LINStateMILPTypeBySheet stateMaskVars = MilpAddStateVars(model);

			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					for (int z = 0; z < 32; z++)
						model.addConstr(stateMaskVars[x][y][z] == xoodooBase.GetBit(outputStateMask, x, y, z));

			vector<LINStateMILPTypeBySheet> inputStateMaskVarsPerRoundChi;
			vector<LINStateMILPTypeBySheet> outputStateMaskVarsPerRoundChi;


			int midR = (startr + endr + 1) / 2;
			GRBLinExpr objExp = 0;

			for (int i = endr; i >= startr + 1; i--)
			{
				outputStateMaskVarsPerRoundChi.insert(outputStateMaskVarsPerRoundChi.begin(), stateMaskVars);
				MilpPropagateThroughRevChi(model, stateMaskVars);
				inputStateMaskVarsPerRoundChi.insert(inputStateMaskVarsPerRoundChi.begin(), stateMaskVars);

				if(i == midR)
				{
					for (int x = 0; x < 4; x++)
						for (int y = 0; y < 3; y++)
							for (int z = 0; z < 32; z++)
								objExp += stateMaskVars[x][y][z];
				}

				MilpPropagateThroughRevLinearLayer(model, stateMaskVars);
			}

			outputStateMaskVarsPerRoundChi.insert(outputStateMaskVarsPerRoundChi.begin(), stateMaskVars);
			MilpPropagateThroughRevChi(model, stateMaskVars);
			inputStateMaskVarsPerRoundChi.insert(inputStateMaskVarsPerRoundChi.begin(), stateMaskVars);

			for(int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					for (int z = 0; z < 32; z++)
						model.addConstr(stateMaskVars[x][y][z] == xoodooBase.GetBit(inputStateMask, x, y, z));

			model.setObjective(objExp, GRB_MINIMIZE);

			model.optimize();

			int solCount = model.get(GRB_IntAttr_SolCount);
			cout << "Get " << solCount << " linear trails." << endl;

			if (solCount == 0)
				return 0.0;
			else
			{
				double approxCor2 = 0.0;
				for (int solIndex = 0; solIndex < solCount; solIndex++)
				{

					model.set(GRB_IntParam_SolutionNumber, solIndex);
					vector<StateTypeBySheet> inputStateMasks, outputStateMasks;
					for (int i = 0; i <= r; i++)
					{
						StateTypeBySheet inputStateMask = MilpReadStateVars(inputStateMaskVarsPerRoundChi[i]);
						StateTypeBySheet outputStateMask = MilpReadStateVars(outputStateMaskVarsPerRoundChi[i]);
						inputStateMasks.push_back(inputStateMask);
						outputStateMasks.push_back(outputStateMask);
						// xoodooBase.OutputState(inputStateMask);
						// xoodooBase.OutputState(outputStateMask);
					}


					approxCor2 += ComputeLinearTrailCorrelation2(inputStateMasks, outputStateMasks, startr, endr, totalRounds);
				}

				return approxCor2;
			}
		}
	}

};
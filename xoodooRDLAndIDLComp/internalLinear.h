
#pragma once
#include<iostream>
#include<array>
#include<vector>
#include"xoodoo.h"
#include"internalDifferential.h"
#include"linear.h"
#include"gurobi_c++.h"

using namespace std;

using ILINLaneMILPType = array<GRBVar, 16>;
using ILINSheetMILPType = array<ILINLaneMILPType, 3>;
using ILINPlaneMILPType = array<ILINLaneMILPType, 4>;
using ILINStateMILPTypeBySheet = array<ILINSheetMILPType, 4>;

ILINLaneMILPType LINRotateLeft(ILINLaneMILPType laneVars, int shift) {
	ILINLaneMILPType newLaneVars;
	for (int i = 0; i < 16; i++)
		newLaneVars[i] = laneVars[(i - shift + 16) % 16];
	return newLaneVars;
}

class InternalLinearBase
{
private:
	const XoodooBase& xoodooBase;
	const InternalDifferentialBase& idBase;
	const LinearBase& linBase;

	int offset;
	int milpThreads;
	array<HalfLaneType, 12> internalDifferencePerRoundConstants;
	array<array<double, 8>, 8> latCor2;
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


public:
	InternalLinearBase(const XoodooBase& xoodooBase, const InternalDifferentialBase& idBase, const LinearBase& linBase, int offset = 0, int milpThreads = 2) :
		xoodooBase(xoodooBase), idBase(idBase), linBase(linBase), offset(offset), milpThreads(milpThreads)
	{
		internalDifferencePerRoundConstants = idBase.getInternalDifferencePerRoundConstants();
		latCor2 = linBase.getLATCor2();
		inputMasksPerOutputMaskForChi = linBase.getInputMasksPerOutputMaskForChi();
		milpChiInequalities = { {0, -1, 0, 0, 1, 1, 1},{0, 0, 0, 1, 1, 1, -1},{0, 0, 1, 0, 1, -1, 1},{0, 1, -1, 0, 0, 1, 1},{0, 1, 1, 0, -1, 0, 1},{0, 1, 1, 1, 0, 0, -1},{3, -1, -1, 0, -1, -1, 1},{3, -1, 0, -1, -1, 1, -1},{3, 0, -1, -1, 1, -1, -1},{0, -1, 1, 1, 1, 0, 0},{0, 0, 0, -1, 1, 1, 1},{0, 1, -1, 1, 0, 1, 0},{3, -1, -1, 1, -1, -1, 0},{3, -1, 1, -1, -1, 0, -1},{3, 1, -1, -1, 0, -1, -1}
		};
	}

	array<array<double, 8>, 8> getLATCor2() const {
		return latCor2;
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

	
	void PropagateThroughRevLinearLayer(HalfStateTypeBySheet& inputStateMask, const HalfLaneType rcDiff) const {
		xoodooBase.RhoWestTranspose(inputStateMask);
		xoodooBase.ThetaTranspose(inputStateMask);
		xoodooBase.RhoEastTranspose(inputStateMask);
	}

	void PropagateThroughLota(const HalfStateTypeBySheet& inputStateMask, double& cor2, const HalfLaneType rcDiff) const {
		HalfLaneType lane00 = inputStateMask[0][0];
		// cout << "Diff " << hex << rcDiff << endl;
		// cout << "lane00 " << hex << lane00 << endl;
		if (HammingWeight(rcDiff & lane00) % 2 == 1)
		{
			// cout << "flipped" << endl;
			cor2 = -cor2;
		}
	}

	void PropagateThroughRevLinearLayerWithCor2(HalfStateTypeBySheet& inputStateMask, double& cor2, const HalfLaneType rcDiff) const
	{
		PropagateThroughLota(inputStateMask, cor2, rcDiff);
		PropagateThroughRevLinearLayer(inputStateMask, rcDiff);
	}

	double ComputeApproximationCorrelationSquareManually(const HalfStateTypeBySheet& inputStateMask, const HalfStateTypeBySheet& outputStateMask, int startr, int endr, int totalRounds) const {
		int r = endr - startr;

		if (r < 0 || startr < 0 || endr < 0 || startr > totalRounds || endr > totalRounds)
		{
			throw MyException(__func__ + string(": Parameters error"));
		}
		else if (r == 0)
		{
			double approxCor2 = 1.0;
			for (int x = 0; x < 4; x++)
				for (int z = 0; z < 16; z++)
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

			HalfStateTypeBySheet prevOutputStateMaskBase = { 0 };
			vector<pair<int, int>> activeColumns;
			vector<int> activeColumnInputMaskIndexes;
			vector<int> activeColumnOutputMasks;

			int internalDifferenceIndexPrevRoundConstant = xoodooBase.ComputeRoundConstantIndex(totalRounds, endr);

			for (int x = 0; x < 4; x++)
				for (int z = 0; z < 16; z++)
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
				HalfStateTypeBySheet prevOutputStateMask = prevOutputStateMaskBase;

				double cor2CurNonLinearLayer = 1.0;

				for (int activeColumnIndex = 0; activeColumnIndex < activeColumnsNum; activeColumnIndex++)
				{
					const pair<int, int>& activeColumnPos = activeColumns[activeColumnIndex];
					const int& activeColumnInputMaskIndex = activeColumnInputMaskIndexes[activeColumnIndex];
					const int& activeColumnOutputMask = activeColumnOutputMasks[activeColumnIndex];
					const int& activeColumnInputMask = inputMasksPerOutputMaskForChi[activeColumnOutputMask][activeColumnInputMaskIndex];

					xoodooBase.SetColumn(prevOutputStateMask, activeColumnPos.first, activeColumnPos.second, activeColumnInputMask);

					cor2CurNonLinearLayer *= latCor2[activeColumnInputMask][activeColumnOutputMask];
				}

				PropagateThroughRevLinearLayerWithCor2(prevOutputStateMask, cor2CurNonLinearLayer, internalDifferencePerRoundConstants[internalDifferenceIndexPrevRoundConstant]);
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

	double ComputeLinearTrailCorrelation2(const vector<HalfStateTypeBySheet>& inputStateMasks, const vector<HalfStateTypeBySheet>& outputStateMasks, int startr, int endr, int totalRounds) const
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
				HalfStateTypeBySheet curInputStateMask = inputStateMasks[i];
				HalfStateTypeBySheet curOutputStateMask = outputStateMasks[i];

				double cor2CurNonLinearLayer = ComputeApproximationCorrelationSquareManually(curInputStateMask, curOutputStateMask, startr + i, startr + i, totalRounds);

				if (i > 0)
				{
					PropagateThroughLota(curInputStateMask, cor2CurNonLinearLayer, internalDifferencePerRoundConstants[xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + i)]);
				}


				trailCor2 *= cor2CurNonLinearLayer;
			}

			double cor2LastNonLinearLayer = ComputeApproximationCorrelationSquareManually(inputStateMasks[r], outputStateMasks[r], startr + r, endr, totalRounds);

			if (r > 0)
				PropagateThroughLota(inputStateMasks[r], cor2LastNonLinearLayer, internalDifferencePerRoundConstants[xoodooBase.ComputeRoundConstantIndex(totalRounds, endr)]);

			trailCor2 *= cor2LastNonLinearLayer;



			return trailCor2;
		}
	}

	bool VerifyLinearTrail(const vector<HalfStateTypeBySheet>& inputStateMasks, const vector<HalfStateTypeBySheet>& outputStateMasks, int startr, int endr, int totalRounds, double toVerifyTrailCor2) const
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
				HalfStateTypeBySheet curInputStateMask = inputStateMasks[i];
				HalfStateTypeBySheet curOutputStateMask = inputStateMasks[i];
				HalfStateTypeBySheet nextInputStateMask = inputStateMasks[i + 1];

				int curInternalDifferenceIndexNextRoundConstant = xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + i + 1);

				PropagateThroughRevLinearLayer(nextInputStateMask, internalDifferencePerRoundConstants[curInternalDifferenceIndexNextRoundConstant]);

				if (curOutputStateMask != nextInputStateMask)
					return false;
			}

			double trailCor2 = ComputeLinearTrailCorrelation2(inputStateMasks, outputStateMasks, startr, endr, totalRounds);

			if (abs(trailCor2 - toVerifyTrailCor2) < 1e-2)
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

	ILINStateMILPTypeBySheet MilpAddStateVars(GRBModel& model) const
	{
		ILINStateMILPTypeBySheet stateMaskVars;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					stateMaskVars[x][y][z] = model.addVar(0, 1, 0, GRB_BINARY);
		return stateMaskVars;
	}

	HalfStateTypeBySheet MilpReadStateVars(ILINStateMILPTypeBySheet& stateMaskVars) const
	{
		HalfStateTypeBySheet stateMask;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					xoodooBase.SetBit(stateMask, x, y, z, round(stateMaskVars[x][y][z].get(GRB_DoubleAttr_Xn)));
		return stateMask;
	}

	ILINPlaneMILPType MilpAddPlaneVars(GRBModel& model) const
	{
		ILINPlaneMILPType planeMaskVars;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
				planeMaskVars[x][z] = model.addVar(0, 1, 0, GRB_BINARY);
		return planeMaskVars;
	}

	void MilpPropagateThroughRevRhoEast(ILINStateMILPTypeBySheet& stateMaskVars) const
	{
		ILINStateMILPTypeBySheet tmpStateMaskVars = stateMaskVars;
		for (int x = 0; x < 4; x++)
			stateMaskVars[x][1] = LINRotateLeft(tmpStateMaskVars[x][1], -1);

		for (int x = 0; x < 4; x++)
			stateMaskVars[x][2] = LINRotateLeft(tmpStateMaskVars[(x + 2) % 4][2], -8);
	}

	void MilpPropagateThroughRevTheta(GRBModel& model, ILINStateMILPTypeBySheet& stateMaskVars) const
	{
		ILINPlaneMILPType parityMaskVars = MilpAddPlaneVars(model);
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
				MilpXor3(model, stateMaskVars[x][0][z], stateMaskVars[x][1][z], stateMaskVars[x][2][z], parityMaskVars[x][z]);

		ILINStateMILPTypeBySheet newStateMaskVars = MilpAddStateVars(model);

		for (int x = 0; x < 4; x++)
		{
			ILINLaneMILPType one = LINRotateLeft(parityMaskVars[(x + 1) % 4], -5);
			ILINLaneMILPType two = LINRotateLeft(parityMaskVars[(x + 1) % 4], -14);
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					MilpXor3(model, stateMaskVars[x][y][z], one[z], two[z], newStateMaskVars[x][y][z]);
		}

		stateMaskVars = newStateMaskVars;
	}

	void MilpPropagateThroughRevRhoWest(ILINStateMILPTypeBySheet& stateMaskVars) const
	{
		ILINStateMILPTypeBySheet tmpStateMaskVars = stateMaskVars;
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

	void MilpPropagateThroughRevChi(GRBModel& model, ILINStateMILPTypeBySheet& stateMaskVars) const
	{
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
				MilpPropagateThroughRevChiColumn(model, stateMaskVars[x][0][z], stateMaskVars[x][1][z], stateMaskVars[x][2][z]);
	}

	void MilpPropagateThroughRevLinearLayer(GRBModel& model, ILINStateMILPTypeBySheet& stateMaskVars) const
	{
		MilpPropagateThroughRevRhoWest(stateMaskVars);
		MilpPropagateThroughRevTheta(model, stateMaskVars);
		MilpPropagateThroughRevRhoEast(stateMaskVars);
	}

	void MilpPropagateThroughRevNonLinearLayer(GRBModel& model, ILINStateMILPTypeBySheet& stateMaskVars) const
	{
		MilpPropagateThroughRevChi(model, stateMaskVars);
	}


	double ComputeApproximationCorrelationSquareAutomatically(const HalfStateTypeBySheet& inputStateMask, const HalfStateTypeBySheet& outputStateMask, int startr, int endr, int totalRounds) const
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
				for (int z = 0; z < 16; z++)
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

			ILINStateMILPTypeBySheet stateMaskVars = MilpAddStateVars(model);

			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					for (int z = 0; z < 16; z++)
						model.addConstr(stateMaskVars[x][y][z] == xoodooBase.GetBit(outputStateMask, x, y, z));

			vector<ILINStateMILPTypeBySheet> inputStateMaskVarsPerRoundChi;
			vector<ILINStateMILPTypeBySheet> outputStateMaskVarsPerRoundChi;


			int midR = (startr + endr + 1) / 2;
			GRBLinExpr objExp = 0;

			for (int i = endr; i >= startr + 1; i--)
			{
				outputStateMaskVarsPerRoundChi.insert(outputStateMaskVarsPerRoundChi.begin(), stateMaskVars);
				MilpPropagateThroughRevChi(model, stateMaskVars);
				inputStateMaskVarsPerRoundChi.insert(inputStateMaskVarsPerRoundChi.begin(), stateMaskVars);

				if (i == midR)
				{
					for (int x = 0; x < 4; x++)
						for (int y = 0; y < 3; y++)
							for (int z = 0; z < 16; z++)
								objExp += stateMaskVars[x][y][z];
				}

				MilpPropagateThroughRevLinearLayer(model, stateMaskVars);
			}

			outputStateMaskVarsPerRoundChi.insert(outputStateMaskVarsPerRoundChi.begin(), stateMaskVars);
			MilpPropagateThroughRevChi(model, stateMaskVars);
			inputStateMaskVarsPerRoundChi.insert(inputStateMaskVarsPerRoundChi.begin(), stateMaskVars);

			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					for (int z = 0; z < 16; z++)
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
					vector<HalfStateTypeBySheet> inputStateMasks, outputStateMasks;
					for (int i = 0; i <= r; i++)
					{
						HalfStateTypeBySheet inputStateMask = MilpReadStateVars(inputStateMaskVarsPerRoundChi[i]);
						HalfStateTypeBySheet outputStateMask = MilpReadStateVars(outputStateMaskVarsPerRoundChi[i]);
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

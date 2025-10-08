#pragma once
#include<iostream>
#include<vector>
#include<iomanip>
#include<map>
#include<random>
#include<omp.h>
#include"xoodoo.h"
#include"linear.h"
#include"callback.h"


using namespace std;

using RDLLaneType = array<int, 32>;
using RDLSheetType = array<RDLLaneType, 3>;
using RDLPlaneType = array<RDLLaneType, 4>;
using RDLStateTypeBySheet = array<RDLSheetType, 4>;

RDLStateTypeBySheet RDLStateTypeBySheetFromStateTypeBySheet(const StateTypeBySheet& state)
{
	RDLStateTypeBySheet rdlState = { 0 };
	for (int x = 0; x < 4; x++)
	{
		for (int y = 0; y < 3; y++)
		{
			for (int z = 0; z < 32; z++)
			{
				if ((state[x][y] >> z) & 0x1)
					rdlState[x][y][z] = 1;
				else
					rdlState[x][y][z] = 0;
			}
		}
	}
	return rdlState;
}

RDLLaneType RotateLeft(const RDLLaneType lane, int shift)
{
	RDLLaneType newLane;
	for (int i = 0; i < 32; i++)
		newLane[i] = lane[(i - shift + 32) % 32];
	return newLane;
}

RDLLaneType RotateRight(const RDLLaneType lane, int shift)
{
	return RotateLeft(lane, 32 - shift);
}

class RotationalDifferentialLinearBase
{
private:
	const XoodooBase& xoodooBase;
	const LinearBase& linearBase;
	const RotationalDifferentialBase& rdBase;
	int milpThreads;

	array<array<int, 8>, 8> dlct;
	array<array<double, 8>, 8> dlctCor;
	array<LaneType, 12> rotationalDifferencePerRoundConstants;

	// first index is the output linear mask, second index is the mask indicating the input unknown difference bits, third index is the input fixed difference concatenating the input linear mask of the nonfixed bits
	array<array<array<double,8>,8>,8> RDLCor2PerOutputLinearMaskForChi;

	// first index is the output linear mask, second index is the mask indicating the input unknown difference bits, third index is the input fixed difference, and the vector contains all the possible input linear masks, which are restricted in the positions of unknown difference bits
	array<array<array<vector<int>, 8>, 8>, 8> RDLInputMasksPerOutputMaskForChi;

	// first index is the input difference, second index is the mask indicating the output unknown linear mask, third index is the output fixed mask concatenating the output difference of the nonfixed bits
	array<array<array<double, 8>, 8>, 8> RDLProbPerInputDifferenceForChi;

	array<array<array<vector<int>, 8>, 8>, 8> RDLOutputDiffsPerInputDifferenceForChi;
	

	// the inequalties for propagation through dlct
	vector<array<int, 7>> milpChiInequalities;

	void GenerateDLCTForChi()
	{
		for (int a = 0; a < 8; a++)
			for (int b = 0; b < 8; b++)
				dlct[a][b] = 0;

		for(int xdiff = 0; xdiff < (1<<3); xdiff ++)
			for (int ymask = 0; ymask < (1 << 3); ymask++)
			{
				for (int x0 = 0; x0 < (1 << 3); x0++)
				{
					int x1 = x0 ^ xdiff;
					int y0 = xoodooBase.ChiColumnByOperation(x0);
					int y1 = xoodooBase.ChiColumnByOperation(x1);
					int ydiff = y0 ^ y1;
					if(HammingWeight(ydiff & ymask) % 2 == 0)
						dlct[xdiff][ymask]++;
				}
			}

		for (int xdiff = 0; xdiff < (1 << 3); xdiff++)
			for (int ymask = 0; ymask < (1 << 3); ymask++)
			{
				dlct[xdiff][ymask] = 2*dlct[xdiff][ymask] - (1 << 3);
				dlctCor[xdiff][ymask] = dlct[xdiff][ymask] / 8.0;
			}
	}


	void GenerateRDLProbPerInputDifferenceForChi()
	{
		vector<vector<int>> outputDiffPerInputDiffForChi = rdBase.getOutputDiffsPerInputDiffForChi();
		array<array<double, 8>, 8> ddtProb = rdBase.getDDTProb();

		for (int inputDiff = 0; inputDiff < (1 << 3); inputDiff++)
		{
			vector<int> possibleOutputDiffs = outputDiffPerInputDiffForChi[inputDiff];
			for (int outputUnknownLinearMaskBitsMask = 0; outputUnknownLinearMaskBitsMask < (1 << 3); outputUnknownLinearMaskBitsMask++)
			{
				RDLProbPerInputDifferenceForChi[inputDiff][outputUnknownLinearMaskBitsMask].fill(0.0);
				int outputFixedLinearMaskBitsMask = (1 << 3) - 1 - outputUnknownLinearMaskBitsMask;
				vector<int> outputFixedLinearMaskBitsPos;
				for (int i = 0; i < 3; i++)
				{
					if ((outputFixedLinearMaskBitsMask >> i) & 1)
					{
						outputFixedLinearMaskBitsPos.push_back(i);
					}
				}

				int outputFixedLinearMaskBitsNum = outputFixedLinearMaskBitsPos.size();

				for (auto& outputDiff : possibleOutputDiffs)
				{
					double prob = ddtProb[inputDiff][outputDiff];
					int outputDiffForFixedLinearMaskBits = outputDiff & outputFixedLinearMaskBitsMask;
					int outputDiffForUnknownLinearMaskBits = outputDiff & outputUnknownLinearMaskBitsMask;

					for (int v = 0; v < (1 << outputFixedLinearMaskBitsNum); v++)
					{
						int outputFixedLinearMaskBits = 0;
						for (int i = 0; i < outputFixedLinearMaskBitsNum; i++)
						{
							if ((v >> i) & 1)
							{
								outputFixedLinearMaskBits |= (1 << outputFixedLinearMaskBitsPos[i]);
							}
						}

						bool isSignFlipped = HammingWeight(outputFixedLinearMaskBits & outputDiffForFixedLinearMaskBits) & 1;
						if (isSignFlipped)
							RDLProbPerInputDifferenceForChi[inputDiff][outputUnknownLinearMaskBitsMask][outputFixedLinearMaskBits | outputDiffForUnknownLinearMaskBits] += -prob;
						else
							RDLProbPerInputDifferenceForChi[inputDiff][outputUnknownLinearMaskBitsMask][outputFixedLinearMaskBits | outputDiffForUnknownLinearMaskBits] += prob;
					}
				}
			}
		}
	}

	void GenerateRDLOutputDiffsPerInputDifferenceForChi()
	{
		for (int inputDiff = 0; inputDiff < (1 << 3); inputDiff++)
		{
			for (int outputUnknownLinearMaskBitsMask = 0; outputUnknownLinearMaskBitsMask < (1 << 3); outputUnknownLinearMaskBitsMask++)
			{
				int outputFixedLinearMaskBitsMask = (1 << 3) - 1 - outputUnknownLinearMaskBitsMask;

				array<double, 8> RDLProbs = RDLProbPerInputDifferenceForChi[inputDiff][outputUnknownLinearMaskBitsMask];
				for (int v = 0; v < (1 << 3); v++)
				{
					if (RDLProbs[v] != 0.0)
					{
						int outputFixedLinearMaskBits = v & outputFixedLinearMaskBitsMask;
						int outputDiffForUnknownLinearMaskBits = v & outputUnknownLinearMaskBitsMask;
						RDLOutputDiffsPerInputDifferenceForChi[inputDiff][outputUnknownLinearMaskBitsMask][outputFixedLinearMaskBits].push_back(outputDiffForUnknownLinearMaskBits);
					}
				}
			}
		}
	}

	void GenerateRDLCor2PerOutputLinearMaskForChi()
	{
		vector<vector<int>> inputLinearMasksPerOutputLinearMaskForChi = linearBase.getInputMasksPerOutputMaskForChi();
		array<array<double, 8>, 8> latCor2 = linearBase.getLATCor2();

		for (int outputLinearMask = 0; outputLinearMask < (1 << 3); outputLinearMask++)
		{
			// cout << "Current output linear mask " << outputLinearMask << endl;
			vector<int> possibleInputLinearMasks = inputLinearMasksPerOutputLinearMaskForChi[outputLinearMask];
			for (int inputUnknownDiffBitsMask = 0; inputUnknownDiffBitsMask < (1 << 3); inputUnknownDiffBitsMask++)
			{
				// cout << "Current unknown difference bits mask " << inputUnknownDiffBitsMask << endl;
				RDLCor2PerOutputLinearMaskForChi[outputLinearMask][inputUnknownDiffBitsMask].fill(0.0);
				int inputFixedDiffBitsMask = (1 << 3) - 1 - inputUnknownDiffBitsMask;
				vector<int> inputFixedDiffBitsPos;
				for (int i = 0; i < 3; i++)
				{
					if ((inputFixedDiffBitsMask >> i) & 1)
					{
						inputFixedDiffBitsPos.push_back(i);
					}
				}
				int inputFixedDiffBitsNum = inputFixedDiffBitsPos.size();

				for (auto& inputLinearMask : possibleInputLinearMasks)
				{
					double cor2 = latCor2[inputLinearMask][outputLinearMask];
					int inputLinearMaskForFixedDiffBits = inputLinearMask & inputFixedDiffBitsMask;
					int inputLinearMaskForUnknownDiffbits = inputLinearMask & inputUnknownDiffBitsMask;

					for (int v = 0; v < (1 << inputFixedDiffBitsNum); v++)
					{

						int inputFixedDiff = 0;
						for (int i = 0; i < inputFixedDiffBitsNum; i++)
						{
							if ((v >> i) & 1)
							{
								inputFixedDiff |= (1 << inputFixedDiffBitsPos[i]);
							}
						}

						bool isSignFlipped = HammingWeight(inputFixedDiff & inputLinearMaskForFixedDiffBits) & 1;
						// cout << "Input fixed difference " << inputFixedDiff << endl;
						// cout << "Input Linear mask for fixed difference " << inputLinearMaskForFixedDiffBits << endl;
						if (isSignFlipped)
							RDLCor2PerOutputLinearMaskForChi[outputLinearMask][inputUnknownDiffBitsMask][inputFixedDiff | inputLinearMaskForUnknownDiffbits] += -cor2;
						else
							RDLCor2PerOutputLinearMaskForChi[outputLinearMask][inputUnknownDiffBitsMask][inputFixedDiff | inputLinearMaskForUnknownDiffbits] += cor2;
					}
				}
				
			}
		}
	}

	void GenerateRDLInputMasksPerOutputMaskForChi()
	{
		for (int outputLinearMask = 0; outputLinearMask < (1 << 3); outputLinearMask++)
		{
			for (int inputUnknownDiffBitsMask = 0; inputUnknownDiffBitsMask < (1 << 3); inputUnknownDiffBitsMask++)
			{
                int inputFixedDiffBitsMask = (1 << 3) - 1 - inputUnknownDiffBitsMask;

				array<double, 8> RDLCor2s = RDLCor2PerOutputLinearMaskForChi[outputLinearMask][inputUnknownDiffBitsMask];
				for (int v = 0; v < (1 << 3); v++)
				{
					if (RDLCor2s[v] != 0.0)
					{
						int inputFixedDiffBits = v & inputFixedDiffBitsMask;
						int inputLinearMaskForUnknownDiffBits = v & inputUnknownDiffBitsMask;
						RDLInputMasksPerOutputMaskForChi[outputLinearMask][inputUnknownDiffBitsMask][inputFixedDiffBits].push_back(inputLinearMaskForUnknownDiffBits);
					}
				}
			}
		}
	}

	int xor2(int a, int b) const
	{
		if (a == -1 || b == -1)
			return -1;

		return a ^ b;
	}

	int xor3(int a, int b, int c) const
	{
		if (a == -1 || b == -1 || c == -1)
			return -1;

		return a ^ b ^ c;
	}

	bool SwitchToNextPossibleOutputDifferencesOrInputLinearMasks(vector<int>& activeColumnsOutputDiffOrInputLinearMaskIndexes, const vector<vector<int>>& activeColumnsPossibleOutputDiffsOrInputLinearMasks) const
	{
		int activeColumnsNum = activeColumnsOutputDiffOrInputLinearMaskIndexes.size();
		for (int i = 0; i < activeColumnsNum; i++)
		{
			activeColumnsOutputDiffOrInputLinearMaskIndexes[i]++;
			if (activeColumnsOutputDiffOrInputLinearMaskIndexes[i] < activeColumnsPossibleOutputDiffsOrInputLinearMasks[i].size())
				return true;
			else
				activeColumnsOutputDiffOrInputLinearMaskIndexes[i] = 0;
		}
		return false;
	}

public:
	RotationalDifferentialLinearBase(const XoodooBase & xoodooBase, const LinearBase& linearBase, const RotationalDifferentialBase& rdBase, int milpThreads = 2) : xoodooBase(xoodooBase), linearBase(linearBase), rdBase(rdBase), milpThreads(milpThreads)
	{
		GenerateRDLCor2PerOutputLinearMaskForChi();
		GenerateRDLProbPerInputDifferenceForChi();
		GenerateRDLInputMasksPerOutputMaskForChi();
		GenerateRDLOutputDiffsPerInputDifferenceForChi();
		GenerateDLCTForChi();

		rotationalDifferencePerRoundConstants = rdBase.getRotationalDifferencePerRoundConstants();

		milpChiInequalities = { {1, -1, 0, 0, 1, -1, 0}, {1, 0, -1, 0, -1, 1, 0}, {1, 0, 0, -1, -1, 0, 1}, {1, -1, 0, 1, 0, 0, -1}, {1, 0, 1, -1, 0, -1, 0}, {1, 0, -1, 0, 0, 1, -1}, {1, 1, -1, 0, -1, 0, 0}, {1, -1, 0, 0, 1, 0, -1}, {1, -1, 1, 0, 0, -1, 0}, {1, 0, -1, 1, 0, 0, -1}, {1, 0, 0, -1, 0, -1, 1}, {1, 1, 0, -1, -1, 0, 0} };

	}

	array<array<array<double, 8>, 8>, 8> getRDLCor2PerOutputLinearMaskForChi() const
	{
		return RDLCor2PerOutputLinearMaskForChi;
	}

	array<array<array<vector<int>, 8>, 8>, 8> getRDLInputMasksPerOutputMaskForChi() const
	{
		return RDLInputMasksPerOutputMaskForChi;
	}

	array<array<array<double, 8>, 8>, 8> getRDLProbPerInputDifferenceForChi() const
	{
		return RDLProbPerInputDifferenceForChi;
	}

	array<array<array<vector<int>, 8>, 8>, 8> getRDLOutputDiffsPerInputDifferenceForChi() const
	{
		return RDLOutputDiffsPerInputDifferenceForChi;
	}

	array<array<double, 8>, 8> getDLCTCor() const
	{
		return dlctCor;
	}

	void OutputRDLState(ostream& os, const RDLStateTypeBySheet& state) const
	{
		for (int y = 0; y < 3; y++)
		{
			for (int x = 0; x < 4; x++)
			{
				for (int z = 0; z < 32; z++)
				{
					os << setw(2) << state[x][y][z];
				}
				os << "|";
			}
			os << endl;
		}

		os << endl;
	}

// =====================================================================================
// ================= Functions below are for differential and linear propagation in the middle part ==================
// =====================================================================================

	void RDLDiffPropagateThroughRhoWest(RDLStateTypeBySheet& state) const
	{
		RDLStateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = tmp[(x + 3) % 4][1];

		for (int x = 0; x < 4; x++)
			state[x][2] = RotateLeft(tmp[x][2], 11);
	}

	void RDLLinPropagateThroughRevRhoWest(RDLStateTypeBySheet& state) const
	{
		RDLStateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = tmp[(x + 1) % 4][1];

		for (int x = 0; x < 4; x++)
			state[x][2] = RotateRight(tmp[x][2], 11);
	}

	void RDLDiffPropagateThroughRhoEast(RDLStateTypeBySheet& state) const
	{
		RDLStateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = RotateLeft(tmp[x][1], 1);

		for (int x = 0; x < 4; x++)
			state[x][2] = RotateLeft(tmp[(x + 2) % 4][2], 8);
	}

	void RDLLinPropagateThroughRevRhoEast(RDLStateTypeBySheet& state) const
	{
		RDLStateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = RotateRight(tmp[x][1], 1);

		for (int x = 0; x < 4; x++)
			state[x][2] = RotateRight(tmp[(x + 2) % 4][2], 8);
	}

	void RDLDiffPropagateThroughLota(RDLStateTypeBySheet& state, const LaneType rcDiff) const
	{
		for (int z = 0; z < 32; z++)
			if (state[0][0][z] != -1)
				state[0][0][z] = xor2(state[0][0][z], (rcDiff >> z) & 1);
	}

	void RDLDiffPropagateThroughTheta(RDLStateTypeBySheet& state) const
	{
		RDLPlaneType parity;
		for(int x = 0; x < 4; x++)
			for(int z = 0; z <32 ;z++)
				parity[x][z] = xor3(state[x][0][z], state[x][1][z], state[x][2][z]);

		for (int x = 0; x < 4; x++) 
		{
			RDLLaneType one = RotateLeft(parity[(x + 3) % 4], 5);
			RDLLaneType two = RotateLeft(parity[(x + 3) % 4], 14);
			for(int y = 0; y < 3; y++)
				for(int z = 0; z <32; z++)
					state[x][y][z] = xor3(state[x][y][z], one[z], two[z]);

		}
	}

	void RDLLinPropagateThroughRevTheta(RDLStateTypeBySheet& state) const
	{
		RDLPlaneType parity;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 32; z++)
				parity[x][z] = xor3(state[x][0][z], state[x][1][z], state[x][2][z]);

		for (int x = 0; x < 4; x++)
		{
			RDLLaneType one = RotateRight(parity[(x + 1) % 4], 5);
			RDLLaneType two = RotateRight(parity[(x + 1) % 4], 14);
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					state[x][y][z] = xor3(state[x][y][z], one[z], two[z]);

		}
	}

	void RDLDiffPropagateThroughChiColumn(int& y0, int& y1, int& y2) const
	{
		if (y0 == 0 && y1 == 0 && y2 == 0)
			;
		else if (y0 == 1 && y1 == 0 && y2 == 0)
		{
			y1 = -1;
			y2 = -1;
		}
		else if (y0 == 0 && y1 == 1 && y2 == 0)
		{
			y0 = -1;
			y2 = -1;
		}
		else if (y0 == 0 && y1 == 0 && y2 == 1)
		{
			y0 = -1;
			y1 = -1;
		}
		else
		{
			y0 = -1;
			y1 = -1;
			y2 = -1;
		}
	}

	void RDLLinPropagateThroughRevChiColumn(int& y0, int& y1, int& y2) const
	{
		RDLDiffPropagateThroughChiColumn(y0, y1, y2);
	}

	void RDLDiffPropagateThroughChi(RDLStateTypeBySheet& state) const
	{
		for(int x = 0; x < 4; x++)
			for(int z = 0; z < 32; z++)
				RDLDiffPropagateThroughChiColumn(state[x][0][z], state[x][1][z], state[x][2][z]);
	}

	void RDLLinPropagateThroughRevChi(RDLStateTypeBySheet& state) const
	{
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 32; z++)
				RDLLinPropagateThroughRevChiColumn(state[x][0][z], state[x][1][z], state[x][2][z]);
	}

	void RDLDiffPropagateThroughLinearLayer(RDLStateTypeBySheet& state, const LaneType rcDiff) const
	{
		RDLDiffPropagateThroughRhoEast(state);
		RDLDiffPropagateThroughTheta(state);
		RDLDiffPropagateThroughRhoWest(state);
		RDLDiffPropagateThroughLota(state, rcDiff);
	}

	void RDLLinPropagateThroughRevLinearLayer(RDLStateTypeBySheet& state) const
	{
		RDLLinPropagateThroughRevRhoWest(state);
		RDLLinPropagateThroughRevTheta(state);
		RDLLinPropagateThroughRevRhoEast(state);
	}

	void RDLDiffPropagateThroughNonLinearLayer(RDLStateTypeBySheet& state) const
	{
		RDLDiffPropagateThroughChi(state);
	}

	void RDLLinPropagateThroughRevNonLinearLayer(RDLStateTypeBySheet& state) const
	{
		RDLLinPropagateThroughRevChi(state);
	}

// =====================================================================================
// ================= Functions below are for computing the correlation of RDL distinguisher in the middle part ==================
// =====================================================================================


	int RDLFowardSingleStateDiffThrouChiEstimation(const StateTypeBySheet& inputStateDiff, const RDLStateTypeBySheet& outputRDLLinMask) const
	{
		array<array<int, 32>, 4> outputUnknownLinearMaskBitsMask = { 0 };
		array<array<int, 32>, 4> outputFixedLinearMaskBits = { 0 };
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 32; z++)
			{
				for (int y = 0; y < 3; y++)
					if (outputRDLLinMask[x][y][z] == -1)
						outputUnknownLinearMaskBitsMask[x][z] |= (1 << y);
					else if (outputRDLLinMask[x][y][z] == 1)
						outputFixedLinearMaskBits[x][z] |= (1 << y);
			}



		double estimatedNumLog2 = 0;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 32; z++)
			{
				int inputDiff = xoodooBase.GetColumn(inputStateDiff, x, z);

				const vector<int>& possibleOutputDiffs = RDLOutputDiffsPerInputDifferenceForChi[inputDiff][outputUnknownLinearMaskBitsMask[x][z]][outputFixedLinearMaskBits[x][z]];

				if (possibleOutputDiffs.size() == 0)
				{
					// cout << "x = " << x << " z = " << z << " has no possible output differences" << endl;
					// xoodooBase.OutputState(cout, inputStateDiff);
					// OutputRDLState(cout, outputRDLLinMask);
					// cout << "inputDiff = " << inputDiff << endl;
					// cout << "outputUnknownLinearMaskBitsMask = " << outputUnknownLinearMaskBitsMask[x][z] << endl;
					return -1;
				}

				if (inputDiff != 0)
				{
					estimatedNumLog2 += log2(possibleOutputDiffs.size());
				}


			}

		return estimatedNumLog2;
	}

	int RDLBackwardSingleStateLinearMaskThroughRevChiEstimation(const RDLStateTypeBySheet& inputStateRDLDiff, const StateTypeBySheet& outputLinMask) const
	{
		array<array<int, 32>, 4> inputUnknownDiffBitsMask = { 0 };
		array<array<int, 32>, 4> inputFixedDiffBits = { 0 };
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 32; z++)
				for (int y = 0; y < 3; y++)
					if (inputStateRDLDiff[x][y][z] == -1)
						inputUnknownDiffBitsMask[x][z] |= (1 << y);
					else if (inputStateRDLDiff[x][y][z] == 1)
						inputFixedDiffBits[x][z] |= (1 << y);




		double estimatedNumLog2 = 0;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 32; z++)
			{
				int outputLinearMask = xoodooBase.GetColumn(outputLinMask, x, z);
				const vector<int>& possibleInputLinearMasks = RDLInputMasksPerOutputMaskForChi[outputLinearMask][inputUnknownDiffBitsMask[x][z]][inputFixedDiffBits[x][z]];

				if (possibleInputLinearMasks.size() == 0)
				{
					/*
					cout << "x = " << x << " z = " << z << " has no possible input linear masks" << endl;
					xoodooBase.OutputState(cout, outputLinMask);
					OutputRDLState(cout, inputStateRDLDiff);
					cout << "outputMask = " << outputLinearMask << endl;
					cout << "inputUnknownDiffBitsMask = " << inputUnknownDiffBitsMask[x][z] << endl;
					*/
					return -1;
				}
				
				if (outputLinearMask != 0)
				{
					estimatedNumLog2 += log2(possibleInputLinearMasks.size());

				}
			}

		return estimatedNumLog2;
	}

	

	vector<pair<StateTypeBySheet, double>> RDLForwardSingleStateDiffThroughChi(const StateTypeBySheet& inputStateDiff, const RDLStateTypeBySheet& outputRDLLinMask) const
	{
		vector<pair<StateTypeBySheet, double>> outputStateDiffAndProbs;

		array<array<int, 32>, 4> outputUnknownLinearMaskBitsMask = { 0 };
		array<array<int, 32>, 4> outputFixedLinearMaskBits = { 0 };
		for(int x = 0; x < 4; x++)
			for (int z = 0; z < 32; z++)
			{
				for (int y = 0; y < 3; y++)
					if (outputRDLLinMask[x][y][z] == -1)
						outputUnknownLinearMaskBitsMask[x][z] |= (1 << y);
					else if(outputRDLLinMask[x][y][z] == 1)
						outputFixedLinearMaskBits[x][z] |= (1 << y);
			}


		vector<pair<int, int>> activeColumns;
		vector<int> activeColumnsOutputDiffIndexes;
		vector<int> activeColumnsInputDiffs;
		vector<vector<int>> activeColumnsPossibleOutputDiffs;
		vector<vector<double>> activeColumnsPossibleProbs;

		double estimatedNumLog2 = 0;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 32; z++)
			{
				int inputDiff = xoodooBase.GetColumn(inputStateDiff, x, z);
				
				const vector<int>& possibleOutputDiffs = RDLOutputDiffsPerInputDifferenceForChi[inputDiff][outputUnknownLinearMaskBitsMask[x][z]][outputFixedLinearMaskBits[x][z]];

				if (possibleOutputDiffs.size() == 0)
					return outputStateDiffAndProbs;

				if (inputDiff != 0)
				{
					estimatedNumLog2 += log2(possibleOutputDiffs.size());
					activeColumns.push_back({ x, z });
					activeColumnsOutputDiffIndexes.push_back(0);
					activeColumnsInputDiffs.push_back(inputDiff);
					activeColumnsPossibleOutputDiffs.push_back(possibleOutputDiffs);
					vector<double> possibleProbs;
					for (auto& outputDiff : possibleOutputDiffs)
						possibleProbs.push_back(RDLProbPerInputDifferenceForChi[inputDiff][outputUnknownLinearMaskBitsMask[x][z]][outputFixedLinearMaskBits[x][z] | outputDiff]);
					activeColumnsPossibleProbs.push_back(possibleProbs);

					/*
					cout << "Output fixed linear mask bits: " << outputFixedLinearMaskBits[x][z] << endl;
					cout << "Output unknown linear mask bits mask: " << outputUnknownLinearMaskBitsMask[x][z] << endl;
					cout << "Input difference: " << inputDiff << endl;
					cout << "Possible output difference size for this column: " << possibleOutputDiffs.size() << endl;
					cout << "Possible probs for this column: ";
					for (auto& prob : possibleProbs)
						cout << prob << " ";
					cout << endl;
					*/
				}

				
			}

		// cout << "Estimated results num log2: " << estimatedNumLog2 << endl;

		
		
		
		

		int activeColumnsNum = activeColumns.size();
		// cout << "Active columns num: " << activeColumnsNum << endl;

		bool isAllActiveColumnsProcessed = false;
		while (!isAllActiveColumnsProcessed)
		{
			StateTypeBySheet outputStateDiff = { 0 };
			double prob = 1.0;

			for (int activeColumnIndex = 0; activeColumnIndex < activeColumnsNum; activeColumnIndex++)
			{
				const pair<int, int>& activeColumnPos = activeColumns[activeColumnIndex];
				const int& activeColumnOutputDiffIndex = activeColumnsOutputDiffIndexes[activeColumnIndex];
				const int& activeColumnInputDiff = activeColumnsInputDiffs[activeColumnIndex];
				const vector<int>& activeColumnPossibleOutputDiffs = activeColumnsPossibleOutputDiffs[activeColumnIndex];
				const vector<double>& activeColumnPossibleProbs = activeColumnsPossibleProbs[activeColumnIndex];
				const int& activeColumnCurOutputDiff = activeColumnPossibleOutputDiffs[activeColumnOutputDiffIndex];
				const double& activeColumnCurProb = activeColumnPossibleProbs[activeColumnOutputDiffIndex];

				xoodooBase.SetColumn(outputStateDiff, activeColumnPos.first, activeColumnPos.second, activeColumnCurOutputDiff);

				prob *= activeColumnCurProb;
			}

			outputStateDiffAndProbs.push_back({ outputStateDiff, prob });

			if (SwitchToNextPossibleOutputDifferencesOrInputLinearMasks(activeColumnsOutputDiffIndexes, activeColumnsPossibleOutputDiffs))
				isAllActiveColumnsProcessed = false;
			else
				isAllActiveColumnsProcessed = true;
		}
		
		return outputStateDiffAndProbs;



	}

	vector<pair<StateTypeBySheet, double>> RDLBackwardSingleStateLinearMaskThroughRevChi(const RDLStateTypeBySheet& inputStateRDLDiff, const StateTypeBySheet& outputLinMask) const
	{
		// OutputRDLState(cout, inputStateRDLDiff);
		// xoodooBase.OutputState(cout, outputLinMask);

		vector<pair<StateTypeBySheet, double>> inputStateRDLLinearMaskAndCor2s;

		array<array<int, 32>, 4> inputUnknownDiffBitsMask = { 0 };
		array<array<int, 32>, 4> inputFixedDiffBits = { 0 };
		for(int x = 0; x < 4;x++)
			for(int z = 0; z < 32; z++)
				for(int y = 0; y < 3; y++)
					if(inputStateRDLDiff[x][y][z] == -1)
						inputUnknownDiffBitsMask[x][z] |= (1 << y);
					else if(inputStateRDLDiff[x][y][z] == 1)
						inputFixedDiffBits[x][z] |= (1 << y);
		


		vector<pair<int, int>> activeColumns;
		vector<int>	activeColumnsInputLinearMaskIndexes;
		vector<int> activeColumnsOutputLinearMasks;
		vector<vector<int>> activeColumnsPossibleInputLinearMasks;
		vector<vector<double>> activeColumnsPossibleCor2s;

		double estimatedNumLog2 = 0;
		for(int x = 0; x < 4;x++)
			for (int z = 0; z < 32; z++)
			{
				int outputLinearMask = xoodooBase.GetColumn(outputLinMask, x, z);
				const vector<int>& possibleInputLinearMasks = RDLInputMasksPerOutputMaskForChi[outputLinearMask][inputUnknownDiffBitsMask[x][z]][inputFixedDiffBits[x][z]];

				if (possibleInputLinearMasks.size() == 0)
					return inputStateRDLLinearMaskAndCor2s;

				if (outputLinearMask != 0)
				{
					estimatedNumLog2 += log2(possibleInputLinearMasks.size());

					activeColumns.push_back({ x, z });
					activeColumnsInputLinearMaskIndexes.push_back(0);
					activeColumnsOutputLinearMasks.push_back(outputLinearMask);
					activeColumnsPossibleInputLinearMasks.push_back(possibleInputLinearMasks);
					vector<double> possibleCor2s;
					for (auto& inputLinearMask : possibleInputLinearMasks)
						possibleCor2s.push_back(RDLCor2PerOutputLinearMaskForChi[outputLinearMask][inputUnknownDiffBitsMask[x][z]][inputFixedDiffBits[x][z] | inputLinearMask]);
					activeColumnsPossibleCor2s.push_back(possibleCor2s);

					/*
					cout << "Input fixed difference bits: " << inputFixedDiffBits[x][z] << endl;
					cout << "Input unknown difference bits mask: " << inputUnknownDiffBitsMask[x][z] << endl;
					cout << "Ouput linear mask: " << outputLinearMask << endl;
					cout << "Possible input linear masks size for this column: " << possibleInputLinearMasks.size() << endl;
					cout << "Possible cor2s for this column: ";
					for (auto& cor2 : possibleCor2s)
						cout << cor2 << " ";
					cout << endl;
					*/
				}
			}


		int activeColumnsNum = activeColumns.size();
	

		bool isAllActiveColumnsProcessed = false;
		while (!isAllActiveColumnsProcessed)
		{
			StateTypeBySheet inputStateLinearMask = { 0 };
			double cor2 = 1.0;

			for (int activeColumnIndex = 0; activeColumnIndex < activeColumnsNum; activeColumnIndex++)
			{
				const pair<int, int>& activeColumnPos = activeColumns[activeColumnIndex];
				const int& activeColumnInputLinearMaskIndex = activeColumnsInputLinearMaskIndexes[activeColumnIndex];
				const int& activeColumnOutputLinearMask = activeColumnsOutputLinearMasks[activeColumnIndex];
				const vector<int>& activeColumnPossibleInputLinearMasks = activeColumnsPossibleInputLinearMasks[activeColumnIndex];
				const vector<double>& activeColumnPossibleCor2s = activeColumnsPossibleCor2s[activeColumnIndex];
				const int& activeColumnCurInputLinearMask = activeColumnPossibleInputLinearMasks[activeColumnInputLinearMaskIndex];
				const double& activeColumnCurCor2 = activeColumnPossibleCor2s[activeColumnInputLinearMaskIndex];

				xoodooBase.SetColumn(inputStateLinearMask, activeColumnPos.first, activeColumnPos.second, activeColumnCurInputLinearMask);

				cor2 *= activeColumnCurCor2;
			}

			inputStateRDLLinearMaskAndCor2s.push_back({ inputStateLinearMask, cor2 });

			if (SwitchToNextPossibleOutputDifferencesOrInputLinearMasks(activeColumnsInputLinearMaskIndexes, activeColumnsPossibleInputLinearMasks))
				isAllActiveColumnsProcessed = false;
			else
				isAllActiveColumnsProcessed = true;
		}

		return inputStateRDLLinearMaskAndCor2s;
	}

	
	
	void RDLForwardExpansionByOneRound(vector<pair<StateTypeBySheet, StateTypeBySheet>>& DiffAndLinearMaskPairs, vector<double>& preComputedCors, vector<pair<RDLStateTypeBySheet, RDLStateTypeBySheet>>& RDLLinearMaskAndDiffPairs, int& startr, int& endr, int totalRounds, bool isCallback) const
	{
		map<pair<StateTypeBySheet, StateTypeBySheet>, double> m0;
		map<pair<StateTypeBySheet, StateTypeBySheet>, pair<RDLStateTypeBySheet, RDLStateTypeBySheet>> m1;

		int rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + 1);
		LaneType rcDiff = rotationalDifferencePerRoundConstants[rcDiffIndex];

		int pairsNum = DiffAndLinearMaskPairs.size();
		cout << __func__ << ": to expand " << pairsNum << " pairs" << endl;
		cout << "Start round: " << startr << endl;
		cout << "End round: " << endr << endl;

		int threadsNum = omp_get_max_threads();
		#pragma omp parallel for num_threads(threadsNum)
		for (int i = 0; i < pairsNum; i++)
		{
			// cout << __func__ << ": pair index " << i << endl;
			const pair<StateTypeBySheet, StateTypeBySheet> diffAndLinearMaskPair = DiffAndLinearMaskPairs[i];
			const StateTypeBySheet& diff = diffAndLinearMaskPair.first;
			const StateTypeBySheet& linearMask = diffAndLinearMaskPair.second;
			const double& preComputedCor = preComputedCors[i];
			const RDLStateTypeBySheet& RDLLinearMask = RDLLinearMaskAndDiffPairs[i].first;
			const RDLStateTypeBySheet& RDLDiff = RDLLinearMaskAndDiffPairs[i].second;

			
			RDLStateTypeBySheet newRDLLinearMask = RDLStateTypeBySheetFromStateTypeBySheet(linearMask);
			for (int j = endr; j > startr + 1; j--)
			{
				RDLLinPropagateThroughRevNonLinearLayer(newRDLLinearMask);
				RDLLinPropagateThroughRevLinearLayer(newRDLLinearMask);
			}
			
			int forwardEstimation = RDLFowardSingleStateDiffThrouChiEstimation(diff, RDLLinearMask);
			if (forwardEstimation == -1)
				continue;
			else
			{
				; // cout << "pair index " << i << " estimation " << forwardEstimation << endl;
			}

			vector<pair<StateTypeBySheet, double>> newDiffAndProbs;
			if(!isCallback)
				newDiffAndProbs = RDLForwardSingleStateDiffThroughChi(diff, RDLLinearMask);
			else
				newDiffAndProbs = RDLForwardSingleStateDiffThroughChiAutomatically(diff, linearMask, startr, endr, totalRounds);

			#pragma omp critical 
			{
				for (auto& diffAndProb : newDiffAndProbs)
				{
					StateTypeBySheet stateDiff = diffAndProb.first;
					double prob = diffAndProb.second;
					rdBase.PropagateThroughLinearLayer(stateDiff, rcDiff);

					if (m0.find({ stateDiff, linearMask }) != m0.end())
						m0[{stateDiff, linearMask}] += prob * preComputedCor;
					else
					{
						m0[{stateDiff, linearMask}] = prob * preComputedCor;

						RDLStateTypeBySheet newRDLDiff = RDLStateTypeBySheetFromStateTypeBySheet(stateDiff);
						for (int j = startr + 1; j < endr; j++)
						{
							int tmpRcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, j + 1);
							LaneType tmpRcDiff = rotationalDifferencePerRoundConstants[tmpRcDiffIndex];
							RDLDiffPropagateThroughNonLinearLayer(newRDLDiff);
							RDLDiffPropagateThroughLinearLayer(newRDLDiff, tmpRcDiff);
						}

						m1[{stateDiff, linearMask}] = { newRDLLinearMask, newRDLDiff };
					}
				}
			}
		}

		DiffAndLinearMaskPairs.clear();
		preComputedCors.clear();
		RDLLinearMaskAndDiffPairs.clear();

		for (auto& item : m0)
		{
			DiffAndLinearMaskPairs.push_back(item.first);
			preComputedCors.push_back(item.second);
			RDLLinearMaskAndDiffPairs.push_back(m1[item.first]);
		}
		
		startr += 1;
	}

	void RDLBackwardExpansionByOneRound(vector<pair<StateTypeBySheet, StateTypeBySheet>>& DiffAndLinearMaskPairs, vector<double>& preComputedCors, vector<pair<RDLStateTypeBySheet, RDLStateTypeBySheet>>& RDLLinearMaskAndDiffPairs, int& startr, int& endr, int totalRounds, bool isCallback) const
	{
		map<pair<StateTypeBySheet, StateTypeBySheet>, double> m0;
		map<pair<StateTypeBySheet, StateTypeBySheet>, pair<RDLStateTypeBySheet, RDLStateTypeBySheet>> m1;

		int rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, endr);
		LaneType rcDiff = rotationalDifferencePerRoundConstants[rcDiffIndex];

		int pairsNum = DiffAndLinearMaskPairs.size();
		cout << __func__ << ": to expand " << pairsNum << " pairs" << endl;
		cout << "Start round: " << startr << endl;
		cout << "End round: " << endr << endl;

		int threadsNum = omp_get_max_threads();
		#pragma omp parallel for num_threads(threadsNum)
		for (int i = 0; i < pairsNum; i++)
		{
			const pair<StateTypeBySheet, StateTypeBySheet> diffAndLinearMaskPair = DiffAndLinearMaskPairs[i];
			const StateTypeBySheet& diff = diffAndLinearMaskPair.first;
			const StateTypeBySheet& linearMask = diffAndLinearMaskPair.second;
			const double& preComputedCor = preComputedCors[i];
			const RDLStateTypeBySheet& RDLLinearMask = RDLLinearMaskAndDiffPairs[i].first;
			const RDLStateTypeBySheet& RDLDiff = RDLLinearMaskAndDiffPairs[i].second;

			RDLStateTypeBySheet newRDLDiff = RDLStateTypeBySheetFromStateTypeBySheet(diff);
			for (int j = startr; j < endr - 1; j++)
			{
				int tmpRcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, j + 1);
				LaneType tmpRcDiff = rotationalDifferencePerRoundConstants[tmpRcDiffIndex];
				RDLDiffPropagateThroughNonLinearLayer(newRDLDiff);
				RDLDiffPropagateThroughLinearLayer(newRDLDiff, tmpRcDiff);
			}

			int backwardEstimation = RDLBackwardSingleStateLinearMaskThroughRevChiEstimation(RDLDiff, linearMask);
			if (backwardEstimation == -1)
				continue;
			else
			{
				cout << "pair index " << i << " estimation " << backwardEstimation << endl;
			}

			vector<pair<StateTypeBySheet, double>> newLinearMaskAndCor2s;
			if(!isCallback)
				newLinearMaskAndCor2s = RDLBackwardSingleStateLinearMaskThroughRevChi(RDLDiff, linearMask);
			else
				newLinearMaskAndCor2s = RDLBackwardSingleStateLinearMaskThroughRevChiAutomatically(diff, linearMask, startr, endr, totalRounds);

			#pragma omp critical 
			{
				for (auto& linearMaskAndCor2 : newLinearMaskAndCor2s)
				{
					StateTypeBySheet stateLinearMask = linearMaskAndCor2.first;
					double cor2 = linearMaskAndCor2.second;
					linearBase.PropagateThroughRevLinearLayerWithCor2(stateLinearMask, cor2, rcDiff);

					if (m0.find({ diff, stateLinearMask }) != m0.end())
						m0[{diff, stateLinearMask}] += preComputedCor * cor2;
					else
					{
						m0[{diff, stateLinearMask}] = preComputedCor * cor2;

						RDLStateTypeBySheet newRDLLinearMask = RDLStateTypeBySheetFromStateTypeBySheet(stateLinearMask);
						for (int j = endr - 1; j > startr; j--)
						{
							RDLLinPropagateThroughRevNonLinearLayer(newRDLLinearMask);
							RDLLinPropagateThroughRevLinearLayer(newRDLLinearMask);
						}

						m1[{diff, stateLinearMask}] = { newRDLLinearMask, newRDLDiff };
					}
				}
			}
		}

		DiffAndLinearMaskPairs.clear();
		preComputedCors.clear();
		RDLLinearMaskAndDiffPairs.clear();

		for (auto& item : m0)
		{
			DiffAndLinearMaskPairs.push_back(item.first);
			preComputedCors.push_back(item.second);
			RDLLinearMaskAndDiffPairs.push_back(m1[item.first]);
		}

		endr -= 1;
	}

	double ExperimentalComputation(const StateTypeBySheet& stateDiff, const StateTypeBySheet& stateLinearMask, int startr, int endr, int totalRounds) const
	{
		random_device rd;
		mt19937 gen(rd());
		uniform_int_distribution<uint32_t> distribution(0, numeric_limits<uint32_t>::max());

		long long totalTestNum = 1ll << 32;
		long long zeroCount = 0;
		int rot = rdBase.getRotationNumber();
		for (long long testIndex = 0; testIndex < totalTestNum; testIndex++)
		{
			StateTypeBySheet state0;
			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					state0[x][y] = distribution(gen);

			StateTypeBySheet state1;
			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
				{
					state1[x][y] = state0[x][y] ^ stateDiff[x][y];
					state1[x][y] = RotateRight(state1[x][y], rot);
				}

			xoodooBase.ReducedRoundForVerification(state0, startr, endr, totalRounds);
			xoodooBase.ReducedRoundForVerification(state1, startr, endr, totalRounds);

			StateTypeBySheet outputStateDiff;
			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					outputStateDiff[x][y] = state0[x][y] ^ RotateLeft(state1[x][y], rot);

			int linearRes = 0;
			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					if ((HammingWeight(outputStateDiff[x][y] & stateLinearMask[x][y]) & 1) == 1)
						linearRes ^= 1;

			if (linearRes == 0)
				zeroCount++;
		}

		

		double zeroProb = zeroCount / (double)totalTestNum;
		// cout << zeroCount << endl;
		// cout << totalTestNum << endl;
		// cout << zeroProb << endl;

		return 2 * zeroProb - 1;
	}
	double ComputeThreeRoundRDLCor(const StateTypeBySheet& inputStateDiff, const StateTypeBySheet& outputStateLinearMask, int startr, int endr, int totalRounds) const
	{
		RDLStateTypeBySheet RDLStateDiff = RDLStateTypeBySheetFromStateTypeBySheet(inputStateDiff);
		for (int i = startr; i < endr; i++)
		{
			int rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(6, i + 1);
			LaneType rcDiff = rotationalDifferencePerRoundConstants[rcDiffIndex];
			RDLDiffPropagateThroughNonLinearLayer(RDLStateDiff);
			RDLDiffPropagateThroughLinearLayer(RDLStateDiff, rcDiff);
		}

		RDLStateTypeBySheet RDLStateLin = RDLStateTypeBySheetFromStateTypeBySheet(outputStateLinearMask);
		for (int i = endr; i > startr; i--)
		{
			RDLLinPropagateThroughRevNonLinearLayer(RDLStateLin);
			RDLLinPropagateThroughRevLinearLayer(RDLStateLin);
		}

		vector<pair<StateTypeBySheet, StateTypeBySheet>> DiffAndLinearMaskPairs = { {inputStateDiff, outputStateLinearMask} };
		vector<pair<RDLStateTypeBySheet, RDLStateTypeBySheet>> RDLLinearMaskAndDiffPairs = { {RDLStateLin, RDLStateDiff} };
		vector<double> preComputedCors = { 1.0 };

		RDLBackwardExpansionByOneRound(DiffAndLinearMaskPairs, preComputedCors, RDLLinearMaskAndDiffPairs, startr, endr, totalRounds, true);
		RDLForwardExpansionByOneRound(DiffAndLinearMaskPairs, preComputedCors, RDLLinearMaskAndDiffPairs, startr, endr, totalRounds, true);
		

		int pairsNum = DiffAndLinearMaskPairs.size();

		cout << __func__ << ": pairs num " << pairsNum << endl;

		double finalCor = 0.0;

		int threadsNum = omp_get_max_threads();

#pragma omp parallel for num_threads(threadsNum)
		for (int i = 0; i < pairsNum; i++)
		{
			cout << __func__ << ": pair index " << i << endl;
			const pair<StateTypeBySheet, StateTypeBySheet>& diffAndLinearMaskPair = DiffAndLinearMaskPairs[i];
			const StateTypeBySheet& diff = diffAndLinearMaskPair.first;
			const StateTypeBySheet& linearMask = diffAndLinearMaskPair.second;
			const RDLStateTypeBySheet& RDLLinearMask = RDLLinearMaskAndDiffPairs[i].first;
			const RDLStateTypeBySheet& RDLDiff = RDLLinearMaskAndDiffPairs[i].second;
			const double& preComputedCor = preComputedCors[i];

			if (preComputedCor == 0.0)
				continue;

			double dlctThroughNonLinearCor = 1.0;

			bool isValid = true;
			for (int x = 0; x < 4; x++)
			{
				for (int z = 0; z < 32; z++)
				{
					int stateDiffColumn = xoodooBase.GetColumn(diff, x, z);
					int stateLinearMaskColumn = xoodooBase.GetColumn(linearMask, x, z);
					dlctThroughNonLinearCor *= dlctCor[stateDiffColumn][stateLinearMaskColumn];
					if (dlctThroughNonLinearCor == 0.0)
					{
						isValid = false;
						break;
					}
				}

				if (!isValid)
					break;
			}

			if (!isValid)
				continue;

			finalCor += preComputedCor * dlctThroughNonLinearCor;
		}

		return finalCor;
	}

	double ComputeFourRoundRDLCor(const StateTypeBySheet& inputStateDiff, const StateTypeBySheet& outputStateLinearMask, int startr, int endr, int totalRounds) const
	{
		RDLStateTypeBySheet RDLStateDiff = RDLStateTypeBySheetFromStateTypeBySheet(inputStateDiff);
		for (int i = startr; i < endr; i++)
		{
			int rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(6, i + 1);
			LaneType rcDiff = rotationalDifferencePerRoundConstants[rcDiffIndex];
			RDLDiffPropagateThroughNonLinearLayer(RDLStateDiff);
			RDLDiffPropagateThroughLinearLayer(RDLStateDiff, rcDiff);
		}

		RDLStateTypeBySheet RDLStateLin = RDLStateTypeBySheetFromStateTypeBySheet(outputStateLinearMask);
		for (int i = endr; i > startr; i--)
		{
			RDLLinPropagateThroughRevNonLinearLayer(RDLStateLin);
			RDLLinPropagateThroughRevLinearLayer(RDLStateLin);
		}

		int initialForwardEstimation = RDLFowardSingleStateDiffThrouChiEstimation(inputStateDiff, RDLStateLin);
		int initialBackwardEstimation = RDLBackwardSingleStateLinearMaskThroughRevChiEstimation(RDLStateDiff, outputStateLinearMask);

		// return ExperimentalComputation(inputStateDiff, outputStateLinearMask, startr, endr, totalRounds);

		/*
		if (initialForwardEstimation > 20 || initialBackwardEstimation > 20)
		{
			cout << "Initial forward estimation: " << initialForwardEstimation << endl;
			cout << "Initial backward estimation: " << initialBackwardEstimation << endl;
			cout << "using experimental computation" << endl;
			return ExperimentalComputation(inputStateDiff, outputStateLinearMask, startr, endr, totalRounds);
		}
		*/



		vector<pair<StateTypeBySheet, StateTypeBySheet>> DiffAndLinearMaskPairs = { {inputStateDiff, outputStateLinearMask} };
		vector<pair<RDLStateTypeBySheet, RDLStateTypeBySheet>> RDLLinearMaskAndDiffPairs = { {RDLStateLin, RDLStateDiff} };
		vector<double> preComputedCors = { 1.0 };

		
		RDLForwardExpansionByOneRound(DiffAndLinearMaskPairs, preComputedCors, RDLLinearMaskAndDiffPairs, startr, endr, totalRounds, true);
		RDLBackwardExpansionByOneRound(DiffAndLinearMaskPairs, preComputedCors, RDLLinearMaskAndDiffPairs, startr, endr, totalRounds, true);
		
		
		
		
		
		
		
		
		
		int pairsNum = DiffAndLinearMaskPairs.size();

		cout << __func__ << ": pairs num " << pairsNum << endl;

		double finalCor = 0.0;

		int threadsNum = omp_get_max_threads();

#pragma omp parallel for num_threads(threadsNum)
		for (int i = 0; i < pairsNum; i++)
		{
			cout << __func__ << ": pair index " << i << endl;
			const pair<StateTypeBySheet, StateTypeBySheet>& diffAndLinearMaskPair = DiffAndLinearMaskPairs[i];
			const StateTypeBySheet& diff = diffAndLinearMaskPair.first;
			const StateTypeBySheet& linearMask = diffAndLinearMaskPair.second;
			const RDLStateTypeBySheet& RDLLinearMask = RDLLinearMaskAndDiffPairs[i].first;
			const RDLStateTypeBySheet& RDLDiff = RDLLinearMaskAndDiffPairs[i].second;
			const double& preComputedCor = preComputedCors[i];

			if (preComputedCor == 0.0)
				continue;

			int forwardEstimation = RDLFowardSingleStateDiffThrouChiEstimation(diff, RDLLinearMask);
			int backwardEstimation = RDLBackwardSingleStateLinearMaskThroughRevChiEstimation(RDLDiff, linearMask);
			

			if (forwardEstimation == -1 || backwardEstimation == -1)
			{
				// cout << __func__ << ": cor is determined as 0" << endl;
				continue;
			}

			cout << __func__ << ": pair index " << i << endl;
			cout << __func__ << ": forward estimation " << forwardEstimation << endl;
			cout << __func__ << ": backward estimation " << backwardEstimation << endl;
			int diffActiveColumnsNum = rdBase.getActiveColumnsNum(diff);
			int linActiveColumnsNum = linearBase.getActiveColumnsNum(linearMask);
			cout << "Difference active columns num " << diffActiveColumnsNum << endl;
			cout << "Linear mask active columns num " << linActiveColumnsNum << endl;

			
			#pragma omp critical 
			{
				
				// cout << "Use milp models to compute." << endl;
				

				
				if(diffActiveColumnsNum > linActiveColumnsNum)
					finalCor += preComputedCor*ComputeTwoRoundRDLCorBackwardAutomatically(diff, linearMask, startr, endr, totalRounds);
				else
					finalCor += preComputedCor*ComputeTwoRoundRDLCorForwardAutomatically(diff, linearMask, startr, endr, totalRounds);
				

				cout << "Milp result " << finalCor << endl;
			}

			continue;
			
			
			vector<pair<StateTypeBySheet, double>> expansionRes;

			bool isForward = false;
			if (forwardEstimation <= backwardEstimation)
			{
				isForward = true;
				expansionRes = RDLForwardSingleStateDiffThroughChi(diff, RDLLinearMask);
			}
			else
			{
				expansionRes = RDLBackwardSingleStateLinearMaskThroughRevChi(RDLDiff, linearMask);
				isForward = false;
			}

			
			if (isForward)
				cout << __func__ << ": forward direction" << endl;
			else
				cout << __func__ << ": backward direction" << endl;

			

			int rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, endr);
			LaneType rcDiff = rotationalDifferencePerRoundConstants[rcDiffIndex];
			if (isForward)
			{
				rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + 1);
				rcDiff = rotationalDifferencePerRoundConstants[rcDiffIndex];
			}

			for(int j = 0; j < expansionRes.size(); j ++)
			{
				StateTypeBySheet& stateDiffOrLinearMask = expansionRes[j].first;
				double& probOrCor2 = expansionRes[j].second;
				if(isForward)
					rdBase.PropagateThroughLinearLayer(stateDiffOrLinearMask, rcDiff);
				else
					linearBase.PropagateThroughRevLinearLayerWithCor2(stateDiffOrLinearMask, probOrCor2, rcDiff);

				const pair<StateTypeBySheet, StateTypeBySheet>& key = isForward ? make_pair(stateDiffOrLinearMask, linearMask) : make_pair(diff, stateDiffOrLinearMask);

				double dlctThroughNonLinearCor = 1.0;
				const auto& nonLinearLayerInputDiff = key.first;
				const auto& nonLinearLayerOutputMask = key.second;

				bool isValid = true;
				for (int x = 0; x < 4; x++)
				{
					for (int z = 0; z < 32; z++)
					{
						int stateDiffColumn = xoodooBase.GetColumn(nonLinearLayerInputDiff, x, z);
						int stateLinearMaskColumn = xoodooBase.GetColumn(nonLinearLayerOutputMask, x, z);
						dlctThroughNonLinearCor *= dlctCor[stateDiffColumn][stateLinearMaskColumn];
						if (dlctThroughNonLinearCor == 0.0)
						{
							isValid = false;
							break;
						}
					}
					
					if (!isValid)
						break;
				}

				if (!isValid)
					continue;

				#pragma omp critical 
				{
					finalCor += probOrCor2 * preComputedCor * dlctThroughNonLinearCor;
				}
			}

			
		}

		return finalCor;
	}


// =====================================================================================
// ================= Functions below are for the milp models of RDL propagation through xoodoo ==================
// =========


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

	void MilpPropagateThroughRevChi(GRBModel& model, RDStateMILPTypeBySheet& stateMaskVars) const
	{
		MilpPropagateThroughChi(model, stateMaskVars);
	}

	double ComputeTwoRoundRDLCorForwardAutomatically(const StateTypeBySheet& inputStateDiff, const StateTypeBySheet& outputStateLinearMask, int startr, int endr, int totalRounds) const
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_PoolSolutions, 200000000);
		env.set(GRB_IntParam_Threads, milpThreads);

		GRBModel model = GRBModel(env);

		RDStateMILPTypeBySheet stateDiffVars = rdBase.MilpAddStateVars(model);

		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					model.addConstr(stateDiffVars[x][y][z] == xoodooBase.GetBit(inputStateDiff, x, y, z));

		RDStateMILPTypeBySheet inputStateDiffVarsForRound0Chi = stateDiffVars;
		rdBase.MilpPropagateThroughChi(model, stateDiffVars);
		RDStateMILPTypeBySheet outputStateDiffVarsForRound0Chi = stateDiffVars;

		int rotationalDifferenceIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + 1);
		rdBase.MilpPropagateThroughLinearLayer(model, stateDiffVars, rotationalDifferencePerRoundConstants[rotationalDifferenceIndex]);

		GRBLinExpr objExp = 0;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					objExp += stateDiffVars[x][y][z];

		RDStateMILPTypeBySheet inputStateDiffVarsForRound1Chi = stateDiffVars;

		MilpPropagateThroughChi(model, stateDiffVars);
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					model.addConstr(stateDiffVars[x][y][z] == xoodooBase.GetBit(outputStateLinearMask, x, y, z));

		LINStateMILPTypeBySheet outputStateLinearMaskVarsForRound1Chi = stateDiffVars;
		model.setObjective(objExp, GRB_MINIMIZE);

		

		model.optimize();

		int solCount = model.get(GRB_IntAttr_SolCount);
		cout << "Get " << solCount << " two round forward RDL trails." << endl;

		if (solCount == 0)
			return 0.0;
		else
		{
			double RDLCor = 0.0;
			for (int solIndex = 0; solIndex < solCount; solIndex++)
			{
				model.set(GRB_IntParam_SolutionNumber, solIndex);

				StateTypeBySheet inputStateDiffForRound0Chi = rdBase.MilpReadStateVars(inputStateDiffVarsForRound0Chi);
				StateTypeBySheet outputStateDiffForRound0Chi = rdBase.MilpReadStateVars(outputStateDiffVarsForRound0Chi);
				StateTypeBySheet inputStateDiffForRound1Chi = rdBase.MilpReadStateVars(inputStateDiffVarsForRound1Chi);
				StateTypeBySheet outputStateLinearMaskForRound1Chi = linearBase.MilpReadStateVars(outputStateLinearMaskVarsForRound1Chi);
				auto ddtProb = rdBase.getDDTProb();

				double RDLCorOneTrail = 1.0;
				for (int x = 0; x < 4; x++)
					for (int z = 0; z < 32; z++)
					{
						int inputDiffForRound0Chi = xoodooBase.GetColumn(inputStateDiffForRound0Chi, x, z);
						int outputDiffForRound0Chi = xoodooBase.GetColumn(outputStateDiffForRound0Chi, x, z);
						int inputDiffForRound1Chi = xoodooBase.GetColumn(inputStateDiffForRound1Chi, x, z);
						int outputLinearMaskForRound1Chi = xoodooBase.GetColumn(outputStateLinearMaskForRound1Chi, x, z);
						RDLCorOneTrail *= ddtProb[inputDiffForRound0Chi][outputDiffForRound0Chi] * dlctCor[inputDiffForRound1Chi][outputLinearMaskForRound1Chi];
					}

				RDLCor += RDLCorOneTrail;
			}


			return RDLCor;
		}
	}

	double ComputeTwoRoundRDLCorBackwardAutomatically(const StateTypeBySheet& inputStateDiff, const StateTypeBySheet& outputStateLinearMask, int startr, int endr, int totalRounds) const
	{
		GRBEnv env = GRBEnv();
		// env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_PoolSolutions, 200000000);
		env.set(GRB_IntParam_Threads, milpThreads);

		GRBModel model = GRBModel(env);

		LINStateMILPTypeBySheet stateMaskVars = linearBase.MilpAddStateVars(model);

		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					model.addConstr(stateMaskVars[x][y][z] == xoodooBase.GetBit(outputStateLinearMask, x, y, z));

		LINStateMILPTypeBySheet outputStateLinearMaskVarsForRound1Chi = stateMaskVars;
		linearBase.MilpPropagateThroughRevNonLinearLayer(model, stateMaskVars);
		LINStateMILPTypeBySheet inputStateLinearMaskVarsForRound1Chi = stateMaskVars;

		linearBase.MilpPropagateThroughRevLinearLayer(model, stateMaskVars);

		GRBLinExpr objExp = 0;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					objExp += stateMaskVars[x][y][z];

		

		LINStateMILPTypeBySheet outputStateLinearMaskVarsForRound0Chi = stateMaskVars;

		MilpPropagateThroughRevChi(model, stateMaskVars);
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					model.addConstr(stateMaskVars[x][y][z] == xoodooBase.GetBit(inputStateDiff, x, y, z));

		RDStateMILPTypeBySheet inputStateDiffVarsForRound0Chi = stateMaskVars;
		model.setObjective(objExp, GRB_MINIMIZE);



		model.optimize();

		int solCount = model.get(GRB_IntAttr_SolCount);
		cout << "Get " << solCount << " two round backward RDL trails." << endl;

		if (solCount == 0)
			return 0.0;
		else
		{
			double RDLCor = 0.0;

			for (int solIndex = 0; solIndex < solCount; solIndex++)
			{
				StateTypeBySheet inputStateDiffForRound0Chi = rdBase.MilpReadStateVars(inputStateDiffVarsForRound0Chi);
				StateTypeBySheet outputStateLinearMaskForRound0Chi = linearBase.MilpReadStateVars( outputStateLinearMaskVarsForRound0Chi);
				StateTypeBySheet inputStateLinearMaskForRound1Chi = linearBase.MilpReadStateVars( inputStateLinearMaskVarsForRound1Chi);
				StateTypeBySheet outputStateLinearMaskForRound1Chi = linearBase.MilpReadStateVars(outputStateLinearMaskVarsForRound1Chi);

				double RDLCorOneTrail = 1.0;
				auto latCor2 = linearBase.getLATCor2();
				for (int x = 0; x < 4; x++)
					for (int z = 0; z < 32; z++)
					{
						int inputDiffForRound0Chi = xoodooBase.GetColumn(inputStateDiffForRound0Chi, x, z);
						int outputLinearMaskForRound0Chi = xoodooBase.GetColumn(outputStateLinearMaskForRound0Chi, x, z);
						int inputLinearMaskForRound1Chi = xoodooBase.GetColumn(inputStateLinearMaskForRound1Chi, x, z);
						int outputLinearMaskForRound1Chi = xoodooBase.GetColumn(outputStateLinearMaskForRound1Chi, x, z);

						RDLCorOneTrail *= dlctCor[inputDiffForRound0Chi][outputLinearMaskForRound0Chi] * latCor2[inputLinearMaskForRound1Chi][outputLinearMaskForRound1Chi];
					}

				int rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, endr);
				LaneType rcDiff = rotationalDifferencePerRoundConstants[rcDiffIndex];
				linearBase.PropagateThroughLota(inputStateLinearMaskForRound1Chi, RDLCorOneTrail, rcDiff);

				RDLCor += RDLCorOneTrail;
			}


			return RDLCor;
		}


	}

	vector<pair<StateTypeBySheet, double>> RDLForwardSingleStateDiffThroughChiAutomatically(const StateTypeBySheet& inputStateDiff, const StateTypeBySheet& outputStateLinMask, int startr, int endr, int totalRounds) const
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_IntParam_Threads, milpThreads);

		GRBModel model = GRBModel(env);

		RDStateMILPTypeBySheet stateDiffVars = rdBase.MilpAddStateVars(model);

		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					model.addConstr(stateDiffVars[x][y][z] == xoodooBase.GetBit(inputStateDiff, x, y, z));

		rdBase.MilpPropagateThroughChi(model, stateDiffVars);

		RDStateMILPTypeBySheet targetSolVars = stateDiffVars;

		int rotationalDifferenceIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + 1);
		rdBase.MilpPropagateThroughLinearLayer(model, stateDiffVars, rotationalDifferencePerRoundConstants[rotationalDifferenceIndex]);

		MilpPropagateThroughChi(model, stateDiffVars);



		LINStateMILPTypeBySheet stateMaskVars = linearBase.MilpAddStateVars(model);
		for(int x = 0; x < 4;x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					model.addConstr(stateMaskVars[x][y][z] == xoodooBase.GetBit(outputStateLinMask, x, y, z));

		for (int i = endr; i > startr + 1; i--)
		{
			linearBase.MilpPropagateThroughRevNonLinearLayer(model, stateMaskVars);
			linearBase.MilpPropagateThroughRevLinearLayer(model, stateMaskVars);
		}

		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					model.addConstr(stateMaskVars[x][y][z] == stateDiffVars[x][y][z]);

		GRBLinExpr objExp = 0;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					objExp += targetSolVars[x][y][z];

		vector<StateTypeBySheet> targetSolutions;
		FullStateForwardExpandCallback cb = FullStateForwardExpandCallback(xoodooBase, targetSolVars, targetSolutions);
		model.setCallback(&cb);

		model.setObjective(objExp, GRB_MAXIMIZE);

		model.optimize();


		vector<pair<StateTypeBySheet, double>> outputStateDiffAndProbs;
		auto ddtProb = rdBase.getDDTProb();

		for (auto& targetSolution : targetSolutions)
		{
			StateTypeBySheet outputStateDiff = targetSolution;
			double prob = 1.0;
			for (int x = 0; x < 4; x++)
				for (int z = 0; z < 32; z++)
				{
					int inputDiffColumn = xoodooBase.GetColumn(inputStateDiff, x, z);
					int outputDiffColumn = xoodooBase.GetColumn(outputStateDiff, x, z);
					prob *= ddtProb[inputDiffColumn][outputDiffColumn];
				}

			outputStateDiffAndProbs.push_back({ outputStateDiff, prob });
		}

		return outputStateDiffAndProbs;
	}

	vector<pair<StateTypeBySheet, double>> RDLBackwardSingleStateLinearMaskThroughRevChiAutomatically(const StateTypeBySheet& inputStateDiff, const StateTypeBySheet& outputStateLinMask, int startr, int endr, int totalRounds) const
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_IntParam_Threads, milpThreads);

		GRBModel model = GRBModel(env);

		LINStateMILPTypeBySheet stateMaskVars = rdBase.MilpAddStateVars(model);

		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					model.addConstr(stateMaskVars[x][y][z] == xoodooBase.GetBit(outputStateLinMask, x, y, z));

		linearBase.MilpPropagateThroughRevChi(model, stateMaskVars);

		RDStateMILPTypeBySheet targetSolVars = stateMaskVars;

		linearBase.MilpPropagateThroughRevLinearLayer(model, stateMaskVars);

		MilpPropagateThroughRevChi(model, stateMaskVars);



		RDStateMILPTypeBySheet stateDiffVars = rdBase.MilpAddStateVars(model);
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					model.addConstr(stateDiffVars[x][y][z] == xoodooBase.GetBit(inputStateDiff, x, y, z));

		for (int i = startr; i < endr - 1; i++)
		{
			rdBase.MilpPropagateThroughNonLinearLayer(model, stateDiffVars);

			int rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, i+1);
			rdBase.MilpPropagateThroughLinearLayer(model, stateDiffVars, rotationalDifferencePerRoundConstants[rcDiffIndex]);
		}

		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					model.addConstr(stateMaskVars[x][y][z] == stateDiffVars[x][y][z]);

		GRBLinExpr objExp = 0;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 32; z++)
					objExp += targetSolVars[x][y][z];

		vector<StateTypeBySheet> targetSolutions;
		FullStateBackwardExpandCallback cb = FullStateBackwardExpandCallback(xoodooBase, targetSolVars, targetSolutions);
		model.setCallback(&cb);

		model.setObjective(objExp, GRB_MAXIMIZE);

		model.optimize();


		vector<pair<StateTypeBySheet, double>> inputStateMaskAndCor2s;
		auto latCor2 = linearBase.getLATCor2();

		for (auto& targetSolution : targetSolutions)
		{
			StateTypeBySheet inputStateMask = targetSolution;
			double cor2 = 1.0;
			for (int x = 0; x < 4; x++)
				for (int z = 0; z < 32; z++)
				{
					int inputMaskColumn = xoodooBase.GetColumn(inputStateMask, x, z);
					int outputMaskColumn = xoodooBase.GetColumn(outputStateLinMask, x, z);
					cor2 *= latCor2[inputMaskColumn][outputMaskColumn];
				}

			inputStateMaskAndCor2s.push_back({ inputStateMask, cor2 });
		}

		return inputStateMaskAndCor2s;
	}
};

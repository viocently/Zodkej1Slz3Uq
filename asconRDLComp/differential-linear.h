#pragma once
#include<iostream>
#include<vector>
#include<iomanip>
#include<map>
#include<random>
#include<omp.h>
#include<set>
#include"ascon.h"
#include"linear.h"
#include"callback.h"

using namespace std;


using RDLRowType = array<int8_t, 64>;
using RDLStateType = array<RDLRowType, 5>;

RDLStateType RDLStateTypeFromStateType(const StateType& state)
{
	RDLStateType rdlState = { 0 };
	for (int row = 0; row < 5; row++)
		for (int col = 0; col < 64; col++)
			if ((state[row] >> col) & 0x1)
				rdlState[row][col] = 1;
			else
				rdlState[row][col] = 0;

	return rdlState;
}


RDLRowType RotateLeft(const RDLRowType& value, int shift) {
	RDLRowType result;
	for (int i = 0; i < 64; i++) {
		result[i] = value[(i - shift + 64) % 64];
	}
	return result;
}

RDLRowType RotateRight(const RDLRowType& value, int shift) {
	return RotateLeft(value, -shift);
}

class RotationalDifferentialLinearBase
{
private:
	const AsconBase& asconBase;
	const LinearBase& linearBase;
	const RotationalDifferentialBase& rdBase;
	int milpThreads;

	array<array<int, 32>, 32> dlct;
	array<array<double, 32>, 32> dlctCor;
	array<RowType, 12> rotationalDifferencePerRoundConstant;

	// first index is the output linear mask, second index is the mask indicating the input unknown difference bits, third index is the input fixed difference concatenating the input linear mask of the nonfixed bits
	array<array<array<double, 32>, 32>, 32> RDLCor2PerOutputLinearMaskForSbox;

	// first index is the output linear mask, second index is the mask indicating the input unknown difference bits, third index is the input fixed difference, and the vector contains all the possible input linear masks, which are restricted in the positions of unknown difference bits
	array<array<array<vector<int>, 32>, 32>, 32> RDLInputMasksPerOutputMaskForSbox;

	// first index is the input difference, second index is the mask indicating the output unknown linear mask, third index is the output fixed mask concatenating the output difference of the nonfixed bits
	array<array<array<double, 32>, 32>, 32> RDLProbPerInputDifferenceForSbox;

	array<array<array<vector<int>, 32>, 32>, 32> RDLOutputDiffsPerInputDifferenceForSbox;

	// the inequalties for propagation through dlct
	vector<array<int, 11>> milpSboxInequalities;


	void GenerateDLCTForSbox()
	{
		for (int a = 0; a < 32; a++)
			for (int b = 0; b < 32; b++)
				dlct[a][b] = 0;

		for (int xdiff = 0; xdiff < (1 << 5); xdiff++)
			for (int ymask = 0; ymask < (1 << 5); ymask++)
			{
				for (int x0 = 0; x0 < (1 << 5); x0++)
				{
					int x1 = x0 ^ xdiff;
					int y0 = asconBase.Sbox(x0);
					int y1 = asconBase.Sbox(x1);
					int ydiff = y0 ^ y1;
					if (HammingWeight(ydiff & ymask) % 2 == 0)
						dlct[xdiff][ymask]++;
				}
			}

		for (int xdiff = 0; xdiff < (1 << 5); xdiff++)
			for (int ymask = 0; ymask < (1 << 5); ymask++)
			{
				dlct[xdiff][ymask] = 2 * dlct[xdiff][ymask] - (1 << 5);
				dlctCor[xdiff][ymask] = dlct[xdiff][ymask] / 32.0;
			}
	}

	void GenerateRDLProbPerInputDifferenceForSbox()
	{
		vector<vector<int>> outputDiffPerInputDiffForSbox = rdBase.GetOutputDiffsPerInputDiffForSbox();
		array<array<double, 32>, 32> ddtProb = rdBase.GetDDTProb();

		for (int inputDiff = 0; inputDiff < (1 << 5); inputDiff++)
		{
			vector<int> possibleOutputDiffs = outputDiffPerInputDiffForSbox[inputDiff];
			for (int outputUnknownLinearMaskBitsMask = 0; outputUnknownLinearMaskBitsMask < (1 << 5); outputUnknownLinearMaskBitsMask++)
			{
				RDLProbPerInputDifferenceForSbox[inputDiff][outputUnknownLinearMaskBitsMask].fill(0.0);
				int outputFixedLinearMaskBitsMask = (1 << 5) - 1 - outputUnknownLinearMaskBitsMask;
				vector<int> outputFixedLinearMaskBitsPos;
				for (int i = 0; i < 5; i++)
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
							RDLProbPerInputDifferenceForSbox[inputDiff][outputUnknownLinearMaskBitsMask][outputFixedLinearMaskBits | outputDiffForUnknownLinearMaskBits] += -prob;
						else
							RDLProbPerInputDifferenceForSbox[inputDiff][outputUnknownLinearMaskBitsMask][outputFixedLinearMaskBits | outputDiffForUnknownLinearMaskBits] += prob;
					}
				}
			}
		}
	}

	void GenerateRDLOutputDiffsPerInputDifferenceForSbox()
	{
		for (int inputDiff = 0; inputDiff < (1 << 5); inputDiff++)
		{
			for (int outputUnknownLinearMaskBitsMask = 0; outputUnknownLinearMaskBitsMask < (1 << 5); outputUnknownLinearMaskBitsMask++)
			{
				int outputFixedLinearMaskBitsMask = (1 << 5) - 1 - outputUnknownLinearMaskBitsMask;

				array<double, 32> RDLProbs = RDLProbPerInputDifferenceForSbox[inputDiff][outputUnknownLinearMaskBitsMask];
				for (int v = 0; v < (1 << 5); v++)
				{
					if (RDLProbs[v] != 0.0)
					{
						int outputFixedLinearMaskBits = v & outputFixedLinearMaskBitsMask;
						int outputDiffForUnknownLinearMaskBits = v & outputUnknownLinearMaskBitsMask;
						RDLOutputDiffsPerInputDifferenceForSbox[inputDiff][outputUnknownLinearMaskBitsMask][outputFixedLinearMaskBits].push_back(outputDiffForUnknownLinearMaskBits);
					}
				}
			}
		}
	}

	void GenerateRDLCor2PerOutputLinearMaskForSbox()
	{
		vector<vector<int>> inputLinearMasksPerOutputLinearMaskForChi = linearBase.GetInputMasksPerOutputMaskForSbox();
		array<array<double, 32>, 32> latCor2 = linearBase.GetLATCor2();

		for (int outputLinearMask = 0; outputLinearMask < (1 << 5); outputLinearMask++)
		{
			// cout << "Current output linear mask " << outputLinearMask << endl;
			vector<int> possibleInputLinearMasks = inputLinearMasksPerOutputLinearMaskForChi[outputLinearMask];
			for (int inputUnknownDiffBitsMask = 0; inputUnknownDiffBitsMask < (1 << 5); inputUnknownDiffBitsMask++)
			{
				// cout << "Current unknown difference bits mask " << inputUnknownDiffBitsMask << endl;
				RDLCor2PerOutputLinearMaskForSbox[outputLinearMask][inputUnknownDiffBitsMask].fill(0.0);
				int inputFixedDiffBitsMask = (1 << 5) - 1 - inputUnknownDiffBitsMask;
				vector<int> inputFixedDiffBitsPos;
				for (int i = 0; i < 5; i++)
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
							RDLCor2PerOutputLinearMaskForSbox[outputLinearMask][inputUnknownDiffBitsMask][inputFixedDiff | inputLinearMaskForUnknownDiffbits] += -cor2;
						else
							RDLCor2PerOutputLinearMaskForSbox[outputLinearMask][inputUnknownDiffBitsMask][inputFixedDiff | inputLinearMaskForUnknownDiffbits] += cor2;
					}
				}

			}
		}
	}

	void GenerateRDLInputMasksPerOutputMaskForSbox()
	{
		for (int outputLinearMask = 0; outputLinearMask < (1 << 5); outputLinearMask++)
		{
			for (int inputUnknownDiffBitsMask = 0; inputUnknownDiffBitsMask < (1 << 5); inputUnknownDiffBitsMask++)
			{
				int inputFixedDiffBitsMask = (1 << 5) - 1 - inputUnknownDiffBitsMask;

				array<double, 32> RDLCor2s = RDLCor2PerOutputLinearMaskForSbox[outputLinearMask][inputUnknownDiffBitsMask];
				for (int v = 0; v < (1 << 5); v++)
				{
					if (RDLCor2s[v] != 0.0)
					{
						int inputFixedDiffBits = v & inputFixedDiffBitsMask;
						int inputLinearMaskForUnknownDiffBits = v & inputUnknownDiffBitsMask;
						RDLInputMasksPerOutputMaskForSbox[outputLinearMask][inputUnknownDiffBitsMask][inputFixedDiffBits].push_back(inputLinearMaskForUnknownDiffBits);
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
	RotationalDifferentialLinearBase(const AsconBase& asconBase, const LinearBase& linearBase, const RotationalDifferentialBase& rdBase, int milpThreads = 2) : asconBase(asconBase), linearBase(linearBase), rdBase(rdBase), milpThreads(milpThreads)
	{
		GenerateRDLCor2PerOutputLinearMaskForSbox();
		GenerateRDLProbPerInputDifferenceForSbox();
		GenerateRDLInputMasksPerOutputMaskForSbox();
		GenerateRDLOutputDiffsPerInputDifferenceForSbox();
		GenerateDLCTForSbox();

		rotationalDifferencePerRoundConstant = rdBase.GetRotationalDifferencePerRoundConstant();

		milpSboxInequalities = {{24, 1, -13, 7, -8, 5, 3, -4, -2, -7, -11}, {11, 1, -3, -5, 5, -2, -2, -3, 3, -6, 0}, {7, 1, -1, 0, -4, -2, -1, -3, 2, 4, -3}, {3, -1, -1, 0, -1, 0, -1, 2, -2, 2, 0}, {7, -7, 4, 2, 2, 3, 1, -4, -2, 0, -1}, {5, -1, 4, -1, -3, 0, 2, 3, -1, -4, 2}, {7, 3, 2, -2, 2, -5, -1, -4, 0, -1, 0}, {2, 0, -2, 1, 0, 0, 0, 1, -1, -1, 0}, {11, -5, -1, 0, -3, 4, -9, -2, 4, 0, 6}, {3, 1, -2, -1, 2, 1, 0, 0, -2, 0, -1}, {3, 1, 0, -1, -2, 1, 3, 0, 1, 0, -3}, {4, 2, 4, 2, -2, -2, -1, -1, 0, -2, 1}, {2, 0, 0, -1, -1, 0, -1, 0, 1, 2, -1}, {2, 2, 0, 1, 0, -2, 1, -1, 0, 0, -1}, {2, -2, 1, 0, 0, 1, 1, 0, -2, 1, 0}, {3, 1, 1, -2, 1, -1, 0, -1, 0, -2, 1}, {2, 2, 1, 0, -1, -1, -2, 0, 0, 1, 1}, {3, -1, -2, 1, 1, 0, 0, 1, -1, -2, 0}, {3, 1, -2, 0, 0, -1, -1, -2, 2, 3, 1}, {1, 0, 0, 1, -1, 0, -1, 1, 1, 0, 0}, {6, -1, 2, -4, 0, -1, 2, -2, 4, -1, -1}, {1, -1, 0, 1, 0, 1, 1, 0, 0, 0, -1}, {1, 0, -1, 0, 0, 1, 1, 0, -1, 0, 0}, {3, -1, 0, -1, 0, -1, 0, 0, 0, 1, -1}};


	}

	array<array<double, 32>, 32> GetDLCTCor() const
	{
		return dlctCor;
	}

	void OutputRDLState(ostream& os, const RDLStateType& state) const
	{
		for (int row = 0; row < 5; row++)
		{
			for (int col = 0; col < 64; col++)
			{
				os << setw(2) << state[row][col];
			}
			os << endl;
		}
	}

	// =====================================================================================
	// ================= Functions below are for differential and linear propagation in the middle part ==================
	// =====================================================================================

	void RDLDiffPropagateThroughRotateRows(RDLStateType& state) const
	{
		for (int row = 0; row < 5; row++)
		{
			auto tmp0 = RotateRight(state[row], asconBase.GetRotationNum(row, 0));
			auto tmp1 = RotateRight(state[row], asconBase.GetRotationNum(row, 1));
			for (int col = 0; col < 64; col++)
			{
				state[row][col] = xor3(state[row][col], tmp0[col], tmp1[col]);
			}
		}
	}

	void RDLLinPropagateThroughInvRotateRows(RDLStateType& state) const
	{
		for (int row = 0; row < 5; row++)
		{
			auto tmp0 = RotateLeft(state[row], asconBase.GetRotationNum(row, 0));
			auto tmp1 = RotateLeft(state[row], asconBase.GetRotationNum(row, 1));
			for (int col = 0; col < 64; col++)
			{
				state[row][col] = xor3(state[row][col], tmp0[col], tmp1[col]);
			}
		}
	}

	void RDLDiffPropagateThroughAddRoundConstant(RDLStateType& state, const RowType rcDiff) const
	{
		for (int col = 0; col < 64; col++)
		{
			state[2][col] = xor2(state[2][col], (rcDiff >> col) & 1ull);
		}
	}

	void RDLDiffPropagateThroughSbox(int& y0, int& y1, int& y2, int& y3, int & y4) const
	{
		if (y4 == 0 && y3 == 0 && y2 == 0 && y1 == 0 && y0 == 0) { 
			y4 = 0; y3 = 0; y2 = 0; y1 = 0; y0 = 0; 
		}
		else if (y4 == 0 && y3 == 0 && y2 == 0 && y1 == 0 && y0 == 1) { 
			y4 = -1; y3 = 1; y2 = -1; y1 = -1; y0 = -1; 
		}
		else if (y4 == 0 && y3 == 0 && y2 == 0 && y1 == 1 && y0 == 0) {
			y4 = 1; y3 = -1; y2 = -1; y1 = -1; y0 = 1;
		}
		else if (y4 == 0 && y3 == 0 && y2 == 0 && y1 == 1 && y0 == 1) {
			y4 = -1; y3 = -1; y2 = -1; y1 = 0; y0 = -1;
		}
		else if (y4 == 0 && y3 == 0 && y2 == 1 && y1 == 0 && y0 == 0) {
			y4 = -1; y3 = -1; y2 = 1; y1 = 1; y0 = 0;
		}
		else if (y4 == 0 && y3 == 0 && y2 == 1 && y1 == 0 && y0 == 1) {
			y4 = 1; y3 = -1; y2 = -1; y1 = -1; y0 = -1;
		}
		else if (y4 == 0 && y3 == 0 && y2 == 1 && y1 == 1 && y0 == 0) {
			y4 = -1;  y3 = -1; y2 = -1; y1 = -1; y0 = 1;
		}
		else if (y4 == 0 && y3 == 0 && y2 == 1 && y1 == 1 && y0 == 1) {
			y4 = 0; y3 = -1; y2 = -1; y1 = 1; y0 = -1;
		}
		else if (y4 == 0 && y3 == 0 && y2 == -1 && y1 == 0 && y0 == 0) {
			y4 = -1; y3 = -1; y2 = -1; y1 = -1; y0 = 0;
		}
		else if (y4 == 0 && y3 == 0 && y2 == -1 && y1 == 1 && y0 == 0) {
			y4 = -1; y3 = -1; y2 = -1; y1 = -1; y0 = 1;
		}
		else if (y4 == 0 && y3 == 1 && y2 == 0 && y1 == 0 && y0 == 0) {
			y4 = -1; y3 = -1; y2 = 1; y1 = 1; y0 = -1;
		}
		else if (y4 == 0 && y3 == 1 && y2 == 0 && y1 == 1 && y0 == 1) {
			y4 = -1; y3 = -1; y2 = -1; y1 = 1; y0 = -1;
		}
		else if (y4 == 0 && y3 == 1 && y2 == 1 && y1 == 0 && y0 == 0) {
			y4 = -1; y3 = -1; y2 = 0; y1 = 0; y0 = -1;
		}
		else if (y4 == 0 && y3 == 1 && y2 == 1 && y1 == 1 && y0 == 0) {
			y4 = -1; y3 = 0; y2 = -1; y1 = -1; y0 = -1;
		}
		else if(y4 == 0 && y3 == 1 && y2 == 1 && y1 == 1 && y0 == 1) {
			y4 = -1 ; y3 = 1 ; y2 = -1 ; y1 = 0 ; y0 = -1;
		}
		else if(y4 == 1 && y3 == 0 && y2 == 0 && y1 == 0 && y0 == 0) {
			y4 = -1 ; y3 = 1 ; y2 = 0 ; y1 = -1 ; y0 = -1;
		}
		else if(y4 == 1 && y3 == 0 && y2 == 0 && y1 == 0 && y0 == 1) {
			y4 = 1 ; y3 = 0 ; y2 = -1 ; y1 = -1 ; y0 = 1;
		}
		else if(y4 == 1 && y3 == 0 && y2 == 0 && y1 == 1 && y0 == 1) {
			y4 = 0 ; y3 = -1 ; y2 = -1 ; y1 = -1 ; y0 = 0;
		}
		else if(y4 == 1 && y3 == 0 && y2 == 1 && y1 == 0 && y0 == 0) {
			y4 = 0 ; y3 = -1 ; y2 = 1 ; y1 = -1 ; y0 = -1;
		}
		else if(y4 == 1 && y3 == 0 && y2 == 1 && y1 == 0 && y0 == 1) {
			y4 = -1 ; y3 = -1 ; y2 = -1 ; y1 = -1 ; y0 = 1;
		}
		else if(y4 == 1 && y3 == 0 && y2 == 1 && y1 == 1 && y0 == 0) {
			y4 = 1 ; y3 = -1 ; y2 = -1 ; y1 = -1 ; y0 = -1;
		}
		else if(y4 == 1 && y3 == 0 && y2 == 1 && y1 == 1 && y0 == 1) {
			y4 = -1 ; y3 = -1 ; y2 = -1 ; y1 = -1 ; y0 = 0;
		}
		else if(y4 == 1 && y3 == 0 && y2 == -1 && y1 == 0 && y0 == 1) {
			y4 = -1 ; y3 = -1 ; y2 = -1 ; y1 = -1 ; y0 = 1;
		}
		else if(y4 == 1 && y3 == 0 && y2 == -1 && y1 == 1 && y0 == 1) {
			y4 = -1 ; y3 = -1 ; y2 = -1 ; y1 = -1 ; y0 = 0;
		}
		else if(y4 == 1 && y3 == 1 && y2 == 0 && y1 == 0 && y0 == 0) {
			y4 = -1 ; y3 = -1 ; y2 = 1 ; y1 = -1 ; y0 = -1;
		}
		else if(y4 == 1 && y3 == 1 && y2 == 1 && y1 == 0 && y0 == 0) {
			y4 = -1 ; y3 = -1 ; y2 = 0 ; y1 = -1 ; y0 = -1;
		}
		else if(y4 == 1 && y3 == 1 && y2 == 1 && y1 == 1 && y0 == 0) {
			y4 = -1 ; y3 = 1 ; y2 = -1 ; y1 = -1 ; y0 = -1;
		}
		else if(y4 == 1 && y3 == 1 && y2 == 1 && y1 == 1 && y0 == 1) {
			y4 = -1 ; y3 = 0 ; y2 = -1 ; y1 = -1 ; y0 = -1;
		}
		else if(y4 == -1 && y3 == 0 && y2 == 0 && y1 == 0 && y0 == 0) {
			y4 = -1 ; y3 = -1 ; y2 = 0 ; y1 = -1 ; y0 = -1;
		}
		else if(y4 == -1 && y3 == 0 && y2 == 1 && y1 == 0 && y0 == 0) {
			y4 = -1 ; y3 = -1 ; y2 = 1 ; y1 = -1 ; y0 = -1;
		}
		else if(y4 == -1 && y3 == 1 && y2 == 0 && y1 == 0 && y0 == 0) {
			y4 = -1 ; y3 = -1 ; y2 = 1 ; y1 = -1 ; y0 = -1;
		}
		else if (y4 == -1 && y3 == 1 && y2 == 1 && y1 == 0 && y0 == 0) {
			y4 = -1; y3 = -1; y2 = 0; y1 = -1; y0 = -1;
		}
		else { y4 = -1; y3 = -1; y2 = -1;  y1 = -1; y0 = -1; }
	}

	void RDLDiffPropagateThroughSboxes(RDLStateType& state) const
	{
		for (int col = 0; col < 64; col++)
		{
			int y0 = state[0][col];
			int y1 = state[1][col];
			int y2 = state[2][col];
			int y3 = state[3][col];
			int y4 = state[4][col];
			RDLDiffPropagateThroughSbox(y4, y3, y2, y1, y0);
			state[0][col] = y0;
			state[1][col] = y1;
			state[2][col] = y2;
			state[3][col] = y3;
			state[4][col] = y4;
		}
	}


	void RDLLinPropagateThroughInvSbox(int& y0, int& y1, int& y2, int& y3, int& y4) const
	{
		if (y4 == 0 && y3 == 0 && y2 == 0 && y1 == 0 && y0 == 0) { y4 = 0; y3 = 0; y2 = 0; y1 = 0; y0 = 0; }
		else if (y4 == 0 && y3 == 0 && y2 == 0 && y1 == 0 && y0 == 1) { y4 = -1; y3 = -1; y2 = 0; y1 = 1; y0 = -1; }
		else if (y4 == 0 && y3 == 0 && y2 == 0 && y1 == 0 && y0 == -1) { y4 = -1; y3 = -1; y2 = 0; y1 = -1; y0 = -1; }
		else if (y4 == 0 && y3 == 0 && y2 == 0 && y1 == 1 && y0 == 0) { y4 = -1; y3 = 1; y2 = 1; y1 = -1; y0 = -1; }
		else if (y4 == 0 && y3 == 0 && y2 == 0 && y1 == 1 && y0 == 1) { y4 = -1; y3 = -1; y2 = 1; y1 = -1; y0 = -1; }
		else if (y4 == 0 && y3 == 0 && y2 == 0 && y1 == 1 && y0 == -1) { y4 = -1; y3 = -1; y2 = 1; y1 = -1; y0 = -1; }
		else if (y4 == 0 && y3 == 0 && y2 == 1 && y1 == 0 && y0 == 0) { y4 = 0; y3 = 1; y2 = 1; y1 = -1; y0 = -1; }
		else if (y4 == 0 && y3 == 0 && y2 == 1 && y1 == 0 && y0 == 1) { y4 = -1; y3 = -1; y2 = 1; y1 = -1; y0 = -1; }
		else if (y4 == 0 && y3 == 0 && y2 == 1 && y1 == 0 && y0 == -1) { y4 = -1; y3 = -1; y2 = 1; y1 = -1; y0 = -1; }
		else if (y4 == 0 && y3 == 0 && y2 == 1 && y1 == 1 && y0 == 0) { y4 = -1; y3 = 0; y2 = 0; y1 = -1; y0 = -1; }
		else if (y4 == 0 && y3 == 0 && y2 == 1 && y1 == 1 && y0 == 1) { y4 = -1; y3 = -1; y2 = 0; y1 = -1; y0 = -1; }
		else if (y4 == 0 && y3 == 0 && y2 == 1 && y1 == 1 && y0 == -1) { y4 = -1; y3 = -1; y2 = 0; y1 = -1; y0 = -1; }
		else if (y4 == 0 && y3 == 0 && y2 == -1 && y1 == 0 && y0 == 0) { y4 = 0; y3 = -1; y2 = -1; y1 = -1; y0 = -1; }
		else if (y4 == 0 && y3 == 1 && y2 == 0 && y1 == 0 && y0 == 0) { y4 = 1; y3 = -1; y2 = -1; y1 = -1; y0 = 1; }
		else if (y4 == 0 && y3 == 1 && y2 == 1 && y1 == 0 && y0 == 0) { y4 = 1; y3 = -1; y2 = -1; y1 = -1; y0 = -1; }
		else if (y4 == 0 && y3 == 1 && y2 == -1 && y1 == 0 && y0 == 0) { y4 = 1; y3 = -1; y2 = -1; y1 = -1; y0 = -1; }
		else if (y4 == 1 && y3 == 0 && y2 == 0 && y1 == 0 && y0 == 0) { y4 = -1; y3 = -1; y2 = -1; y1 = 1; y0 = -1; }
		else if (y4 == 1 && y3 == 0 && y2 == 0 && y1 == 0 && y0 == 1) { y4 = 1; y3 = -1; y2 = -1; y1 = 0; y0 = 1; }
		else if (y4 == 1 && y3 == 0 && y2 == 1 && y1 == 0 && y0 == 1) { y4 = 1; y3 = -1; y2 = -1; y1 = -1; y0 = -1; }
		else if (y4 == 1 && y3 == 0 && y2 == -1 && y1 == 0 && y0 == 1) { y4 = 1; y3 = -1; y2 = -1; y1 = -1; y0 = -1; }
		else if (y4 == 1 && y3 == 1 && y2 == 0 && y1 == 0 && y0 == 1) { y4 = 0; y3 = -1; y2 = -1; y1 = -1; y0 = 0; }
		else if (y4 == 1 && y3 == 1 && y2 == 1 && y1 == 0 && y0 == 1) { y4 = 0; y3 = -1; y2 = -1; y1 = -1; y0 = -1; }
		else if (y4 == 1 && y3 == 1 && y2 == -1 && y1 == 0 && y0 == 1) { y4 = 0; y3 = -1; y2 = -1; y1 = -1; y0 = -1; }
		else {y4 = -1; y3 = -1; y2 = -1; y1 = -1; y0 = -1;}
	}

	void RDLLinPropagateThroughInvSboxes(RDLStateType& state) const
	{
		for (int col = 0; col < 64; col++)
		{
			int y0 = state[0][col];
			int y1 = state[1][col];
			int y2 = state[2][col];
			int y3 = state[3][col];
			int y4 = state[4][col];
			RDLLinPropagateThroughInvSbox(y4, y3, y2, y1, y0);
			state[0][col] = y0;
			state[1][col] = y1;
			state[2][col] = y2;
			state[3][col] = y3;
			state[4][col] = y4;
		}
	}

	void RDLDiffPropagateThroughLinearLayer(RDLStateType& state, const RowType rcDiff) const
	{
		RDLDiffPropagateThroughRotateRows(state);
		RDLDiffPropagateThroughAddRoundConstant(state, rcDiff);
	}

	void RDLLinPropagateThroughInvLinearLayer(RDLStateType& state) const
	{
		RDLLinPropagateThroughInvRotateRows(state);
	}

	void RDLDiffPropagateThroughNonLinearLayer(RDLStateType& state) const
	{
		RDLDiffPropagateThroughSboxes(state);
	}

	void RDLLinPropagateThroughInvNonLinearLayer(RDLStateType& state) const
	{
		RDLLinPropagateThroughInvSboxes(state);
	}


	// =====================================================================================
	// ================= Functions below are for computing the correlation of RDL distinguisher in the middle part ==================
	// =====================================================================================

	double ExperimentalComputation(const StateType& inputStateDiff, const StateType& outputStateMask, int startr, int endr) const
	{
		random_device rd;
		mt19937 gen(rd());
		uniform_int_distribution<uint64_t> distribution(0, numeric_limits<uint64_t>::max());

		long long totalTestNum = 1ll << 27;
		long long zeroCount = 0;

		int rot = rdBase.GetRotationNumber();

		for (long long testIndex = 0; testIndex < totalTestNum; testIndex++)
		{
			StateType state0;
			for (int row = 0; row < 5; row++)
				state0[row] = distribution(gen);

			StateType state1;
			for (int row = 0; row < 5; row++)
			{
				state1[row] = state0[row] ^ inputStateDiff[row];
				state1[row] = RotateRight(state1[row], rot);
			}

			asconBase.ReducedRoundsForVerification(state0, startr, endr);
			asconBase.ReducedRoundsForVerification(state1, startr, endr);

			StateType outputStateDiff;
			for (int row = 0; row < 5; row++)
				outputStateDiff[row] = state0[row] ^ RotateLeft(state1[row], rot);


			int outputHW = 0;
			for (int row = 0; row < 5; row++)
				outputHW += HammingWeight(outputStateDiff[row] & outputStateMask[row]);

			if (outputHW % 2 == 0)
				zeroCount++;
		}

		double zeroProb = zeroCount / (double)totalTestNum;

		return 2 * zeroProb - 1;

	}


	double ComputeOneRoundRDLCor(const StateType& inputStateDiff, const StateType& outputStateMask, int r) const
	{
		double RDLCor = 1.0;
		StateType sboxesOutputMask = outputStateMask;
		linearBase.PropagateThroughInvLinearLayerWithCor2(sboxesOutputMask, RDLCor, rotationalDifferencePerRoundConstant[r + 1]);

		for (int col = 0; col < 64; col++)
		{
			uint8_t inputColumnDiff = asconBase.GetColumn(inputStateDiff, col);
			uint8_t outputColumnMask = asconBase.GetColumn(sboxesOutputMask, col);
			RDLCor *= dlctCor[inputColumnDiff][outputColumnMask];
		}

		return RDLCor;
	}

	// the inputStateDiff corresponds the input to the round function, and the outputRDLMask correponds to the output of the Sbox layer
	double RDLForwardSingleStateDiffThroughSboxEstimation(const StateType& inputStateDiff, const RDLStateType& outputRDLMask) const
	{
		array<int, 64> outputUnknownLinearMaskBitsMask = { 0 };
		array<int, 64> outputFixedLinearMaskBits = { 0 };

		for (int col = 0; col < 64; col++)
		{
			for (int row = 0; row < 5; row++)
				if (outputRDLMask[row][col] == -1)
					outputUnknownLinearMaskBitsMask[col] |= (1 << (4-row));
				else if (outputRDLMask[row][col] == 1)
					outputFixedLinearMaskBits[col] |= (1 << (4-row));
		}


		double estimatedNumLog2 = 0;
		for (int col = 0; col < 64; col++)
		{
			int inputDiff = asconBase.GetColumn(inputStateDiff, col);

			const vector<int>& possibleOutputDiffs = RDLOutputDiffsPerInputDifferenceForSbox[inputDiff][outputUnknownLinearMaskBitsMask[col]][outputFixedLinearMaskBits[col]];
			// cout << inputDiff << "," << outputUnknownLinearMaskBitsMask[col] << "," << outputFixedLinearMaskBits[col] << " => ";

			if (possibleOutputDiffs.size() == 0)
				return -1;

			if (inputDiff != 0)
				estimatedNumLog2 += log2(possibleOutputDiffs.size());
		}

		return estimatedNumLog2;
	}

	// the inputRDLDiff corresponds to the input of the round function (also the Sbox layer), and the outputMask corresponds to the output of the round function
	double RDLBackwardSingleStateMaskThroughInvSboxEstimation(const RDLStateType& inputRDLDiff, const StateType& outputMask) const
	{
		StateType sboxesOutputMask = outputMask;
		linearBase.PropagateThroughInvLinearLayer(sboxesOutputMask);

		array<int, 64> inputUnknownDiffBitsMask = { 0 };
		array<int, 64> inputFixedDiffBits = { 0 };
		for (int col = 0; col < 64; col++)
		{
			for (int row = 0; row < 5; row++)
				if (inputRDLDiff[row][col] == -1)
					inputUnknownDiffBitsMask[col] |= (1 << (4-row));
				else if (inputRDLDiff[row][col] == 1)
					inputFixedDiffBits[col] |= (1 << (4-row));
		}

		double estimatedNumLog2 = 0;
		for (int col = 0; col < 64; col++)
		{
			int outputMaskCol = asconBase.GetColumn(sboxesOutputMask, col);
			const vector<int>& possibleInputLinearMasks = RDLInputMasksPerOutputMaskForSbox[outputMaskCol][inputUnknownDiffBitsMask[col]][inputFixedDiffBits[col]];
			if (possibleInputLinearMasks.size() == 0)
				return -1;
			if (outputMaskCol != 0)
				estimatedNumLog2 += log2(possibleInputLinearMasks.size());
		}

		return estimatedNumLog2;
	}

	void RDLBackwardExpansionByOneRound(vector<pair<StateType, StateType>>& diffAndMaskPairs, vector<double>& preComputedCors, vector<pair<RDLStateType, RDLStateType>>& rdlMaskAndDiffPairs, int& startr, int& endr, bool isCallBack) const
	{
		map<pair<StateType, StateType>, double> m0;
		map<pair<StateType, StateType>, pair<RDLStateType, RDLStateType>> m1;
		int pairsNum = diffAndMaskPairs.size();
		cout << __func__ << ": to expand " << pairsNum << " pairs" << endl;
		cout << "Start round: " << startr << endl;
		cout << "End round: " << endr << endl;
		int threadsNum = omp_get_max_threads();

#pragma omp parallel for num_threads(threadsNum)
		for (int i = 0; i < pairsNum; i++)
		{
			// cout << "Processing pair " << i + 1 << " / " << pairsNum << "\n";

			const pair<StateType, StateType> diffAndMaskPair = diffAndMaskPairs[i];
			const StateType& diff = diffAndMaskPair.first;
			const StateType& mask = diffAndMaskPair.second;
			const double& preComputedCor = preComputedCors[i];
			const RDLStateType& rdlMask = rdlMaskAndDiffPairs[i].first;
			const RDLStateType& rdlDiff = rdlMaskAndDiffPairs[i].second;
			RDLStateType newRDLDiff = RDLStateTypeFromStateType(diff);
			for (int j = startr; j < endr - 1; j++)
			{
				RDLDiffPropagateThroughNonLinearLayer(newRDLDiff);
				RDLDiffPropagateThroughLinearLayer(newRDLDiff, rotationalDifferencePerRoundConstant[j + 1]);
			}


			int backwardEstimation = RDLBackwardSingleStateMaskThroughInvSboxEstimation(rdlDiff, mask);
			if (backwardEstimation == -1)
				continue;
			else
			{
				;
			}

			vector<pair<StateType, double>> newMaskAndCor2s;

			if (isCallBack)
			{
				int midr = ((startr + endr) / 2) > (endr - 1) ? (endr - 1) : ((startr + endr) / 2);
				newMaskAndCor2s = RDLBackwardSinlgeStateMaskThroughOneRoundAutomatically(diff, mask, startr, endr, midr);
			}
			else
				newMaskAndCor2s = RDLBackwardSingleStateMaskThroughOneRound(rdlDiff, mask);

#pragma omp critical 
			{
				for (auto& maskAndCor2 : newMaskAndCor2s)
				{
					StateType stateMask = maskAndCor2.first;
					double cor2 = maskAndCor2.second;
					if (m0.find({ diff, stateMask }) != m0.end())
						m0[{diff, stateMask}] += cor2 * preComputedCor;
					else
					{
						m0[{diff, stateMask}] = cor2 * preComputedCor;
						RDLStateType newRDLMask = RDLStateTypeFromStateType(stateMask);
						for (int j = endr - 1; j > startr; j--)
						{
							RDLLinPropagateThroughInvLinearLayer(newRDLMask);
							RDLLinPropagateThroughInvNonLinearLayer(newRDLMask);
						}

						RDLLinPropagateThroughInvLinearLayer(newRDLMask);

						m1[{diff, stateMask}] = { newRDLMask, newRDLDiff };
					}
				}
			}

		}

		diffAndMaskPairs.clear();
		preComputedCors.clear();
		rdlMaskAndDiffPairs.clear();
		for (auto& it : m0)
		{
			diffAndMaskPairs.push_back(it.first);
			preComputedCors.push_back(it.second);
			rdlMaskAndDiffPairs.push_back(m1[it.first]);
		}

		endr -= 1;
	}


	void RDLForwardExpansionByOneRound(vector<pair<StateType, StateType>>& diffAndMaskPairs, vector<double>& preComputedCors, vector<pair<RDLStateType, RDLStateType>>& rdlMaskAndDiffPairs, int& startr, int& endr, bool isCallBack) const
	{
		map<pair<StateType, StateType>, double> m0;
		map<pair<StateType, StateType>, pair<RDLStateType, RDLStateType>> m1;

		int pairsNum = diffAndMaskPairs.size();
		cout << __func__ << ": to expand " << pairsNum << " pairs" << endl;
		cout << "Start round: " << startr << endl;
		cout << "End round: " << endr << endl;

		int threadsNum = omp_get_max_threads();
#pragma omp parallel for num_threads(threadsNum)
		for (int i = 0; i < pairsNum; i++)
		{

			// cout << "Processing pair " << i + 1 << " / " << pairsNum << "\n";
			const pair<StateType, StateType> diffAndMaskPair = diffAndMaskPairs[i];
			const StateType& diff = diffAndMaskPair.first;
			const StateType& mask = diffAndMaskPair.second;
			const double& preComputedCor = preComputedCors[i];
			const RDLStateType& rdlMask = rdlMaskAndDiffPairs[i].first;
			const RDLStateType& rdlDiff = rdlMaskAndDiffPairs[i].second;

			RDLStateType newRDLMask = RDLStateTypeFromStateType(mask);
			for (int j = endr; j > startr + 1; j--)
			{
				RDLLinPropagateThroughInvLinearLayer(newRDLMask);
				RDLLinPropagateThroughInvNonLinearLayer(newRDLMask);
			}
			RDLLinPropagateThroughInvLinearLayer(newRDLMask);


			int forwardEstimation = RDLForwardSingleStateDiffThroughSboxEstimation(diff, rdlMask);
			
			if (forwardEstimation == -1)
				continue;
			else
			{
				;
			}
			

			vector<pair<StateType, double>> newDiffAndProbs;

			if (isCallBack)
			{
				int midr = ((startr + endr) / 2) < (startr + 1) ? (startr + 1) : ((startr + endr) / 2);
				newDiffAndProbs = RDLForwardSingleStateDiffThroughOneRoundAutomatically(diff, mask, startr, endr, midr);
			}
			else
				newDiffAndProbs = RDLForwardSingleStateDiffThroughOneRound(diff, rdlMask, startr);

#pragma omp critical 
			{
				for (auto& diffAndProb : newDiffAndProbs)
				{
					StateType stateDiff = diffAndProb.first;
					double prob = diffAndProb.second;


					if (m0.find({ stateDiff, mask }) != m0.end())
						m0[{stateDiff, mask}] += prob * preComputedCor;
					else
					{
						m0[{stateDiff, mask}] = prob * preComputedCor;

						RDLStateType newRDLDiff = RDLStateTypeFromStateType(stateDiff);
						for (int j = startr + 1; j < endr; j++)
						{
							RDLDiffPropagateThroughNonLinearLayer(newRDLDiff);
							RDLDiffPropagateThroughLinearLayer(newRDLDiff, rotationalDifferencePerRoundConstant[j + 1]);
						}

						m1[{stateDiff, mask}] = { newRDLMask, newRDLDiff };

					}
				}
			}
		}

		diffAndMaskPairs.clear();
		preComputedCors.clear();
		rdlMaskAndDiffPairs.clear();

		for (auto& it : m0)
		{
			diffAndMaskPairs.push_back(it.first);
			preComputedCors.push_back(it.second);
			rdlMaskAndDiffPairs.push_back(m1[it.first]);
		}

		startr += 1;

	}

	void FilterLowCorPairs(vector<pair<StateType, StateType>>& diffAndMaskPairs, vector<double>& preComputedCors, vector<pair<RDLStateType, RDLStateType>>& rdlMaskAndDiffPairs) const
	{
		cout << __func__ << ": before filtering, pairs num " << diffAndMaskPairs.size() << endl;

		double bound = 0;
		if (diffAndMaskPairs.size() < 10000)
			bound = 0.0002;
		else if (diffAndMaskPairs.size() < 100000)
			bound = 0.00002;
		else if (diffAndMaskPairs.size() < 1000000)
			bound = 0.000002;
		else
			bound = 0.0000002;

		vector<pair<StateType, StateType>> newDiffAndMaskPairs;
		vector<double> newPreComputedCors;
		vector<pair<RDLStateType, RDLStateType>> newRdlMaskAndDiffPairs;

		int pairsNum = diffAndMaskPairs.size();
		for (int i = 0; i < pairsNum; i++)
		{
			if (abs(preComputedCors[i]) >= bound)
			{
				newDiffAndMaskPairs.push_back(diffAndMaskPairs[i]);
				newPreComputedCors.push_back(preComputedCors[i]);
				newRdlMaskAndDiffPairs.push_back(rdlMaskAndDiffPairs[i]);
			}
		}

		diffAndMaskPairs = newDiffAndMaskPairs;
		preComputedCors = newPreComputedCors;
		rdlMaskAndDiffPairs = newRdlMaskAndDiffPairs;

		cout << __func__ << ": after filtering, pairs num " << diffAndMaskPairs.size() << endl;
	}


	double ComputeMultiRoundsRDLCor(const StateType& inputStateDiff, const StateType& outputStateMask, int startr, int endr) const
	{
		RDLStateType RDLStateDiff = RDLStateTypeFromStateType(inputStateDiff);
		for (int i = startr; i < endr; i++)
		{
			RDLDiffPropagateThroughNonLinearLayer(RDLStateDiff);
			RDLDiffPropagateThroughLinearLayer(RDLStateDiff, rotationalDifferencePerRoundConstant[i + 1]);

			// OutputState(cout, RDLStateDiff);
		}

		RDLStateType RDLStateMask = RDLStateTypeFromStateType(outputStateMask);
		for (int i = endr; i > startr; i--)
		{
			RDLLinPropagateThroughInvLinearLayer(RDLStateMask);
			RDLLinPropagateThroughInvNonLinearLayer(RDLStateMask);
			// OutputState(cout, RDLStateMask);
		}
		RDLLinPropagateThroughInvLinearLayer(RDLStateMask);

		// OutputState(cout, RDLStateMask);

		double initialForwardEstimation = RDLForwardSingleStateDiffThroughSboxEstimation(inputStateDiff, RDLStateMask);
		double initialBackwardEstimation = RDLBackwardSingleStateMaskThroughInvSboxEstimation(RDLStateDiff, outputStateMask);

		cout << __func__ << ": initial forward estimation " << initialForwardEstimation << endl;
		cout << __func__ << ": initial backward estimation " << initialBackwardEstimation << endl;

		if (initialForwardEstimation == -1 || initialBackwardEstimation == -1)
			return 0.0;

		vector<pair<StateType, StateType>> diffAndMaskPairs = { {inputStateDiff, outputStateMask} };
		vector<pair<RDLStateType, RDLStateType>> rdlMaskAndDiffPairs = { {RDLStateMask, RDLStateDiff} };
		vector<double> preComputedCors = { 1.0 };


		// this is the expansion process
		// here, we call backward divisions twice, so that the RDL correlation computation of 3-round EM is reduced to 1-round RDL correlation computation; if you found the program stuck at here, try to change the division methods to forward divisions; for example, you can also call RDLBackwardExpansionByOneRound(diffAndMaskPairs, preComputedCors, rdlMaskAndDiffPairs, startr, endr, false), and then call RDLForwardExpansionByOueRound(diffAndMaskPairs, preComputedCors, rdlMaskAndDiffPairs, startr, endr, false). 
		// If you would like to compute the RDL correlation of 4-round EM, then you should adjust here to call forward or backward divisions three times to reduce the computation into multiple 1-round RDL correlation computations
		RDLBackwardExpansionByOneRound(diffAndMaskPairs, preComputedCors, rdlMaskAndDiffPairs, startr, endr, false);
		RDLBackwardExpansionByOneRound(diffAndMaskPairs, preComputedCors, rdlMaskAndDiffPairs, startr, endr, false);


		// end of the expansion process

		int pairsNum = diffAndMaskPairs.size();
		cout << __func__ << ": pairs num " << pairsNum << endl;

		double finalCor = 0.0;

		int threadsNum = omp_get_max_threads();
#pragma omp parallel for num_threads(threadsNum)
		for (int i = 0; i < pairsNum; i++)
		{
			const pair<StateType, StateType>& diffAndMaskPair = diffAndMaskPairs[i];
			const StateType& diff = diffAndMaskPair.first;
			const StateType& mask = diffAndMaskPair.second;
			const RDLStateType& rdlMask = rdlMaskAndDiffPairs[i].first;
			const RDLStateType& rdlDiff = rdlMaskAndDiffPairs[i].second;
			const double& preComputedCor = preComputedCors[i];


			if (preComputedCor == 0.0)
				continue;

			// asconBase.OutputState(cout, diff);

			double forwardEstimation = RDLForwardSingleStateDiffThroughSboxEstimation(diff, rdlMask);
			double backwardEstimation = RDLBackwardSingleStateMaskThroughInvSboxEstimation(rdlDiff, mask);
			// cout << __func__ << ": pair index " << i << endl;
			// cout << __func__ << ": forward estimation " << forwardEstimation << endl;
			// cout << __func__ << ": backward estimation " << backwardEstimation << endl;
			// cout << __func__ << ": precomputed cor " << preComputedCor << endl;



			if (forwardEstimation == -1 || backwardEstimation == -1)
				continue;


			int diffActiveColumnsNum = rdBase.GetActiveColumnsNum(diff);
			int linActiveColumnsNum = linearBase.GetActiveColumnsNum(mask);
			// cout << "Difference active columns num " << diffActiveColumnsNum << endl;
			// cout << "Linear mask active columns num " << linActiveColumnsNum << endl;

			double theoreticalCor = ComputeMultipleRoundsRDLCorAutomatically(diff, mask, startr, endr, (startr + endr + 1) / 2);
			// double experimentalCor = ExperimentalComputation(diff, mask, startr, endr);
			// cout << "theoretical cor " << theoreticalCor << endl;
			// cout << "experimental cor " << experimentalCor << endl;

			/*
			{
			   if (abs(theoreticalCor - experimentalCor) / abs(experimentalCor) > 0.25)
			   {
				   cout << "Warning: theoretical cor and experimental cor differ a lot!" << endl;
				   cout << "Start Nr " << startr << endl;
				   cout << "End Nr " << endr << endl;
				   cout << "Input diff: " << endl;
				   asconBase.OutputState(cout, diff);
				   asconBase.OutputStateAs3DArr(cout, diff);
				   cout << "Output mask: " << endl;
				   asconBase.OutputState(cout, mask);
				   asconBase.OutputStateAs3DArr(cout, mask);
			   }
			}
			*/


			finalCor += preComputedCor * theoreticalCor;

			cout << "Update final Cor to " << finalCor << endl;
		}

		return finalCor;

	}

	// =====================================================================================
	// ================= Functions below are for the milp models of RDL propagation through xoodoo ==================
	// =========

	void MilpPropagateThroughSbox(GRBModel& model, GRBVar& v0, GRBVar& v1, GRBVar& v2, GRBVar& v3, GRBVar & v4) const
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

	void MilpPropagateThroughInvSbox(GRBModel& model, GRBVar& u0, GRBVar& u1, GRBVar& u2, GRBVar& u3, GRBVar& u4) const
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

	double ComputeMultipleRoundsRDLCorAutomatically(const StateType& inputStateDiff, const StateType& outputStateMask, int startr, int endr, int midr) const
	{
		if (endr == startr)
		{
			double RDLCor = 1.0;

			StateType sboxesInputDiff = inputStateDiff;
			StateType sboxesOutputMask = outputStateMask;
			linearBase.PropagateThroughInvLinearLayerWithCor2(sboxesOutputMask, RDLCor, rotationalDifferencePerRoundConstant[startr + 1]);

			for (int col = 0; col < 64; col++)
			{
				uint8_t inputColumnDiff = asconBase.GetColumn(sboxesInputDiff, col);
				uint8_t outputColumnMask = asconBase.GetColumn(sboxesOutputMask, col);
				RDLCor *= dlctCor[inputColumnDiff][outputColumnMask];
			}

			return RDLCor;
		}



		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_PoolSolutions, 200000000);
		env.set(GRB_IntParam_Threads, milpThreads);

		GRBModel  model = GRBModel(env);

		RDStateMILPType stateDiffVars = rdBase.MilpAddStateVars(model);
		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
				model.addConstr(stateDiffVars[row][col] == asconBase.GetBit(inputStateDiff, row, col));

		vector<RDStateMILPType> stateDiffsVars = { stateDiffVars };

		for (int i = startr; i < midr; i++)
		{
			rdBase.MilpPropagateThroughNonLinearLayer(model, stateDiffVars);
			rdBase.MilpPropagateThroughLinearLayer(model, stateDiffVars, rotationalDifferencePerRoundConstant[i + 1]);
			stateDiffsVars.emplace_back(stateDiffVars);
		}

		RDStateMILPType objStateVars = stateDiffVars;

		MilpPropagateThroughSboxes(model, stateDiffVars);

		LINStateMILPType stateMaskVars = linearBase.MilpAddStateVars(model);
		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
				model.addConstr(stateMaskVars[row][col] == asconBase.GetBit(outputStateMask, row, col));

		vector<LINStateMILPType> stateMasksVars = { stateMaskVars };
		for (int i = endr; i >= midr + 1; i--)
		{
			linearBase.MilpPropagateThroughInvLinearLayer(model, stateMaskVars);
			linearBase.MilpPropagateThroughInvNonLinearLayer(model, stateMaskVars);
			stateMasksVars.insert(stateMasksVars.begin(), stateMaskVars);
		}

		linearBase.MilpPropagateThroughInvLinearLayer(model, stateMaskVars);

		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
				model.addConstr(stateMaskVars[row][col] == stateDiffVars[row][col]);

		GRBLinExpr objExp = 0;
		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
			{
				objExp += objStateVars[row][col];
			}

		model.setObjective(objExp, GRB_MAXIMIZE);

		model.optimize();

		int solCount = model.get(GRB_IntAttr_SolCount);

		cout << "Get " << solCount << " RDL trails." << endl;

		if (solCount == 0)
			return 0.0;
		else
		{
			double RDLCor = 0.0;
			for (int solIndex = 0; solIndex < solCount; solIndex++)
			{
				model.set(GRB_IntParam_SolutionNumber, solIndex);

				vector<StateType> stateDiffs;
				vector<StateType> stateMasks;

				for (auto& stateDiffVars : stateDiffsVars)
				{
					StateType stateDiff = rdBase.MilpReadStateVars(stateDiffVars);
					stateDiffs.push_back(stateDiff);
				}

				// asconBase.OutputState(cout, stateDiffs[1]);

				for (auto& stateMaskVars : stateMasksVars)
				{
					StateType stateMask = linearBase.MilpReadStateVars(stateMaskVars);
					stateMasks.push_back(stateMask);
				}

				double difProb = rdBase.ComputeDifferentialTrailProbability(stateDiffs, startr, midr - 1);
				double linCor2 = linearBase.ComputeLinearTrailCorrelationSquare(stateMasks, midr + 1, endr);
				double dlctCor = ComputeOneRoundRDLCor(stateDiffs.back(), stateMasks.front(), midr);

				// cout << difProb << "," << linCor2 << "," << dlctCor << endl;

				RDLCor += difProb * linCor2 * dlctCor;
			}

			return RDLCor;
		}
	}

	vector<pair<StateType, double>> RDLForwardSingleStateDiffThroughOneRound(const StateType& inputStateDiff, const RDLStateType& outputRDLMask, int startr) const
	{
		vector<pair<StateType, double>> outputStateDiffAndProbs;

		array<int, 64> outputUnknownLinearMaskBitsMask = { 0 };
		array<int, 64> outputFixedLinearMaskBits = { 0 };

		for (int col = 0; col < 64; col++)
			for (int row = 0; row < 5; row++)
				if (outputRDLMask[row][col] == -1)
					outputUnknownLinearMaskBitsMask[col] |= (1 << (4-row));
				else if (outputRDLMask[row][col] == 1)
					outputFixedLinearMaskBits[col] |= (1 << (4-row));

		vector<int> activeColumns;
		vector<int> activeColumnsOutputDiffIndexes;
		vector<int> activeColumnsInputDiffs;
		vector<vector<int>> activeColumnsPossibleOutputDiffs;
		vector<vector<double>> activeColumnsPossibleProbs;

		for (int col = 0; col < 64; col++)
		{
			int inputDiff = asconBase.GetColumn(inputStateDiff, col);

			const vector<int>& possibleOutputDiffs = RDLOutputDiffsPerInputDifferenceForSbox[inputDiff][outputUnknownLinearMaskBitsMask[col]][outputFixedLinearMaskBits[col]];

			if (possibleOutputDiffs.size() == 0)
				return outputStateDiffAndProbs;

			if (inputDiff != 0)
			{
				activeColumns.push_back(col);
				activeColumnsOutputDiffIndexes.push_back(0);
				activeColumnsInputDiffs.push_back(inputDiff);
				activeColumnsPossibleOutputDiffs.push_back(possibleOutputDiffs);

				vector<double> possibleProbs;
				for (auto& outputDiff : possibleOutputDiffs)
					possibleProbs.push_back(RDLProbPerInputDifferenceForSbox[inputDiff][outputUnknownLinearMaskBitsMask[col]][outputFixedLinearMaskBits[col] | outputDiff]);
				activeColumnsPossibleProbs.push_back(possibleProbs);
			}
		}

		int activeColumnsNum = activeColumns.size();
		bool isAllActiveColumnsProcessed = false;
		while (!isAllActiveColumnsProcessed)
		{
			StateType sboxesOutputDiff = { 0 };
			double prob = 1.0;

			for (int activeColumnIndex = 0; activeColumnIndex < activeColumnsNum; activeColumnIndex++)
			{
				const int& activeColumnPos = activeColumns[activeColumnIndex];
				const int& activeColumnOutputDiffIndex = activeColumnsOutputDiffIndexes[activeColumnIndex];
				const int& activeColumnInputDiff = activeColumnsInputDiffs[activeColumnIndex];
				const vector<int>& activeColumnPossibleOutputDiffs = activeColumnsPossibleOutputDiffs[activeColumnIndex];
				const vector<double>& activeColumnPossibleProbs = activeColumnsPossibleProbs[activeColumnIndex];
				const int& activeColumnCurOutputDiff = activeColumnPossibleOutputDiffs[activeColumnOutputDiffIndex];
				const double& activeColumnCurProb = activeColumnPossibleProbs[activeColumnOutputDiffIndex];

				asconBase.SetColumn(sboxesOutputDiff, activeColumnPos, activeColumnCurOutputDiff);

				prob *= activeColumnCurProb;
			}

			StateType outputStateDiff = sboxesOutputDiff;
			rdBase.PropagateThroughLinearLayer(outputStateDiff, rotationalDifferencePerRoundConstant[startr + 1]);

			outputStateDiffAndProbs.push_back({ outputStateDiff, prob });

			if (SwitchToNextPossibleOutputDifferencesOrInputLinearMasks(activeColumnsOutputDiffIndexes, activeColumnsPossibleOutputDiffs))
				isAllActiveColumnsProcessed = false;
			else
				isAllActiveColumnsProcessed = true;
		}

		return outputStateDiffAndProbs;
	}

	// midr is the number of rounds where you put the DLCT in
	vector<pair<StateType, double>> RDLForwardSingleStateDiffThroughOneRoundAutomatically(const StateType& inputStateDiff, const StateType& outputStateMask, int startr, int endr, int midr) const
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_IntParam_Threads, milpThreads);

		GRBModel  model = GRBModel(env);

		RDStateMILPType stateDiffVars = rdBase.MilpAddStateVars(model);
		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
				model.addConstr(stateDiffVars[row][col] == asconBase.GetBit(inputStateDiff, row, col));

		RDStateMILPType targetSolVars;
		for (int i = startr; i < midr; i++)
		{
			rdBase.MilpPropagateThroughNonLinearLayer(model, stateDiffVars);
			rdBase.MilpPropagateThroughLinearLayer(model, stateDiffVars, rotationalDifferencePerRoundConstant[i + 1]);

			if (i == startr)
				targetSolVars = stateDiffVars;
		}

		MilpPropagateThroughSboxes(model, stateDiffVars);

		LINStateMILPType stateMaskVars = linearBase.MilpAddStateVars(model);
		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
				model.addConstr(stateMaskVars[row][col] == asconBase.GetBit(outputStateMask, row, col));

		for (int i = endr; i >= midr + 1; i--)
		{
			linearBase.MilpPropagateThroughInvLinearLayer(model, stateMaskVars);
			linearBase.MilpPropagateThroughInvNonLinearLayer(model, stateMaskVars);
		}

		linearBase.MilpPropagateThroughInvLinearLayer(model, stateMaskVars);

		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
				model.addConstr(stateMaskVars[row][col] == stateDiffVars[row][col]);

		GRBLinExpr objExp = 0;
		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
			{
				objExp += targetSolVars[row][col];
			}

		set<StateType> targetSolutions;
		FullStateForwardExpandCallback cb = FullStateForwardExpandCallback(asconBase, targetSolVars, targetSolutions);
		model.setCallback(&cb);

		model.setObjective(objExp, GRB_MAXIMIZE);

		model.optimize();

		vector<pair<StateType, double>> outputStateDiffAndProbs;
		auto ddtProb = rdBase.GetDDTProb();

		for (auto& targetSolution : targetSolutions)
		{
			StateType outputStateDiff = targetSolution;
			double prob = rdBase.ComputeDifferentialProbabilityManually(inputStateDiff, outputStateDiff, startr, startr);

			outputStateDiffAndProbs.push_back({ outputStateDiff, prob });
		}

		return outputStateDiffAndProbs;
	}

	vector<pair<StateType, double>> RDLBackwardSinlgeStateMaskThroughOneRoundAutomatically(const StateType& inputStateDiff, const StateType& outputStateMask, int startr, int endr, int midr) const
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_IntParam_Threads, milpThreads);

		GRBModel  model = GRBModel(env);

		RDStateMILPType stateDiffVars = rdBase.MilpAddStateVars(model);
		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
				model.addConstr(stateDiffVars[row][col] == asconBase.GetBit(inputStateDiff, row, col));

		RDStateMILPType targetSolVars;
		for (int i = startr; i < midr; i++)
		{
			rdBase.MilpPropagateThroughNonLinearLayer(model, stateDiffVars);
			rdBase.MilpPropagateThroughLinearLayer(model, stateDiffVars, rotationalDifferencePerRoundConstant[i + 1]);

		}

		MilpPropagateThroughSboxes(model, stateDiffVars);

		LINStateMILPType stateMaskVars = linearBase.MilpAddStateVars(model);
		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
				model.addConstr(stateMaskVars[row][col] == asconBase.GetBit(outputStateMask, row, col));

		for (int i = endr; i >= midr + 1; i--)
		{
			linearBase.MilpPropagateThroughInvLinearLayer(model, stateMaskVars);
			linearBase.MilpPropagateThroughInvNonLinearLayer(model, stateMaskVars);

			if (i == endr)
				targetSolVars = stateMaskVars;
		}

		linearBase.MilpPropagateThroughInvLinearLayer(model, stateMaskVars);

		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
				model.addConstr(stateMaskVars[row][col] == stateDiffVars[row][col]);

		GRBLinExpr objExp = 0;
		for (int row = 0; row < 5; row++)
			for (int col = 0; col < 64; col++)
			{
				objExp += targetSolVars[row][col];
			}

		set<StateType> targetSolutions;
		FullStateBackwardExpandCallback cb = FullStateBackwardExpandCallback(asconBase, targetSolVars, targetSolutions);
		model.setCallback(&cb);

		model.setObjective(objExp, GRB_MAXIMIZE);

		model.optimize();

		vector<pair<StateType, double>> inputStateMaskAndCor2s;
		auto latCor2 = linearBase.GetLATCor2();
		for (auto& targetSolution : targetSolutions)
		{
			StateType inputStateMask = targetSolution;
			double cor2 = linearBase.ComputeApproximationCorrelationSquareManually(inputStateMask, outputStateMask, endr, endr);
			inputStateMaskAndCor2s.push_back({ inputStateMask, cor2 });
		}

		return inputStateMaskAndCor2s;
	}

	vector<pair<StateType, double>> RDLBackwardSingleStateMaskThroughOneRound(const RDLStateType& inputRDLDiff, const StateType& outputMask) const
	{
		StateType sboxesOutputMask = outputMask;
		linearBase.PropagateThroughInvLinearLayer(sboxesOutputMask);

		vector<pair<StateType, double>> inputStateRDLMaskAndCors;

		array<int, 64> inputUnknownDiffBitsMask = { 0 };
		array<int, 64> inputFixedDiffBits = { 0 };
		for (int col = 0; col < 64; col++)
		{
			for (int row = 0; row < 5; row++)
				if (inputRDLDiff[row][col] == -1)
					inputUnknownDiffBitsMask[col] |= (1 << (4-row));
				else if (inputRDLDiff[row][col] == 1)
					inputFixedDiffBits[col] |= (1 << (4-row));
		}

		vector<int> activeColumns;
		vector<int> activeColumnsInputMaskIndexes;
		vector<int> activeColumnsOutputMasks;
		vector<vector<int>> activeColumnsPossibleInputMasks;
		vector<vector<double>> activeColumnsPossibleCor2s;

		for (int col = 0; col < 64; col++)
		{
			int outputMask = asconBase.GetColumn(sboxesOutputMask, col);
			const vector<int>& possibleInputMasks = RDLInputMasksPerOutputMaskForSbox[outputMask][inputUnknownDiffBitsMask[col]][inputFixedDiffBits[col]];

			if (possibleInputMasks.size() == 0)
				return inputStateRDLMaskAndCors;

			if (outputMask != 0)
			{
				activeColumns.push_back(col);
				activeColumnsInputMaskIndexes.push_back(0);
				activeColumnsOutputMasks.push_back(outputMask);
				activeColumnsPossibleInputMasks.push_back(possibleInputMasks);
				vector<double> possibleCor2s;
				for (auto& inputMask : possibleInputMasks)
					possibleCor2s.push_back(RDLCor2PerOutputLinearMaskForSbox[outputMask][inputUnknownDiffBitsMask[col]][inputFixedDiffBits[col] | inputMask]);
				activeColumnsPossibleCor2s.push_back(possibleCor2s);
			}
		}

		int activeColumnsNum = activeColumns.size();

		bool isAllActiveColumnsProcessed = false;
		while (!isAllActiveColumnsProcessed)
		{
			StateType sboxesInputMask = { 0 };
			double cor2 = 1.0;
			for (int activeColumnIndex = 0; activeColumnIndex < activeColumnsNum; activeColumnIndex++)
			{
				const int& activeColumnPos = activeColumns[activeColumnIndex];
				const int& activeColumnInputMaskIndex = activeColumnsInputMaskIndexes[activeColumnIndex];
				const int& activeColumnOutputMask = activeColumnsOutputMasks[activeColumnIndex];
				const vector<int>& activeColumnPossibleInputMasks = activeColumnsPossibleInputMasks[activeColumnIndex];
				const vector<double>& activeColumnPossibleCor2s = activeColumnsPossibleCor2s[activeColumnIndex];
				const int& activeColumnCurInputMask = activeColumnPossibleInputMasks[activeColumnInputMaskIndex];
				const double& activeColumnCurCor2 = activeColumnPossibleCor2s[activeColumnInputMaskIndex];

				asconBase.SetColumn(sboxesInputMask, activeColumnPos, activeColumnCurInputMask);
				cor2 *= activeColumnCurCor2;
			}


			inputStateRDLMaskAndCors.push_back({ sboxesInputMask, cor2 });

			if (SwitchToNextPossibleOutputDifferencesOrInputLinearMasks(activeColumnsInputMaskIndexes, activeColumnsPossibleInputMasks))
				isAllActiveColumnsProcessed = false;
			else
				isAllActiveColumnsProcessed = true;
		}

		return inputStateRDLMaskAndCors;

	}

	void OutputState(ostream& os, const RDLStateType& state) const
	{
		for (int row = 0; row < 5; row++)
		{
			os << "Row " << row << ": ";
			for(int col = 0; col < 64; col++)
				os << int(state[row][col]) << " ";
			os << endl;
		}
	}



};


#pragma once
#include"differential-linear.h"
#include"internalDifferential.h"
#include"internalLinear.h"


using namespace std;

using IDLLaneType = array<int, 16>;
using IDLSheetType = array<IDLLaneType, 3>;
using IDLPlaneType = array<IDLLaneType, 4>;
using IDLStateTypeBySheet = array<IDLSheetType, 4>;

IDLStateTypeBySheet IDLStateTypeBySheetFromStateTypeBySheet(const HalfStateTypeBySheet& state)
{
	IDLStateTypeBySheet rdlState = { 0 };
	for (int x = 0; x < 4; x++)
	{
		for (int y = 0; y < 3; y++)
		{
			for (int z = 0; z < 16; z++)
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

IDLLaneType RotateLeft(const IDLLaneType lane, int shift)
{
	IDLLaneType newLane;
	for (int i = 0; i < 16; i++)
		newLane[i] = lane[(i - shift + 16) % 16];
	return newLane;
}

IDLLaneType RotateRight(const IDLLaneType lane, int shift)
{
	return RotateLeft(lane, 16 - shift);
}

class InternalDifferentialLinearBase
{
private:
	const XoodooBase& xoodooBase;
	const RotationalDifferentialLinearBase& rdlBase;
	const InternalDifferentialBase& idBase;
	const InternalLinearBase& iLinearBase;

	int milpThreads;

	array<array<double, 8>, 8> dlctCor;
	array<HalfLaneType, 12> internalDifferencePerRoundConstants;

	array<array<array<double, 8>, 8>, 8> IDLCor2PerOutputLinearMaskForChi;
	array<array<array<vector<int>, 8>, 8>, 8> IDLInputMasksPerOutputMaskForChi;
	array<array<array<double, 8>, 8>, 8> IDLProbPerInputDifferenceForChi;

	array<array<array<vector<int>, 8>, 8>, 8> IDLOutputDiffsPerInputDifferenceForChi;

	vector<array<int, 7>> milpChiInequalities;

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
	InternalDifferentialLinearBase(const XoodooBase& xoodooBase, const InternalLinearBase& iLinearBase, const InternalDifferentialBase& idBase, const RotationalDifferentialLinearBase& rdlBase, int milpThreads = 2) : xoodooBase(xoodooBase), iLinearBase(iLinearBase), idBase(idBase), rdlBase(rdlBase), milpThreads(milpThreads)
	{
		internalDifferencePerRoundConstants = idBase.getInternalDifferencePerRoundConstants();
		dlctCor = rdlBase.getDLCTCor();
		IDLCor2PerOutputLinearMaskForChi = rdlBase.getRDLCor2PerOutputLinearMaskForChi();
		IDLInputMasksPerOutputMaskForChi = rdlBase.getRDLInputMasksPerOutputMaskForChi();
		IDLProbPerInputDifferenceForChi = rdlBase.getRDLProbPerInputDifferenceForChi();
		IDLOutputDiffsPerInputDifferenceForChi = rdlBase.getRDLOutputDiffsPerInputDifferenceForChi();

		milpChiInequalities = { {1, -1, 0, 0, 1, -1, 0}, {1, 0, -1, 0, -1, 1, 0}, {1, 0, 0, -1, -1, 0, 1}, {1, -1, 0, 1, 0, 0, -1}, {1, 0, 1, -1, 0, -1, 0}, {1, 0, -1, 0, 0, 1, -1}, {1, 1, -1, 0, -1, 0, 0}, {1, -1, 0, 0, 1, 0, -1}, {1, -1, 1, 0, 0, -1, 0}, {1, 0, -1, 1, 0, 0, -1}, {1, 0, 0, -1, 0, -1, 1}, {1, 1, 0, -1, -1, 0, 0} };
	}

	void OutputIDLState(ostream& os, const RDLStateTypeBySheet& state) const
	{
		for (int y = 0; y < 3; y++)
		{
			for (int x = 0; x < 4; x++)
			{
				for (int z = 0; z < 16; z++)
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
	void IDLDiffPropagateThroughRhoWest(IDLStateTypeBySheet& state) const
	{
		IDLStateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = tmp[(x + 3) % 4][1];

		for (int x = 0; x < 4; x++)
			state[x][2] = RotateLeft(tmp[x][2], 11);
	}

	void IDLLinPropagateThroughRevRhoWest(IDLStateTypeBySheet& state) const
	{
		IDLStateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = tmp[(x + 1) % 4][1];

		for (int x = 0; x < 4; x++)
			state[x][2] = RotateRight(tmp[x][2], 11);
	}

	void IDLDiffPropagateThroughRhoEast(IDLStateTypeBySheet& state) const
	{
		IDLStateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = RotateLeft(tmp[x][1], 1);

		for (int x = 0; x < 4; x++)
			state[x][2] = RotateLeft(tmp[(x + 2) % 4][2], 8);
	}

	void IDLLinPropagateThroughRevRhoEast(IDLStateTypeBySheet& state) const
	{
		IDLStateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = RotateRight(tmp[x][1], 1);

		for (int x = 0; x < 4; x++)
			state[x][2] = RotateRight(tmp[(x + 2) % 4][2], 8);
	}

	void IDLDiffPropagateThroughLota(IDLStateTypeBySheet& state, const LaneType rcDiff) const
	{
		for (int z = 0; z < 16; z++)
			if (state[0][0][z] != -1)
				state[0][0][z] = xor2(state[0][0][z], (rcDiff >> z) & 1);
	}

	void IDLDiffPropagateThroughTheta(IDLStateTypeBySheet& state) const
	{
		IDLPlaneType parity;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
				parity[x][z] = xor3(state[x][0][z], state[x][1][z], state[x][2][z]);

		for (int x = 0; x < 4; x++)
		{
			IDLLaneType one = RotateLeft(parity[(x + 3) % 4], 5);
			IDLLaneType two = RotateLeft(parity[(x + 3) % 4], 14);
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					state[x][y][z] = xor3(state[x][y][z], one[z], two[z]);

		}
	}

	void IDLLinPropagateThroughRevTheta(IDLStateTypeBySheet& state) const
	{
		IDLPlaneType parity;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
				parity[x][z] = xor3(state[x][0][z], state[x][1][z], state[x][2][z]);

		for (int x = 0; x < 4; x++)
		{
			IDLLaneType one = RotateRight(parity[(x + 1) % 4], 5);
			IDLLaneType two = RotateRight(parity[(x + 1) % 4], 14);
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					state[x][y][z] = xor3(state[x][y][z], one[z], two[z]);

		}
	}

	void IDLDiffPropagateThroughChiColumn(int& y0, int& y1, int& y2) const
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

	void IDLLinPropagateThroughRevChiColumn(int& y0, int& y1, int& y2) const
	{
		IDLDiffPropagateThroughChiColumn(y0, y1, y2);
	}

	void IDLDiffPropagateThroughChi(IDLStateTypeBySheet& state) const
	{
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
				IDLDiffPropagateThroughChiColumn(state[x][0][z], state[x][1][z], state[x][2][z]);
	}

	void IDLLinPropagateThroughRevChi(IDLStateTypeBySheet& state) const
	{
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
				IDLLinPropagateThroughRevChiColumn(state[x][0][z], state[x][1][z], state[x][2][z]);
	}

	void IDLDiffPropagateThroughLinearLayer(IDLStateTypeBySheet& state, const HalfLaneType rcDiff) const
	{
		IDLDiffPropagateThroughRhoEast(state);
		IDLDiffPropagateThroughTheta(state);
		IDLDiffPropagateThroughRhoWest(state);
		IDLDiffPropagateThroughLota(state, rcDiff);
	}

	void IDLLinPropagateThroughRevLinearLayer(IDLStateTypeBySheet& state) const
	{
		IDLLinPropagateThroughRevRhoWest(state);
		IDLLinPropagateThroughRevTheta(state);
		IDLLinPropagateThroughRevRhoEast(state);
	}

	void IDLDiffPropagateThroughNonLinearLayer(IDLStateTypeBySheet& state) const
	{
		IDLDiffPropagateThroughChi(state);
	}

	void IDLLinPropagateThroughRevNonLinearLayer(IDLStateTypeBySheet& state) const
	{
		IDLLinPropagateThroughRevChi(state);
	}

	// =====================================================================================
	// ================= Functions below are for computing the correlation of RDL distinguisher in the middle part ==================
	// =====================================================================================

	int IDLFowardSingleStateDiffThrouChiEstimation(const HalfStateTypeBySheet& inputStateDiff, const IDLStateTypeBySheet& outputIDLLinMask) const
	{
		array<array<int, 16>, 4> outputUnknownLinearMaskBitsMask = { 0 };
		array<array<int, 16>, 4> outputFixedLinearMaskBits = { 0 };
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
			{
				for (int y = 0; y < 3; y++)
					if (outputIDLLinMask[x][y][z] == -1)
						outputUnknownLinearMaskBitsMask[x][z] |= (1 << y);
					else if (outputIDLLinMask[x][y][z] == 1)
						outputFixedLinearMaskBits[x][z] |= (1 << y);
			}



		double estimatedNumLog2 = 0;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
			{
				int inputDiff = xoodooBase.GetColumn(inputStateDiff, x, z);

				const vector<int>& possibleOutputDiffs = IDLOutputDiffsPerInputDifferenceForChi[inputDiff][outputUnknownLinearMaskBitsMask[x][z]][outputFixedLinearMaskBits[x][z]];

				if (possibleOutputDiffs.size() == 0)
				{
					return -1;
				}

				if (inputDiff != 0)
				{
					estimatedNumLog2 += log2(possibleOutputDiffs.size());
				}


			}

		return estimatedNumLog2;
	}

	int IDLBackwardSingleStateLinearMaskThroughRevChiEstimation(const IDLStateTypeBySheet& inputStateIDLDiff, const HalfStateTypeBySheet& outputLinMask) const
	{
		array<array<int, 16>, 4> inputUnknownDiffBitsMask = { 0 };
		array<array<int, 16>, 4> inputFixedDiffBits = { 0 };
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
				for (int y = 0; y < 3; y++)
					if (inputStateIDLDiff[x][y][z] == -1)
						inputUnknownDiffBitsMask[x][z] |= (1 << y);
					else if (inputStateIDLDiff[x][y][z] == 1)
						inputFixedDiffBits[x][z] |= (1 << y);




		double estimatedNumLog2 = 0;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
			{
				int outputLinearMask = xoodooBase.GetColumn(outputLinMask, x, z);
				const vector<int>& possibleInputLinearMasks = IDLInputMasksPerOutputMaskForChi[outputLinearMask][inputUnknownDiffBitsMask[x][z]][inputFixedDiffBits[x][z]];

				if (possibleInputLinearMasks.size() == 0)
				{
					return -1;
				}

				if (outputLinearMask != 0)
				{
					estimatedNumLog2 += log2(possibleInputLinearMasks.size());

				}
			}

		return estimatedNumLog2;
	}


	vector<pair<HalfStateTypeBySheet, double>> IDLForwardSingleStateDiffThroughChi(const HalfStateTypeBySheet& inputStateDiff, const IDLStateTypeBySheet& outputIDLLinMask) const
	{
		vector<pair<HalfStateTypeBySheet, double>> outputStateDiffAndProbs;

		array<array<int, 16>, 4> outputUnknownLinearMaskBitsMask = { 0 };
		array<array<int, 16>, 4> outputFixedLinearMaskBits = { 0 };
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
			{
				for (int y = 0; y < 3; y++)
					if (outputIDLLinMask[x][y][z] == -1)
						outputUnknownLinearMaskBitsMask[x][z] |= (1 << y);
					else if (outputIDLLinMask[x][y][z] == 1)
						outputFixedLinearMaskBits[x][z] |= (1 << y);
			}


		vector<pair<int, int>> activeColumns;
		vector<int> activeColumnsOutputDiffIndexes;
		vector<int> activeColumnsInputDiffs;
		vector<vector<int>> activeColumnsPossibleOutputDiffs;
		vector<vector<double>> activeColumnsPossibleProbs;

		double estimatedNumLog2 = 0;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
			{
				int inputDiff = xoodooBase.GetColumn(inputStateDiff, x, z);

				const vector<int>& possibleOutputDiffs = IDLOutputDiffsPerInputDifferenceForChi[inputDiff][outputUnknownLinearMaskBitsMask[x][z]][outputFixedLinearMaskBits[x][z]];

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
						possibleProbs.push_back(IDLProbPerInputDifferenceForChi[inputDiff][outputUnknownLinearMaskBitsMask[x][z]][outputFixedLinearMaskBits[x][z] | outputDiff]);
					activeColumnsPossibleProbs.push_back(possibleProbs);

				}


			}

		// cout << "Estimated results num log2: " << estimatedNumLog2 << endl;






		int activeColumnsNum = activeColumns.size();
		// cout << "Active columns num: " << activeColumnsNum << endl;

		bool isAllActiveColumnsProcessed = false;
		while (!isAllActiveColumnsProcessed)
		{
			HalfStateTypeBySheet outputStateDiff = { 0 };
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

	vector<pair<HalfStateTypeBySheet, double>> IDLBackwardSingleStateLinearMaskThroughRevChi(const IDLStateTypeBySheet& inputStateIDLDiff, const HalfStateTypeBySheet& outputLinMask) const
	{
		// OutputRDLState(cout, inputStateRDLDiff);
		// xoodooBase.OutputState(cout, outputLinMask);

		vector<pair<HalfStateTypeBySheet, double>> inputStateLinearMaskAndCor2s;

		array<array<int, 16>, 4> inputUnknownDiffBitsMask = { 0 };
		array<array<int, 16>, 4> inputFixedDiffBits = { 0 };
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
				for (int y = 0; y < 3; y++)
					if (inputStateIDLDiff[x][y][z] == -1)
						inputUnknownDiffBitsMask[x][z] |= (1 << y);
					else if (inputStateIDLDiff[x][y][z] == 1)
						inputFixedDiffBits[x][z] |= (1 << y);



		vector<pair<int, int>> activeColumns;
		vector<int>	activeColumnsInputLinearMaskIndexes;
		vector<int> activeColumnsOutputLinearMasks;
		vector<vector<int>> activeColumnsPossibleInputLinearMasks;
		vector<vector<double>> activeColumnsPossibleCor2s;

		double estimatedNumLog2 = 0;
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
			{
				int outputLinearMask = xoodooBase.GetColumn(outputLinMask, x, z);
				const vector<int>& possibleInputLinearMasks = IDLInputMasksPerOutputMaskForChi[outputLinearMask][inputUnknownDiffBitsMask[x][z]][inputFixedDiffBits[x][z]];

				if (possibleInputLinearMasks.size() == 0)
					return inputStateLinearMaskAndCor2s;

				if (outputLinearMask != 0)
				{
					estimatedNumLog2 += log2(possibleInputLinearMasks.size());

					activeColumns.push_back({ x, z });
					activeColumnsInputLinearMaskIndexes.push_back(0);
					activeColumnsOutputLinearMasks.push_back(outputLinearMask);
					activeColumnsPossibleInputLinearMasks.push_back(possibleInputLinearMasks);
					vector<double> possibleCor2s;
					for (auto& inputLinearMask : possibleInputLinearMasks)
						possibleCor2s.push_back(IDLCor2PerOutputLinearMaskForChi[outputLinearMask][inputUnknownDiffBitsMask[x][z]][inputFixedDiffBits[x][z] | inputLinearMask]);
					activeColumnsPossibleCor2s.push_back(possibleCor2s);

				}
			}


		int activeColumnsNum = activeColumns.size();


		bool isAllActiveColumnsProcessed = false;
		while (!isAllActiveColumnsProcessed)
		{
			HalfStateTypeBySheet inputStateLinearMask = { 0 };
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

			inputStateLinearMaskAndCor2s.push_back({ inputStateLinearMask, cor2 });

			if (SwitchToNextPossibleOutputDifferencesOrInputLinearMasks(activeColumnsInputLinearMaskIndexes, activeColumnsPossibleInputLinearMasks))
				isAllActiveColumnsProcessed = false;
			else
				isAllActiveColumnsProcessed = true;
		}

		return inputStateLinearMaskAndCor2s;
	}

	void IDLForwardExpansionByOneRound(vector<pair<HalfStateTypeBySheet, HalfStateTypeBySheet>>& DiffAndLinearMaskPairs, vector<double>& preComputedCors, vector<pair<IDLStateTypeBySheet, IDLStateTypeBySheet>>& IDLLinearMaskAndDiffPairs, int& startr, int& endr, int totalRounds, bool isCallback) const
	{
		map<pair<HalfStateTypeBySheet, HalfStateTypeBySheet>, double> m0;
		map<pair<HalfStateTypeBySheet, HalfStateTypeBySheet>, pair<IDLStateTypeBySheet, IDLStateTypeBySheet>> m1;

		int rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + 1);
		HalfLaneType rcDiff = internalDifferencePerRoundConstants[rcDiffIndex];

		int pairsNum = DiffAndLinearMaskPairs.size();
		cout << __func__ << ": to expand " << pairsNum << " pairs" << endl;
		cout << "Start round: " << startr << endl;
		cout << "End round: " << endr << endl;

		int threadsNum = omp_get_max_threads();
#pragma omp parallel for num_threads(threadsNum)
		for (int i = 0; i < pairsNum; i++)
		{
			// cout << __func__ << ": pair index " << i << endl;
			const pair<HalfStateTypeBySheet, HalfStateTypeBySheet> diffAndLinearMaskPair = DiffAndLinearMaskPairs[i];
			const HalfStateTypeBySheet& diff = diffAndLinearMaskPair.first;
			const HalfStateTypeBySheet& linearMask = diffAndLinearMaskPair.second;
			const double& preComputedCor = preComputedCors[i];
			const IDLStateTypeBySheet& IDLLinearMask = IDLLinearMaskAndDiffPairs[i].first;
			const IDLStateTypeBySheet& IDLDiff = IDLLinearMaskAndDiffPairs[i].second;


			IDLStateTypeBySheet newIDLLinearMask = IDLStateTypeBySheetFromStateTypeBySheet(linearMask);
			for (int j = endr; j > startr + 1; j--)
			{
				IDLLinPropagateThroughRevNonLinearLayer(newIDLLinearMask);
				IDLLinPropagateThroughRevLinearLayer(newIDLLinearMask);
			}

			int forwardEstimation = IDLFowardSingleStateDiffThrouChiEstimation(diff, IDLLinearMask);
			if (forwardEstimation == -1)
				continue;
			else
			{
				; // cout << "pair index " << i << " estimation " << forwardEstimation << endl;
			}

			vector<pair<HalfStateTypeBySheet, double>> newDiffAndProbs;
			if (!isCallback)
				newDiffAndProbs = IDLForwardSingleStateDiffThroughChi(diff, IDLLinearMask);
			else
				newDiffAndProbs = IDLForwardSingleStateDiffThroughChiAutomatically(diff, linearMask, startr, endr, totalRounds);

#pragma omp critical 
			{
				for (auto& diffAndProb : newDiffAndProbs)
				{
					HalfStateTypeBySheet stateDiff = diffAndProb.first;
					double prob = diffAndProb.second;
					idBase.PropagateThroughLinearLayer(stateDiff, rcDiff);

					if (m0.find({ stateDiff, linearMask }) != m0.end())
						m0[{stateDiff, linearMask}] += prob * preComputedCor;
					else
					{
						m0[{stateDiff, linearMask}] = prob * preComputedCor;

						IDLStateTypeBySheet newIDLDiff = IDLStateTypeBySheetFromStateTypeBySheet(stateDiff);
						for (int j = startr + 1; j < endr; j++)
						{
							int tmpRcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, j + 1);
							HalfLaneType tmpRcDiff = internalDifferencePerRoundConstants[tmpRcDiffIndex];
							IDLDiffPropagateThroughNonLinearLayer(newIDLDiff);
							IDLDiffPropagateThroughLinearLayer(newIDLDiff, tmpRcDiff);
						}

						m1[{stateDiff, linearMask}] = { newIDLLinearMask, newIDLDiff };
					}
				}
			}
		}

		DiffAndLinearMaskPairs.clear();
		preComputedCors.clear();
		IDLLinearMaskAndDiffPairs.clear();

		for (auto& item : m0)
		{
			DiffAndLinearMaskPairs.push_back(item.first);
			preComputedCors.push_back(item.second);
			IDLLinearMaskAndDiffPairs.push_back(m1[item.first]);
		}

		startr += 1;
	}


	void IDLBackwardExpansionByOneRound(vector<pair<HalfStateTypeBySheet, HalfStateTypeBySheet>>& DiffAndLinearMaskPairs, vector<double>& preComputedCors, vector<pair<IDLStateTypeBySheet, IDLStateTypeBySheet>>& IDLLinearMaskAndDiffPairs, int& startr, int& endr, int totalRounds, bool isCallback) const
	{
		map<pair<HalfStateTypeBySheet, HalfStateTypeBySheet>, double> m0;
		map<pair<HalfStateTypeBySheet, HalfStateTypeBySheet>, pair<IDLStateTypeBySheet, IDLStateTypeBySheet>> m1;

		int rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, endr);
		HalfLaneType rcDiff = internalDifferencePerRoundConstants[rcDiffIndex];

		int pairsNum = DiffAndLinearMaskPairs.size();
		cout << __func__ << ": to expand " << pairsNum << " pairs" << endl;
		cout << "Start round: " << startr << endl;
		cout << "End round: " << endr << endl;

		int threadsNum = omp_get_max_threads();
#pragma omp parallel for num_threads(threadsNum)
		for (int i = 0; i < pairsNum; i++)
		{
			const pair<HalfStateTypeBySheet, HalfStateTypeBySheet> diffAndLinearMaskPair = DiffAndLinearMaskPairs[i];
			const HalfStateTypeBySheet& diff = diffAndLinearMaskPair.first;
			const HalfStateTypeBySheet& linearMask = diffAndLinearMaskPair.second;
			const double& preComputedCor = preComputedCors[i];
			const IDLStateTypeBySheet& IDLLinearMask = IDLLinearMaskAndDiffPairs[i].first;
			const IDLStateTypeBySheet& IDLDiff = IDLLinearMaskAndDiffPairs[i].second;

			IDLStateTypeBySheet newIDLDiff = IDLStateTypeBySheetFromStateTypeBySheet(diff);
			for (int j = startr; j < endr - 1; j++)
			{
				int tmpRcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, j + 1);
				HalfLaneType tmpRcDiff = internalDifferencePerRoundConstants[tmpRcDiffIndex];
				IDLDiffPropagateThroughNonLinearLayer(newIDLDiff);
				IDLDiffPropagateThroughLinearLayer(newIDLDiff, tmpRcDiff);
			}

			int backwardEstimation = IDLBackwardSingleStateLinearMaskThroughRevChiEstimation(IDLDiff, linearMask);
			if (backwardEstimation == -1)
				continue;
			else
			{
				cout << "pair index " << i << " estimation " << backwardEstimation << endl;
			}

			vector<pair<HalfStateTypeBySheet, double>> newLinearMaskAndCor2s;
			if (!isCallback)
				newLinearMaskAndCor2s = IDLBackwardSingleStateLinearMaskThroughRevChi(IDLDiff, linearMask);
			else
				newLinearMaskAndCor2s = IDLBackwardSingleStateLinearMaskThroughRevChiAutomatically(diff, linearMask, startr, endr, totalRounds);

#pragma omp critical 
			{
				for (auto& linearMaskAndCor2 : newLinearMaskAndCor2s)
				{
					HalfStateTypeBySheet stateLinearMask = linearMaskAndCor2.first;
					double cor2 = linearMaskAndCor2.second;
					iLinearBase.PropagateThroughRevLinearLayerWithCor2(stateLinearMask, cor2, rcDiff);

					if (m0.find({ diff, stateLinearMask }) != m0.end())
						m0[{diff, stateLinearMask}] += preComputedCor * cor2;
					else
					{
						m0[{diff, stateLinearMask}] = preComputedCor * cor2;

						IDLStateTypeBySheet newIDLLinearMask = IDLStateTypeBySheetFromStateTypeBySheet(stateLinearMask);
						for (int j = endr - 1; j > startr; j--)
						{
							IDLLinPropagateThroughRevNonLinearLayer(newIDLLinearMask);
							IDLLinPropagateThroughRevLinearLayer(newIDLLinearMask);
						}

						m1[{diff, stateLinearMask}] = { newIDLLinearMask, newIDLDiff };
					}
				}
			}
		}

		DiffAndLinearMaskPairs.clear();
		preComputedCors.clear();
		IDLLinearMaskAndDiffPairs.clear();

		for (auto& item : m0)
		{
			DiffAndLinearMaskPairs.push_back(item.first);
			preComputedCors.push_back(item.second);
			IDLLinearMaskAndDiffPairs.push_back(m1[item.first]);
		}

		endr -= 1;
	}

	double ExperimentalComputation(const HalfStateTypeBySheet& stateDiff, const HalfStateTypeBySheet& stateLinearMask, int startr, int endr, int totalRounds) const
	{
		random_device rd;
		mt19937 gen(rd());
		uniform_int_distribution<uint16_t> distribution(0, numeric_limits<uint16_t>::max());

		long long totalTestNum = 1ll << 28;
		long long zeroCount = 0;
		for (long long testIndex = 0; testIndex < totalTestNum; testIndex++)
		{
			StateTypeBySheet state0;
			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
				{
					state0[x][y] = distribution(gen);
					state0[x][y] ^= ((state0[x][y] ^ stateDiff[x][y]) << 16);
				}



			xoodooBase.ReducedRoundForRotDiffVerification(state0, startr, endr, totalRounds);

			HalfStateTypeBySheet outputStateDiff;
			for (int x = 0; x < 4; x++)
				for (int y = 0; y < 3; y++)
					outputStateDiff[x][y] = state0[x][y] ^ (state0[x][y] >> 16);

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

	double ComputeThreeRoundIDLCor(const HalfStateTypeBySheet& inputStateDiff, const HalfStateTypeBySheet& outputStateLinearMask, int startr, int endr, int totalRounds) const
	{
		IDLStateTypeBySheet IDLStateDiff = IDLStateTypeBySheetFromStateTypeBySheet(inputStateDiff);
		for (int i = startr; i < endr; i++)
		{
			int rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(6, i + 1);
			HalfLaneType rcDiff = internalDifferencePerRoundConstants[rcDiffIndex];
			IDLDiffPropagateThroughNonLinearLayer(IDLStateDiff);
			IDLDiffPropagateThroughLinearLayer(IDLStateDiff, rcDiff);
		}

		IDLStateTypeBySheet IDLStateLin = IDLStateTypeBySheetFromStateTypeBySheet(outputStateLinearMask);
		for (int i = endr; i > startr; i--)
		{
			IDLLinPropagateThroughRevNonLinearLayer(IDLStateLin);
			IDLLinPropagateThroughRevLinearLayer(IDLStateLin);
		}

		vector<pair<HalfStateTypeBySheet, HalfStateTypeBySheet>> DiffAndLinearMaskPairs = { {inputStateDiff, outputStateLinearMask} };
		vector<pair<IDLStateTypeBySheet, IDLStateTypeBySheet>> IDLLinearMaskAndDiffPairs = { {IDLStateLin, IDLStateDiff} };
		vector<double> preComputedCors = { 1.0 };

		// here, we actually call the backward division first, then we call the forward division. If the program stucks here, try to switch to call the backward divisions twice or call the forward divisions twice.
		IDLBackwardExpansionByOneRound(DiffAndLinearMaskPairs, preComputedCors, IDLLinearMaskAndDiffPairs, startr, endr, totalRounds, true);
		IDLForwardExpansionByOneRound(DiffAndLinearMaskPairs, preComputedCors, IDLLinearMaskAndDiffPairs, startr, endr, totalRounds, true);


		int pairsNum = DiffAndLinearMaskPairs.size();

		cout << __func__ << ": pairs num " << pairsNum << endl;

		double finalCor = 0.0;

		int threadsNum = omp_get_max_threads();

#pragma omp parallel for num_threads(threadsNum)
		for (int i = 0; i < pairsNum; i++)
		{
			cout << __func__ << ": pair index " << i << endl;
			const pair<HalfStateTypeBySheet, HalfStateTypeBySheet>& diffAndLinearMaskPair = DiffAndLinearMaskPairs[i];
			const HalfStateTypeBySheet& diff = diffAndLinearMaskPair.first;
			const HalfStateTypeBySheet& linearMask = diffAndLinearMaskPair.second;
			const IDLStateTypeBySheet& IDLLinearMask = IDLLinearMaskAndDiffPairs[i].first;
			const IDLStateTypeBySheet& IDLDiff = IDLLinearMaskAndDiffPairs[i].second;
			const double& preComputedCor = preComputedCors[i];

			if (preComputedCor == 0.0)
				continue;

			double dlctThroughNonLinearCor = 1.0;

			bool isValid = true;
			for (int x = 0; x < 4; x++)
			{
				for (int z = 0; z < 16; z++)
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

	double ComputeFourRoundIDLCor(const HalfStateTypeBySheet& inputStateDiff, const HalfStateTypeBySheet& outputStateLinearMask, int startr, int endr, int totalRounds) const
	{
		IDLStateTypeBySheet IDLStateDiff = IDLStateTypeBySheetFromStateTypeBySheet(inputStateDiff);
		for (int i = startr; i < endr; i++)
		{
			int rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(6, i + 1);
			HalfLaneType rcDiff = internalDifferencePerRoundConstants[rcDiffIndex];
			IDLDiffPropagateThroughNonLinearLayer(IDLStateDiff);
			IDLDiffPropagateThroughLinearLayer(IDLStateDiff, rcDiff);
		}

		IDLStateTypeBySheet IDLStateLin = IDLStateTypeBySheetFromStateTypeBySheet(outputStateLinearMask);
		for (int i = endr; i > startr; i--)
		{
			IDLLinPropagateThroughRevNonLinearLayer(IDLStateLin);
			IDLLinPropagateThroughRevLinearLayer(IDLStateLin);
		}

		int initialForwardEstimation = IDLFowardSingleStateDiffThrouChiEstimation(inputStateDiff, IDLStateLin);
		int initialBackwardEstimation = IDLBackwardSingleStateLinearMaskThroughRevChiEstimation(IDLStateDiff, outputStateLinearMask);



		vector<pair<HalfStateTypeBySheet, HalfStateTypeBySheet>> DiffAndLinearMaskPairs = { {inputStateDiff, outputStateLinearMask} };
		vector<pair<IDLStateTypeBySheet, IDLStateTypeBySheet>> IDLLinearMaskAndDiffPairs = { {IDLStateLin, IDLStateDiff} };
		vector<double> preComputedCors = { 1.0 };

		// here we call forward division once, and then call backward divisions once; so that the 4-round IDL correlation computation is reduced to multiple 2-round IDL correlation computations.
		IDLForwardExpansionByOneRound(DiffAndLinearMaskPairs, preComputedCors, IDLLinearMaskAndDiffPairs, startr, endr, totalRounds, true);
		IDLBackwardExpansionByOneRound(DiffAndLinearMaskPairs, preComputedCors, IDLLinearMaskAndDiffPairs, startr, endr, totalRounds, true);









		int pairsNum = DiffAndLinearMaskPairs.size();

		cout << __func__ << ": pairs num " << pairsNum << endl;

		double finalCor = 0.0;

		int threadsNum = omp_get_max_threads();

#pragma omp parallel for num_threads(threadsNum)
		for (int i = 0; i < pairsNum; i++)
		{
			cout << __func__ << ": pair index " << i << endl;
			const pair<HalfStateTypeBySheet, HalfStateTypeBySheet>& diffAndLinearMaskPair = DiffAndLinearMaskPairs[i];
			const HalfStateTypeBySheet& diff = diffAndLinearMaskPair.first;
			const HalfStateTypeBySheet& linearMask = diffAndLinearMaskPair.second;
			const IDLStateTypeBySheet& IDLLinearMask = IDLLinearMaskAndDiffPairs[i].first;
			const IDLStateTypeBySheet& IDLDiff = IDLLinearMaskAndDiffPairs[i].second;
			const double& preComputedCor = preComputedCors[i];

			if (preComputedCor == 0.0)
				continue;

			int forwardEstimation = IDLFowardSingleStateDiffThrouChiEstimation(diff, IDLLinearMask);
			int backwardEstimation = IDLBackwardSingleStateLinearMaskThroughRevChiEstimation(IDLDiff, linearMask);


			if (forwardEstimation == -1 || backwardEstimation == -1)
			{
				// cout << __func__ << ": cor is determined as 0" << endl;
				continue;
			}

			cout << __func__ << ": pair index " << i << endl;
			cout << __func__ << ": forward estimation " << forwardEstimation << endl;
			cout << __func__ << ": backward estimation " << backwardEstimation << endl;
			int diffActiveColumnsNum = idBase.getActiveColumnsNum(diff);
			int linActiveColumnsNum = iLinearBase.getActiveColumnsNum(linearMask);
			cout << "Difference active columns num " << diffActiveColumnsNum << endl;
			cout << "Linear mask active columns num " << linActiveColumnsNum << endl;




			cout << "Use milp models to compute." << endl;


#pragma omp critical 
			{
				if (diffActiveColumnsNum > linActiveColumnsNum)
					finalCor += preComputedCor * ComputeTwoRoundIDLCorBackwardAutomatically(diff, linearMask, startr, endr, totalRounds);
				else
					finalCor += preComputedCor * ComputeTwoRoundIDLCorForwardAutomatically(diff, linearMask, startr, endr, totalRounds);
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

	void MilpPropagateThroughChi(GRBModel& model, IDStateMILPTypeBySheet& stateDiffVars) const
	{
		for (int x = 0; x < 4; x++)
			for (int z = 0; z < 16; z++)
				MilpPropagateThroughChiColumn(model, stateDiffVars[x][0][z], stateDiffVars[x][1][z], stateDiffVars[x][2][z]);
	}

	void MilpPropagateThroughRevChi(GRBModel& model, IDStateMILPTypeBySheet& stateMaskVars) const
	{
		MilpPropagateThroughChi(model, stateMaskVars);
	}

	double ComputeTwoRoundIDLCorForwardAutomatically(const HalfStateTypeBySheet& inputStateDiff, const HalfStateTypeBySheet& outputStateLinearMask, int startr, int endr, int totalRounds) const
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_PoolSolutions, 200000000);
		env.set(GRB_IntParam_Threads, milpThreads);

		GRBModel model = GRBModel(env);

		IDStateMILPTypeBySheet stateDiffVars = idBase.MilpAddStateVars(model);

		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					model.addConstr(stateDiffVars[x][y][z] == xoodooBase.GetBit(inputStateDiff, x, y, z));

		IDStateMILPTypeBySheet inputStateDiffVarsForRound0Chi = stateDiffVars;
		idBase.MilpPropagateThroughChi(model, stateDiffVars);
		IDStateMILPTypeBySheet outputStateDiffVarsForRound0Chi = stateDiffVars;

		int internalDifferenceIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + 1);
		idBase.MilpPropagateThroughLinearLayer(model, stateDiffVars, internalDifferencePerRoundConstants[internalDifferenceIndex]);

		GRBLinExpr objExp = 0;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					objExp += stateDiffVars[x][y][z];

		IDStateMILPTypeBySheet inputStateDiffVarsForRound1Chi = stateDiffVars;

		MilpPropagateThroughChi(model, stateDiffVars);
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					model.addConstr(stateDiffVars[x][y][z] == xoodooBase.GetBit(outputStateLinearMask, x, y, z));

		ILINStateMILPTypeBySheet outputStateLinearMaskVarsForRound1Chi = stateDiffVars;
		model.setObjective(objExp, GRB_MINIMIZE);



		model.optimize();

		int solCount = model.get(GRB_IntAttr_SolCount);
		cout << "Get " << solCount << " two round forward RDL trails." << endl;

		if (solCount == 0)
			return 0.0;
		else
		{
			double IDLCor = 0.0;
			for (int solIndex = 0; solIndex < solCount; solIndex++)
			{
				model.set(GRB_IntParam_SolutionNumber, solIndex);

				HalfStateTypeBySheet inputStateDiffForRound0Chi = idBase.MilpReadStateVars(inputStateDiffVarsForRound0Chi);
				HalfStateTypeBySheet outputStateDiffForRound0Chi = idBase.MilpReadStateVars(outputStateDiffVarsForRound0Chi);
				HalfStateTypeBySheet inputStateDiffForRound1Chi = idBase.MilpReadStateVars(inputStateDiffVarsForRound1Chi);
				HalfStateTypeBySheet outputStateLinearMaskForRound1Chi = iLinearBase.MilpReadStateVars(outputStateLinearMaskVarsForRound1Chi);
				auto ddtProb = idBase.getDDTProb();

				double IDLCorOneTrail = 1.0;
				for (int x = 0; x < 4; x++)
					for (int z = 0; z < 16; z++)
					{
						int inputDiffForRound0Chi = xoodooBase.GetColumn(inputStateDiffForRound0Chi, x, z);
						int outputDiffForRound0Chi = xoodooBase.GetColumn(outputStateDiffForRound0Chi, x, z);
						int inputDiffForRound1Chi = xoodooBase.GetColumn(inputStateDiffForRound1Chi, x, z);
						int outputLinearMaskForRound1Chi = xoodooBase.GetColumn(outputStateLinearMaskForRound1Chi, x, z);
						IDLCorOneTrail *= ddtProb[inputDiffForRound0Chi][outputDiffForRound0Chi] * dlctCor[inputDiffForRound1Chi][outputLinearMaskForRound1Chi];
					}

				IDLCor += IDLCorOneTrail;
			}


			return IDLCor;
		}
	}

	double ComputeTwoRoundIDLCorBackwardAutomatically(const HalfStateTypeBySheet& inputStateDiff, const HalfStateTypeBySheet& outputStateLinearMask, int startr, int endr, int totalRounds) const
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_PoolSolutions, 200000000);
		env.set(GRB_IntParam_Threads, milpThreads);

		GRBModel model = GRBModel(env);

		ILINStateMILPTypeBySheet stateMaskVars = iLinearBase.MilpAddStateVars(model);

		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					model.addConstr(stateMaskVars[x][y][z] == xoodooBase.GetBit(outputStateLinearMask, x, y, z));

		ILINStateMILPTypeBySheet outputStateLinearMaskVarsForRound1Chi = stateMaskVars;
		iLinearBase.MilpPropagateThroughRevNonLinearLayer(model, stateMaskVars);
		ILINStateMILPTypeBySheet inputStateLinearMaskVarsForRound1Chi = stateMaskVars;

		iLinearBase.MilpPropagateThroughRevLinearLayer(model, stateMaskVars);

		GRBLinExpr objExp = 0;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					objExp += stateMaskVars[x][y][z];



		ILINStateMILPTypeBySheet outputStateLinearMaskVarsForRound0Chi = stateMaskVars;

		MilpPropagateThroughRevChi(model, stateMaskVars);
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					model.addConstr(stateMaskVars[x][y][z] == xoodooBase.GetBit(inputStateDiff, x, y, z));

		IDStateMILPTypeBySheet inputStateDiffVarsForRound0Chi = stateMaskVars;
		model.setObjective(objExp, GRB_MINIMIZE);



		model.optimize();

		int solCount = model.get(GRB_IntAttr_SolCount);
		cout << "Get " << solCount << " two round backward RDL trails." << endl;

		if (solCount == 0)
			return 0.0;
		else
		{
			double IDLCor = 0.0;

			for (int solIndex = 0; solIndex < solCount; solIndex++)
			{
				HalfStateTypeBySheet inputStateDiffForRound0Chi = idBase.MilpReadStateVars(inputStateDiffVarsForRound0Chi);
				HalfStateTypeBySheet outputStateLinearMaskForRound0Chi = iLinearBase.MilpReadStateVars(outputStateLinearMaskVarsForRound0Chi);
				HalfStateTypeBySheet inputStateLinearMaskForRound1Chi = iLinearBase.MilpReadStateVars(inputStateLinearMaskVarsForRound1Chi);
				HalfStateTypeBySheet outputStateLinearMaskForRound1Chi = iLinearBase.MilpReadStateVars(outputStateLinearMaskVarsForRound1Chi);

				double IDLCorOneTrail = 1.0;
				auto latCor2 = iLinearBase.getLATCor2();
				for (int x = 0; x < 4; x++)
					for (int z = 0; z < 16; z++)
					{
						int inputDiffForRound0Chi = xoodooBase.GetColumn(inputStateDiffForRound0Chi, x, z);
						int outputLinearMaskForRound0Chi = xoodooBase.GetColumn(outputStateLinearMaskForRound0Chi, x, z);
						int inputLinearMaskForRound1Chi = xoodooBase.GetColumn(inputStateLinearMaskForRound1Chi, x, z);
						int outputLinearMaskForRound1Chi = xoodooBase.GetColumn(outputStateLinearMaskForRound1Chi, x, z);

						IDLCorOneTrail *= dlctCor[inputDiffForRound0Chi][outputLinearMaskForRound0Chi] * latCor2[inputLinearMaskForRound1Chi][outputLinearMaskForRound1Chi];
					}

				int rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, endr);
				HalfLaneType rcDiff = internalDifferencePerRoundConstants[rcDiffIndex];
				iLinearBase.PropagateThroughLota(inputStateLinearMaskForRound1Chi, IDLCorOneTrail, rcDiff);

				IDLCor += IDLCorOneTrail;
			}


			return IDLCor;
		}


	}

	vector<pair<HalfStateTypeBySheet, double>> IDLForwardSingleStateDiffThroughChiAutomatically(const HalfStateTypeBySheet& inputStateDiff, const HalfStateTypeBySheet& outputStateLinMask, int startr, int endr, int totalRounds) const
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_IntParam_Threads, milpThreads);

		GRBModel model = GRBModel(env);

		IDStateMILPTypeBySheet stateDiffVars = idBase.MilpAddStateVars(model);

		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					model.addConstr(stateDiffVars[x][y][z] == xoodooBase.GetBit(inputStateDiff, x, y, z));

		idBase.MilpPropagateThroughChi(model, stateDiffVars);

		IDStateMILPTypeBySheet targetSolVars = stateDiffVars;

		int internalDifferenceIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, startr + 1);
		idBase.MilpPropagateThroughLinearLayer(model, stateDiffVars, internalDifferencePerRoundConstants[internalDifferenceIndex]);

		MilpPropagateThroughChi(model, stateDiffVars);



		ILINStateMILPTypeBySheet stateMaskVars = iLinearBase.MilpAddStateVars(model);
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					model.addConstr(stateMaskVars[x][y][z] == xoodooBase.GetBit(outputStateLinMask, x, y, z));

		for (int i = endr; i > startr + 1; i--)
		{
			iLinearBase.MilpPropagateThroughRevNonLinearLayer(model, stateMaskVars);
			iLinearBase.MilpPropagateThroughRevLinearLayer(model, stateMaskVars);
		}

		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					model.addConstr(stateMaskVars[x][y][z] == stateDiffVars[x][y][z]);

		GRBLinExpr objExp = 0;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					objExp += targetSolVars[x][y][z];

		vector<HalfStateTypeBySheet> targetSolutions;
		HalfStateForwardExpandCallback cb = HalfStateForwardExpandCallback(xoodooBase, targetSolVars, targetSolutions);
		model.setCallback(&cb);

		model.setObjective(objExp, GRB_MAXIMIZE);

		model.optimize();


		vector<pair<HalfStateTypeBySheet, double>> outputStateDiffAndProbs;
		auto ddtProb = idBase.getDDTProb();

		for (auto& targetSolution : targetSolutions)
		{
			HalfStateTypeBySheet outputStateDiff = targetSolution;
			double prob = 1.0;
			for (int x = 0; x < 4; x++)
				for (int z = 0; z < 16; z++)
				{
					int inputDiffColumn = xoodooBase.GetColumn(inputStateDiff, x, z);
					int outputDiffColumn = xoodooBase.GetColumn(outputStateDiff, x, z);
					prob *= ddtProb[inputDiffColumn][outputDiffColumn];
				}

			outputStateDiffAndProbs.push_back({ outputStateDiff, prob });
		}

		return outputStateDiffAndProbs;
	}

	vector<pair<HalfStateTypeBySheet, double>> IDLBackwardSingleStateLinearMaskThroughRevChiAutomatically(const HalfStateTypeBySheet& inputStateDiff, const HalfStateTypeBySheet& outputStateLinMask, int startr, int endr, int totalRounds) const
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_IntParam_Threads, milpThreads);

		GRBModel model = GRBModel(env);

		ILINStateMILPTypeBySheet stateMaskVars = idBase.MilpAddStateVars(model);

		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					model.addConstr(stateMaskVars[x][y][z] == xoodooBase.GetBit(outputStateLinMask, x, y, z));

		iLinearBase.MilpPropagateThroughRevChi(model, stateMaskVars);

		IDStateMILPTypeBySheet targetSolVars = stateMaskVars;

		iLinearBase.MilpPropagateThroughRevLinearLayer(model, stateMaskVars);

		MilpPropagateThroughRevChi(model, stateMaskVars);



		IDStateMILPTypeBySheet stateDiffVars = idBase.MilpAddStateVars(model);
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					model.addConstr(stateDiffVars[x][y][z] == xoodooBase.GetBit(inputStateDiff, x, y, z));

		for (int i = startr; i < endr - 1; i++)
		{
			idBase.MilpPropagateThroughNonLinearLayer(model, stateDiffVars);

			int rcDiffIndex = xoodooBase.ComputeRoundConstantIndex(totalRounds, i + 1);
			idBase.MilpPropagateThroughLinearLayer(model, stateDiffVars, internalDifferencePerRoundConstants[rcDiffIndex]);
		}

		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					model.addConstr(stateMaskVars[x][y][z] == stateDiffVars[x][y][z]);

		GRBLinExpr objExp = 0;
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				for (int z = 0; z < 16; z++)
					objExp += targetSolVars[x][y][z];

		vector<HalfStateTypeBySheet> targetSolutions;
		HalfStateBackwardExpandCallback cb = HalfStateBackwardExpandCallback(xoodooBase, targetSolVars, targetSolutions);
		model.setCallback(&cb);

		model.setObjective(objExp, GRB_MAXIMIZE);

		model.optimize();


		vector<pair<HalfStateTypeBySheet, double>> inputStateMaskAndCor2s;
		auto latCor2 = iLinearBase.getLATCor2();

		for (auto& targetSolution : targetSolutions)
		{
			HalfStateTypeBySheet inputStateMask = targetSolution;
			double cor2 = 1.0;
			for (int x = 0; x < 4; x++)
				for (int z = 0; z < 16; z++)
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


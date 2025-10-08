#pragma once

#include<iostream>
#include<array>
#include<vector>
#include<iomanip>

using namespace std;

class BinaryMatrix64 {
public:
	std::array<uint64_t, 64> rows;

	BinaryMatrix64() {
		rows.fill(0);
	}

	// Portable popcount (count number of 1 bits)
	static int popcount(uint64_t x) {
		int count = 0;
		while (x) {
			count += x & 1;
			x >>= 1;
		}
		return count;
	}

	// Generate the matrix L for ASCON's linear layer
	// For operation: y = x ^ (x >> a) ^ (x >> b)
	static BinaryMatrix64 generateLinearMatrix(int a, int b) {
		BinaryMatrix64 L;

		// Identity part (from x)
		for (int i = 0; i < 64; i++) {
			L.rows[i] |= (1ULL << i);
		}

		// Rotation by a part (from x >> a)
		for (int i = 0; i < 64; i++) {
			int j = (i + a) % 64;  // Source bit position
			L.rows[i] |= (1ULL << j);
		}

		// Rotation by b part (from x >> b)
		for (int i = 0; i < 64; i++) {
			int j = (i + b) % 64;  // Source bit position
			L.rows[i] |= (1ULL << j);
		}

		return L;
	}

	// Compute inverse matrix using Gaussian elimination
	BinaryMatrix64 inverse() const {
		// Create augmented matrix [L | I]
		std::array<uint64_t, 64> left = rows;  // Copy of L
		std::array<uint64_t, 64> right;         // Identity matrix

		// Initialize right side as identity matrix
		for (int i = 0; i < 64; i++) {
			right[i] = 1ULL << i;
		}

		// Gaussian elimination
		for (int col = 0; col < 64; col++) {
			// Find pivot row
			int pivot = -1;
			for (int row = col; row < 64; row++) {
				if ((left[row] >> col) & 1) {
					pivot = row;
					break;
				}
			}

			if (pivot == -1) {
				throw std::runtime_error("Matrix is not invertible");
			}

			// Swap rows if necessary
			if (pivot != col) {
				std::swap(left[pivot], left[col]);
				std::swap(right[pivot], right[col]);
			}

			// Eliminate column in all other rows
			for (int row = 0; row < 64; row++) {
				if (row != col && ((left[row] >> col) & 1)) {
					left[row] ^= left[col];
					right[row] ^= right[col];
				}
			}
		}

		// The right side now contains the inverse
		BinaryMatrix64 result;
		result.rows = right;
		return result;
	}

	// Apply matrix to a 64-bit value
	uint64_t apply(uint64_t x) const {
		uint64_t result = 0;

		for (int i = 0; i < 64; i++) {
			// Count parity of (rows[i] & x)
			uint64_t temp = rows[i] & x;
			int parity = popcount(temp) & 1;
			result |= (uint64_t(parity) << i);
		}

		return result;
	}

	// Alternative apply using parity trick (no popcount needed)
	uint64_t applyFast(uint64_t x) const {
		uint64_t result = 0;

		for (int i = 0; i < 64; i++) {
			// Compute parity using XOR reduction
			uint64_t temp = rows[i] & x;
			temp ^= temp >> 32;
			temp ^= temp >> 16;
			temp ^= temp >> 8;
			temp ^= temp >> 4;
			temp ^= temp >> 2;
			temp ^= temp >> 1;
			result |= ((temp & 1) << i);
		}

		return result;
	}

	// Print matrix (for debugging, shows first 8x8 block)
	void print() const {
		std::cout << "First 8x8 block:" << std::endl;
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 8; j++) {
				std::cout << ((rows[i] >> j) & 1);
			}
			std::cout << "..." << std::endl;
		}
		std::cout << "..." << std::endl;
	}
};

using RowType = uint64_t;
using StateType = array<RowType, 5>;

class MyException : public std::exception {
public:
	MyException(const std::string& message) : msg_(message) {}


	virtual const char* what() const noexcept override {
		return msg_.c_str();
	}

private:
	std::string msg_;
};

RowType RotateLeft(const RowType value, int shift) {
	return (value << shift) | (value >> (64 - shift));
}

RowType RotateRight(const RowType value, int shift) {
	return (value >> shift) | (value << (64 - shift));
}

class AsconBase {
private:
	const array<RowType, 12> roundConstants = { 0xf0, 0xe1, 0xd2, 0xc3, 0xd4, 0xa5, 0x96, 0x87, 0x78, 0x69, 0x5a, 0x4b };

	const std::array<std::array<int, 2>, 5> rotateNums = { {
	{19, 28},
	{61, 39},
	{1, 6},
	{10, 17},
	{7, 41}} };

	array<BinaryMatrix64, 5> inverseLinearMatrices;

	const array<int, 32> sBoxLUT = { 0x4, 0xb, 0x1f, 0x14, 0x1a, 0x15, 0x9, 0x2, 0x1b, 0x5, 0x8, 0x12, 0x1d, 0x3, 0x6, 0x1c, 0x1e, 0x13, 0x7, 0xe, 0x0, 0xd, 0x11, 0x18, 0x10,0xc,0x1,0x19,0x16,0xa,0xf,0x17 };

	void GenerateInverseLinearMatrices() 
	{
		for (int i = 0; i < 5; i++)
		{
			BinaryMatrix64 L = BinaryMatrix64::generateLinearMatrix(rotateNums[i][0], rotateNums[i][1]);
			BinaryMatrix64 L_inv = L.inverse();
			inverseLinearMatrices[i] = L_inv;
		}
	}

public:

	AsconBase() {
		GenerateInverseLinearMatrices();
	};

	RowType GetRoundConstant(int index) const
	{
		return roundConstants[index];
	}

	int GetRotationNum(int index0, int index1) const
	{
		return rotateNums[index0][index1];
	}

	void SetBit(StateType& state, int row, int col, bool value) const
	{
		if (value)
		{
			state[row] |= (1ull << col);
		}
		else
			state[row] &= ~(1ull << col);
	}

	bool GetBit(const StateType& state, int row, int col) const
	{
		return (state[row] >> col) & 1ull;
	}

	void SetColumn(StateType& state, int col, uint8_t value) const
	{
		for (int row = 0; row < 5; row++)
		{
			bool bit = (value >> row) & 1ull;
			SetBit(state, 4-row, col, bit);
		}
	}

	uint8_t GetColumn(const StateType& state, int col) const
	{
		uint8_t value = 0;
		for (int row = 0; row < 5; row++)
		{
			bool bit = GetBit(state, row, col);
			value |= (bit << (4-row));
		}
		return value;
	}

	void RotateRows(StateType& state) const
	{
		for (int row = 0; row < 5; row++)
		{
			state[row] = state[row] ^ RotateRight(state[row], rotateNums[row][0]) ^ RotateRight(state[row], rotateNums[row][1]);
		}
	}

	void RotateRowsTranspose(StateType& state) const
	{
		for (int row = 0; row < 5; row++)
		{
			state[row] = state[row] ^ RotateLeft(state[row], rotateNums[row][0]) ^ RotateLeft(state[row], rotateNums[row][1]);
		}
	}

	void RotateRowsInverse(StateType& state) const
	{
		for (int row = 0; row < 5; row++)
		{
			state[row] = inverseLinearMatrices[row].applyFast(state[row]);
		}
	}

	void LinearLayer(StateType& state, const RowType rc) const
	{
		RotateRows(state);
		state[2] ^= rc;
	}

	void LinearLayerInverse(StateType& state, const RowType rc) const
	{
		state[2] ^= rc;
		RotateRowsInverse(state);
	}



	uint8_t Sbox(uint8_t input) const
	{
		return sBoxLUT[input];
	}

	void SboxLayer(StateType& state) const
	{
		for (int col = 0; col < 64; col++)
		{
			uint8_t input = GetColumn(state, col);
			uint8_t output = Sbox(input);
			SetColumn(state, col, output);
		}
	}

	void FullRound(StateType& state, const RowType rc) const
	{
		// first is sbox layer, then is linear layer, and final is the add round constant
		SboxLayer(state);
		LinearLayer(state, rc);
	}

	void Permute(StateType& state, int rounds) const
	{
		for (int i = 0; i < rounds; i++)
		{
			FullRound(state, roundConstants[i + 1]);
		}
	}

	void ReducedRoundsForVerification(StateType& state, int startr, int endr) const
	{
		for (int i = startr + 1; i <= endr + 1; i++)
		{
			FullRound(state, roundConstants[i]);
		}
	}


	void OutputState(ostream& os, const StateType& state) const
	{
		for (int row = 0; row < 5; row++)
		{
			os << "Row " << row << ": ";
			for (int col = 63; col >= 0; col--)
			{
				os << ((state[row] >> col) & 1ull);
			}
			os << endl;
		}
	}

	void OutputStateAs3DArr(ostream& os, const StateType& state) const
	{
		for (int row = 0; row < 5; row++)
		{
			for (int col = 0; col < 64; col++)
			{
				os << ((state[row] >> col) & 1ull);
				os << ", ";
			}
			os << endl;
		}
	}

	virtual ~AsconBase() = default;
};

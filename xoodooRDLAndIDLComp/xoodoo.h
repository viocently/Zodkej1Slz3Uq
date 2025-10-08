#pragma once
#include<iostream>
#include<array>
#include<vector>
#include<iomanip>

using namespace std;

using LaneType = uint32_t;
using SheetType = array<LaneType, 3>;
using PlaneType = array<LaneType, 4>;
using StateTypeBySheet = array<SheetType, 4>;
using StateTypeByPlane = array<PlaneType, 3>;

using HalfLaneType = uint16_t;
using HalfSheetType = array<HalfLaneType, 3>;
using HalfPlaneType = array<HalfLaneType, 4>;
using HalfStateTypeBySheet = array<HalfSheetType, 4>;
using HalfStateTypeByPlane = array<HalfPlaneType, 3>;


class MyException : public std::exception {
public:
	MyException(const std::string& message) : msg_(message) {}


	virtual const char* what() const noexcept override {
		return msg_.c_str();
	}

private:
	std::string msg_;
};


LaneType RotateLeft(const LaneType lane, int shift) {
	return (shift == 0) ? lane : ((lane << shift) | (lane >> (32 - shift)));
}

HalfLaneType RotateLeft(const HalfLaneType lane, int shift) {
	return (shift == 0) ? lane : ((lane << shift) | (lane >> (16 - shift)));
}

LaneType RotateRight(const LaneType lane, int shift) {
	return (shift == 0) ? lane : ((lane >> shift) | (lane << (32 - shift)));
}

HalfLaneType RotateRight(const HalfLaneType lane, int shift) {
	return (shift == 0) ? lane : ((lane >> shift) | (lane << (16 - shift)));
}

class XoodooBase
{
private:


	const array<LaneType, 12> RoundConstants = {
		0x00000058, 0x00000038, 0x000003C0, 0x000000D0,
		0x00000120, 0x00000014, 0x00000060, 0x0000002C,
		0x00000380, 0x000000F0, 0x000001A0, 0x00000012
	};

	const array<int, 8> sboxLUT = {0, 3, 6, 1, 5, 4, 2, 7};

public:
	LaneType GetRoundConstant(int index) const {
		return RoundConstants[index];
	}
	
	int ComputeRoundConstantIndex(int rounds, int round) const {
		return 12 - (rounds - round);
	}

	void SetBit(StateTypeBySheet& state, int x, int y, int z, bool value) const {
		if(value)
			state[x][y] |= (1u << z);
		else
			state[x][y] &= ~(1u << z);
	}

	void SetBit(HalfStateTypeBySheet& state, int x, int y, int z, bool value) const {
		if (value)
			state[x][y] |= (1u << z);
		else
			state[x][y] &= ~(1u << z);
	}

	bool GetBit(const StateTypeBySheet& state, int x, int y, int z) const {
		return (state[x][y] >> z) & 1;
	}

	bool GetBit(const HalfStateTypeBySheet& state, int x, int y, int z) const {
		return (state[x][y] >> z) & 1;
	}

	void SetColumn(StateTypeBySheet& state, int x, int z, uint8_t value) const {
		SetBit(state, x, 0, z, (value >> 0) & 1);
		SetBit(state, x, 1, z, (value >> 1) & 1);
		SetBit(state, x, 2, z, (value >> 2) & 1);
	}

	void SetColumn(HalfStateTypeBySheet& state, int x, int z, uint8_t value) const {
		SetBit(state, x, 0, z, (value >> 0) & 1);
		SetBit(state, x, 1, z, (value >> 1) & 1);
		SetBit(state, x, 2, z, (value >> 2) & 1);
	}
	
	uint8_t GetColumn(const StateTypeBySheet& state, int x, int z) const {
		return (GetBit(state, x, 0, z) << 0) | (GetBit(state, x, 1, z) << 1) | (GetBit(state, x, 2, z) << 2);
	}

	uint8_t GetColumn(const HalfStateTypeBySheet& state, int x, int z) const {
		return (GetBit(state, x, 0, z) << 0) | (GetBit(state, x, 1, z) << 1) | (GetBit(state, x, 2, z) << 2);
	}
	



	void Theta(StateTypeBySheet& state) const {
		PlaneType parity = { 0 };
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				parity[x] ^= state[x][y];

		for (int x = 0; x < 4; x++) {
			LaneType one = RotateLeft(parity[(x + 3) % 4], 5);
			LaneType two = RotateLeft(parity[(x + 3) % 4], 14);
			for (int y = 0; y < 3; y++)
				state[x][y] ^= one ^ two;
		}
	}

	void ThetaTranspose(StateTypeBySheet& state) const
	{
		PlaneType parity = { 0 };
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				parity[x] ^= state[x][y];

		for (int x = 0; x < 4; x++) {
			LaneType one = RotateRight(parity[(x + 1) % 4], 5);
			LaneType two = RotateRight(parity[(x + 1) % 4], 14);
			for (int y = 0; y < 3; y++)
				state[x][y] ^= one ^ two;
		}
	}

	void ThetaTranspose(HalfStateTypeBySheet& state) const
	{
		HalfPlaneType parity = { 0 };
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				parity[x] ^= state[x][y];

		for (int x = 0; x < 4; x++) {
			HalfLaneType one = RotateRight(parity[(x + 1) % 4], 5);
			HalfLaneType two = RotateRight(parity[(x + 1) % 4], 14);
			for (int y = 0; y < 3; y++)
				state[x][y] ^= one ^ two;
		}
	}

	void RhoWest(StateTypeBySheet& state) const {
		StateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = tmp[(x + 3) % 4][1];

		for (int x = 0; x < 4; x++)
			state[x][2] = RotateLeft(tmp[x][2], 11);
	}

	void RhoWestTranspose(StateTypeBySheet & state) const {
		StateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = tmp[(x + 1) % 4][1];

		for(int x = 0; x < 4; x++)
			state[x][2] = RotateRight(tmp[x][2], 11);
	}

	void RhoWestTranspose(HalfStateTypeBySheet& state) const {
		HalfStateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = tmp[(x + 1) % 4][1];

		for (int x = 0; x < 4; x++)
			state[x][2] = RotateRight(tmp[x][2], 11);
	}

	void Lota(StateTypeBySheet& state, const LaneType rc) const {
		state[0][0] ^= rc;
	}

	

	void Chi(StateTypeBySheet& state) const {
		StateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			for (int y = 0; y < 3; y++)
				state[x][y] = tmp[x][y] ^ ((~tmp[x][(y + 1) % 3]) & tmp[x][(y + 2) % 3]);
	}

	uint8_t ChiColumnByOperation(uint8_t col) const {
		int a0 = (col >> 0) & 1;
		int a1 = (col >> 1) & 1;
		int a2 = (col >> 2) & 1;
		int b0 = a0 ^ ((a1 ^ 1) & a2);
		int b1 = a1 ^ ((a2 ^ 1) & a0);
		int b2 = a2 ^ ((a0 ^ 1) & a1);
		return (b2 << 2) | (b1 << 1) | b0;
	}
	
	uint8_t ChiColumnByLUT(uint8_t& col) const {
		return sboxLUT[col];
	}

	void RhoEast(StateTypeBySheet& state) const {
		StateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = RotateLeft(tmp[x][1], 1);

		for (int x = 0; x < 4; x++)
			state[x][2] = RotateLeft(tmp[(x + 2) % 4][2], 8);
	}

	void RhoEastTranspose(StateTypeBySheet& state) const {
		StateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = RotateRight(tmp[x][1], 1);

		for (int x = 0; x < 4; x++)
			state[x][2] = RotateRight(tmp[(x + 2) % 4][2], 8);
	}

	void RhoEastTranspose(HalfStateTypeBySheet& state) const {
		HalfStateTypeBySheet tmp(state);
		for (int x = 0; x < 4; x++)
			state[x][1] = RotateRight(tmp[x][1], 1);

		for (int x = 0; x < 4; x++)
			state[x][2] = RotateRight(tmp[(x + 2) % 4][2], 8);
	}


	void FullRound(StateTypeBySheet& state, const LaneType rc) const {
		Theta(state);
		RhoWest(state);
		Lota(state, rc);
		Chi(state);
		RhoEast(state);
	}

	void ReducedRoundForVerification(StateTypeBySheet& state, int startr, int endr, int totalRounds) const {
		NonLinearLayer(state);
		for (int i = startr+1; i <= endr; i++)
		{
			LaneType rc = GetRoundConstant(ComputeRoundConstantIndex(totalRounds, i));
			LinearLayer(state, rc);
			NonLinearLayer(state);
		}
	}

	// Through Chi, and then linear layer, stop at the input of Chi
	void ReducedRoundForRotDiffVerification(StateTypeBySheet& state, int startr, int endr, int totalRounds) const {
		NonLinearLayer(state);
		for (int i = startr + 1; i <= endr; i++)
		{
			LaneType rc = GetRoundConstant(ComputeRoundConstantIndex(totalRounds, i));
			LinearLayer(state, rc);
			NonLinearLayer(state);
		}

		// LaneType rc = GetRoundConstant(ComputeRoundConstantIndex(totalRounds, endr));
		// LinearLayer(state, rc);
	}

	

	void HalfRoundFromChi(StateTypeBySheet& state) const {
		Chi(state);
		RhoEast(state);
	}

	void HalftRoundToChi(StateTypeBySheet& state, const LaneType rc) const {
		Theta(state);
		RhoWest(state);
		Lota(state, rc);
		Chi(state);
	}
	
	void LinearLayer(StateTypeBySheet& state, const LaneType rc) const {
		RhoEast(state);
		Theta(state);
		RhoWest(state);
		Lota(state, rc);
	}


	void NonLinearLayer(StateTypeBySheet& state) const {
		Chi(state);
	}

	void Permute(StateTypeBySheet& state, int rounds) const {
		for (int round = 0; round < rounds; round++) {
			FullRound(state, RoundConstants[ComputeRoundConstantIndex(rounds, round)]);
		}
	}

	void OutputState(ostream & os, const StateTypeBySheet& state) const {
		for(int y = 0; y < 3; y++)
			{
			for(int x = 0; x < 4; x++)
			{
				for(int z = 0; z < 32; z++)
				{
					os << setw(2) << GetBit(state, x, y, z);
				}
				os << "|";
			}
			os << endl;
		}

		os << endl;
	}

	void OutputState(ostream& os, const HalfStateTypeBySheet& state) const {
		for (int y = 0; y < 3; y++)
		{
			for (int x = 0; x < 4; x++)
			{
				for (int z = 0; z < 16; z++)
				{
					os << setw(2) << GetBit(state, x, y, z);
				}
				os << "|";
			}
			os << endl;
		}

		os << endl;
	}

};


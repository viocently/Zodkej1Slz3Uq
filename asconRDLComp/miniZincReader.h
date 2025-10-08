#pragma once
#include<iostream>
#include<vector>
#include"ascon.h"

using namespace std;

class MiniZincReader
{
private:
	const AsconBase& asconBase;

public:
	MiniZincReader(const AsconBase& asconBase) : asconBase(asconBase) {}

	StateType ReadStateFrom3DArray(const vector<int>& arr) const
	{
		StateType state;
		for (int row = 0; row < 5; row++)
		{
			for (int col = 0; col < 64; col++)
			{

				asconBase.SetBit(state, row, col, arr[64 * row + col]);
			}
		}

		return state;
	}

	vector<StateType> ReadStatesFrom3DArrays(const vector<vector<int>>& arrs) const
	{
		vector<StateType> states;
		for (const auto& arr : arrs)
		{
			states.push_back(ReadStateFrom3DArray(arr));
		}
		return states;
	}

};

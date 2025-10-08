#pragma once
#include<iostream>
#include<vector>
#include"xoodoo.h"
#include"differential-linear.h"

using namespace std;

class MiniZincReader
{
private:
	const XoodooBase& xoodooBase;
	
public:
	MiniZincReader(const XoodooBase& xoodooBase) : xoodooBase(xoodooBase) {}
	
	StateTypeBySheet ReadStateFrom3DArray(const vector<int> & arr) const
	{
		StateTypeBySheet state;
		for(int x = 0; x < 4; x++)
		{
			for(int y = 0; y < 3; y++)
			{
				for(int z = 0; z < 32; z++)
				{
					xoodooBase.SetBit(state, x, y, z, arr[x * 3 * 32 + y * 32 + z]);
				}
			}
		}
		return state;
	}

	HalfStateTypeBySheet ReadHalfStateFrom3DArray(const vector<int>& arr) const
	{
		HalfStateTypeBySheet state;
		for (int x = 0; x < 4; x++)
		{
			for (int y = 0; y < 3; y++)
			{
				for (int z = 0; z < 16; z++)
				{
					xoodooBase.SetBit(state, x, y, z, arr[x * 3 * 16 + y * 16 + z]);
				}
			}
		}
		return state;
	}

	RDLStateTypeBySheet ReadRDLStateFrom3DArray(const vector<int>& arr) const
	{
		RDLStateTypeBySheet state;
		for (int x = 0; x < 4; x++)
		{
			for (int y = 0; y < 3; y++)
			{
				for (int z = 0; z < 32; z++)
				{
					state[x][y][z] = arr[x * 3 * 32 + y * 32 + z];
				}
			}
		}
		return state;
	}

	vector<StateTypeBySheet> ReadStatesFrom3DArrays(const vector<vector<int>>& arrs) const
	{
		vector<StateTypeBySheet> states;
		for(const auto& arr : arrs)
		{
			states.push_back(ReadStateFrom3DArray(arr));
		}
		return states;
	}

	vector<RDLStateTypeBySheet> ReadRDLStatesFrom3DArrays(const vector<vector<int>>& arrs) const
	{
		vector<RDLStateTypeBySheet> states;
		for (const auto& arr : arrs)
		{
			states.push_back(ReadRDLStateFrom3DArray(arr));
		}
		return states;
	}
};

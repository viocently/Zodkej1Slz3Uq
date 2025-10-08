#pragma once
#include<vector>
#include"xoodoo.h"
#include"internalDifferential.h"
#include"differential.h"
#include"linear.h"
#include"internalLinear.h"
#include"gurobi_c++.h"


class FullStateBackwardExpandCallback : public GRBCallback
{
public:
	const XoodooBase& xoodooBase;
	vector<StateTypeBySheet>* pSolutions = NULL;
	LINStateMILPTypeBySheet solVars;


	FullStateBackwardExpandCallback(const XoodooBase& xoodooBase, LINStateMILPTypeBySheet& solVars, vector<StateTypeBySheet>& solutions) : xoodooBase(xoodooBase)
	{
		this->solVars = solVars;
		this->pSolutions = &solutions;
	}

protected:
	void callback()
	{
		try {
			if (where == GRB_CB_MIPSOL)
			{
				// read solution
				StateTypeBySheet sol = { 0 };
				for (int x = 0; x < 4; x++)
					for (int y = 0; y < 3; y++)
						for (int z = 0; z < 32; z++)
							if (round(getSolution(solVars[x][y][z])) == 1)
								xoodooBase.SetBit(sol, x, y, z, true);
							else
								xoodooBase.SetBit(sol, x, y, z, false);


				// add lazy constr
				GRBLinExpr excludeCon = 0;
				for (int x = 0; x < 4; x++)
					for (int y = 0; y < 3; y++)
						for (int z = 0; z < 32; z++)
							if ((sol[x][y] >> z) & 0x1)
								excludeCon += (1 - solVars[x][y][z]);
							else
								excludeCon += solVars[x][y][z];
				addLazy(excludeCon >= 1);

				(*pSolutions).emplace_back(sol);

				// cout << "find one backward solution by callback: ";
				// xoodooBase.OutputState(cout, sol);
			}
		}
		catch (GRBException e) {
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Error during callback" << endl;
		}
	}
};

class HalfStateBackwardExpandCallback : public GRBCallback
{
public:
	const XoodooBase& xoodooBase;
	vector<HalfStateTypeBySheet>* pSolutions = NULL;
	ILINStateMILPTypeBySheet solVars;


	HalfStateBackwardExpandCallback(const XoodooBase& xoodooBase, ILINStateMILPTypeBySheet& solVars, vector<HalfStateTypeBySheet>& solutions) : xoodooBase(xoodooBase)
	{
		this->solVars = solVars;
		this->pSolutions = &solutions;
	}

protected:
	void callback()
	{
		try {
			if (where == GRB_CB_MIPSOL)
			{
				// read solution
				HalfStateTypeBySheet sol = { 0 };
				for (int x = 0; x < 4; x++)
					for (int y = 0; y < 3; y++)
						for (int z = 0; z < 16; z++)
							if (round(getSolution(solVars[x][y][z])) == 1)
								xoodooBase.SetBit(sol, x, y, z, true);
							else
								xoodooBase.SetBit(sol, x, y, z, false);


				// add lazy constr
				GRBLinExpr excludeCon = 0;
				for (int x = 0; x < 4; x++)
					for (int y = 0; y < 3; y++)
						for (int z = 0; z < 16; z++)
							if ((sol[x][y] >> z) & 0x1)
								excludeCon += (1 - solVars[x][y][z]);
							else
								excludeCon += solVars[x][y][z];
				addLazy(excludeCon >= 1);

				(*pSolutions).emplace_back(sol);

				// cout << "find one backward solution by callback: ";
				// xoodooBase.OutputState(cout, sol);
			}
		}
		catch (GRBException e) {
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Error during callback" << endl;
		}
	}
};



class FullStateForwardExpandCallback : public GRBCallback
{
public:
	const XoodooBase& xoodooBase;
	vector<StateTypeBySheet>* pSolutions = NULL;
	RDStateMILPTypeBySheet solVars;


	FullStateForwardExpandCallback(const XoodooBase& xoodooBase, RDStateMILPTypeBySheet& solVars, vector<StateTypeBySheet>& solutions) : xoodooBase(xoodooBase)
	{
		this->solVars = solVars;
		this->pSolutions = &solutions;
	}

protected:
	void callback()
	{
		try {
			if (where == GRB_CB_MIPSOL)
			{
				// read solution
				StateTypeBySheet sol = { 0 };
				for (int x = 0; x < 4; x++)
					for (int y = 0; y < 3; y++)
						for (int z = 0; z < 32; z++)
							if (round(getSolution(solVars[x][y][z])) == 1)
								xoodooBase.SetBit(sol, x, y, z, true);
							else
								xoodooBase.SetBit(sol, x, y, z, false);
							

				// add lazy constr
				GRBLinExpr excludeCon = 0;
				for (int x = 0; x < 4; x++)
					for (int y = 0; y < 3; y++)
						for (int z = 0; z < 32; z++)
							if ((sol[x][y] >> z) & 0x1)
								excludeCon += (1 - solVars[x][y][z]);
							else
								excludeCon += solVars[x][y][z];
				addLazy(excludeCon >= 1);

				(*pSolutions).emplace_back(sol);

				// cout << "find one forward solution by callback: ";
				// xoodooBase.OutputState(cout, sol);
			}
		}
		catch (GRBException e) {
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Error during callback" << endl;
		}
	}
};


class HalfStateForwardExpandCallback : public GRBCallback
{
public:
	const XoodooBase& xoodooBase;
	vector<HalfStateTypeBySheet>* pSolutions = NULL;
	IDStateMILPTypeBySheet solVars;


	HalfStateForwardExpandCallback(const XoodooBase& xoodooBase, IDStateMILPTypeBySheet& solVars, vector<HalfStateTypeBySheet>& solutions) : xoodooBase(xoodooBase)
	{
		this->solVars = solVars;
		this->pSolutions = &solutions;
	}

protected:
	void callback()
	{
		try {
			if (where == GRB_CB_MIPSOL)
			{
				// read solution
				HalfStateTypeBySheet sol = { 0 };
				for (int x = 0; x < 4; x++)
					for (int y = 0; y < 3; y++)
						for (int z = 0; z < 16; z++)
							if (round(getSolution(solVars[x][y][z])) == 1)
								xoodooBase.SetBit(sol, x, y, z, true);
							else
								xoodooBase.SetBit(sol, x, y, z, false);


				// add lazy constr
				GRBLinExpr excludeCon = 0;
				for (int x = 0; x < 4; x++)
					for (int y = 0; y < 3; y++)
						for (int z = 0; z < 16; z++)
							if ((sol[x][y] >> z) & 0x1)
								excludeCon += (1 - solVars[x][y][z]);
							else
								excludeCon += solVars[x][y][z];
				addLazy(excludeCon >= 1);

				(*pSolutions).emplace_back(sol);

				// cout << "find one forward solution by callback: ";
				// xoodooBase.OutputState(cout, sol);
			}
		}
		catch (GRBException e) {
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Error during callback" << endl;
		}
	}
};

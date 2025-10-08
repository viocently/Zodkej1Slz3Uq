#pragma once
#include<vector>
#include<set>
#include"ascon.h"
#include"differential.h"
#include"linear.h"
#include"gurobi_c++.h"


class FullStateBackwardExpandCallback : public GRBCallback
{
public:
	const AsconBase& asconBase;
	std::set<StateType>* pSolutions = NULL;
	LINStateMILPType solVars;


	FullStateBackwardExpandCallback(const AsconBase& asconBase, LINStateMILPType& solVars, std::set<StateType>& solutions) : asconBase(asconBase)
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
				StateType sol = { 0 };
				for (int x = 0; x < 5; x++)
					for (int y = 0; y < 64; y++)
						if (round(getSolution(solVars[x][y])) == 1)
							asconBase.SetBit(sol, x, y, true);
						else
							asconBase.SetBit(sol, x, y, false);


				// add lazy constr
				GRBLinExpr excludeCon = 0;
				for (int x = 0; x < 5; x++)
					for (int y = 0; y < 64; y++)
						if ((sol[x] >> y) & 0x1)
							excludeCon += (1 - solVars[x][y]);
						else
							excludeCon += solVars[x][y];
				addLazy(excludeCon >= 1);

				(*pSolutions).emplace(sol);

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


using FullStateForwardExpandCallback = FullStateBackwardExpandCallback;

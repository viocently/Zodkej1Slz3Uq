Each folder of a cipher name contains the minizinc model file for RDL or IDL distinguisher search on the cipher. For example, the folder 'alzette' contains the models of SDL and RDL distinguishers search on Alzette. Here is an example for how to solve the models using minizinc command line in ubuntu: 

`nohup ~/MiniZincIDE-2.8.5-bundle-linux-x86_64/bin/minizinc   --solver com.google.ortools.sat -O2 -v -a  -f -p 64 asconRDL.mzn &`

Remember to use Or-Tools with multi-threads in Minizinc, because this will bring a significant speed up for solving the models.

The folders 'asconRDLComp' and 'xoodooRDLAndIDLComp' contains the C++ files for computing the RDL (or IDL) correlation of 3-round or 4-round EM, when the input rotational (or internal) difference and the linear mask on the output rotational (or internal) difference of EM are known. Note that to compile these files, you need to install gurobi and open the option -fopenmp in g++, because we will use openmp to compute the correlation with multi-threads. Here is an example for how to compile these files in ubuntu: 
`g++ asconRDL.cpp -o asconRDL -std=c++17 -O2 -lm -lpthread -fopenmp -I/$GUROBI_HOME/include/ -L/$GUROBI_HOME/lib -lgurobi_c++ -lgurobi91 -lm`

include "globals.mzn";
int: n = 32;
int: rot = 20;
int: bits = 7;
int: precision = pow(2, bits);

array[0..n-1] of 0..1: rc = array1d(0..n-1, [0,1,0,0,0,1,1,0,1,0,0,0,1,0,1,0,1,0,0,0,0,1,1,1,1,1,1,0,1,1,0,1]);
%array[0..n-1] of 0..1: rc = array1d(0..n-1, [0,1,1,0,1,0,1,0,0,1,0,1,1,0,1,1,0,0,1,0,1,1,0,1,0,0,0,1,1,1,0,0]);

predicate difModularAddition(array[int] of var 0..1: x, array[int] of var 0..1: y, array[int] of var 0..1: z, int: m, var int: pLog2) = 
    let {
      array[0..m-2] of var 0..1: neq;
      
      constraint forall(i in 0..m-2)
      (
        if x[i] == y[i] /\ y[i] == z[i] then
          neq[i] = 0
        else
          neq[i] = 1
        endif
      );
    }
    in
      x[0]+y[0]+z[0] in {0,2} 
      /\
      (
        forall(i in 1..m-1)
        (
          if neq[i-1] == 0 then 
            x[i]+y[i]+z[i]+x[i-1] in {0,2,4} 
          endif
        ) 
      ) 
      /\
      pLog2 = sum(i in 0..m-2)(neq[i]);
      
predicate rDifModularAddition(array[int] of var 0..1: x, array[int] of var 0..1: y, array[int] of var 0..1: z, var int: pLog2) =
    let{
      array[0..rot-1] of var 0..1: xLow;
      array[0..n-rot-1] of var 0..1: xHigh;
      array[0..rot-1] of var 0..1: yLow;
      array[0..n-rot-1] of var 0..1: yHigh;
      array[0..rot-1] of var 0..1: zLow;
      array[0..n-rot-1] of var 0..1: zHigh;
      
      constraint forall(i in 0..rot-1)(xLow[i] = x[i]);
      constraint forall(i in 0..rot-1)(yLow[i] = y[i]);
      constraint forall(i in 0..rot-1)(zLow[i] = z[i]);
      
      constraint forall(i in 0..n-rot-1)(xHigh[i] = x[rot+i]);
      constraint forall(i in 0..n-rot-1)(yHigh[i] = y[rot+i]);
      constraint forall(i in 0..n-rot-1)(zHigh[i] = z[rot+i]);
      
      var int: pLog2Low;
      var int: pLog2High;
    }
    in
      difModularAddition(xLow, yLow, zLow, rot, pLog2Low)
      /\
      difModularAddition(xHigh,yHigh, zHigh, n-rot, pLog2High)
      /\
      pLog2 = pLog2Low + pLog2High + 2;
    


predicate linModularAddition(array[int] of var 0..1: x, array[int] of var 0..1: y, array[int] of var 0..1: z, var int: cor2Log2) =
    let {
      array[1..n] of var 0..1: c;
      
      constraint c[n] = 0;
    }
    in
      (
        forall(i in 1..n-1)
        (
          if c[i+1] == 0 then
            x[i] = z[i] /\ y[i] = z[i] /\ c[i] = z[i]
          else
            x[i]+y[i]+c[i]+z[i] in {1,3}
          endif
        )
      )
      /\
      (
        if c[1] == 0 then
          x[0] = z[0] /\ y[0] = z[0]
        endif
      )
      /\
      cor2Log2 = sum(i in 1..n)(c[i]);

% this function is specified for rotational difference
predicate DLModularAdditionV2(array[int] of var 0..1: x, array[int] of var 0..1: y, array[int] of var 0..1: z, int: m, var int: corLog2) =
    let {
      array[0..m-1] of var -precision..precision: negParts;
      array[0..m-1] of var -precision..precision: constParts;
      array[0..m-1] of var -precision..precision: negHelper;
      array[0..m-1] of var -precision..precision: constHelper;
      array[0..m-1] of var 1..bits: corsLog2;

      constraint negParts[m-1] = 0;
      constraint constParts[m-1] = precision;
    }
    in
      forall(i in 0..m-1)
      (
        if z[i] == 1 then
          (negHelper[i] = constParts[i]) /\ (constHelper[i] = negParts[i])
        else
          (negHelper[i] = negParts[i]) /\ (constHelper[i] = constParts[i])
        endif
      )
      /\
      forall(i in 1..m-1)
      (
        if x[i-1] == y[i-1] then
          (2*negParts[i-1] = negHelper[i]) /\ (corsLog2[i] = bits)
          /\
          if x[i-1] == 1 then
            2*constParts[i-1] = 2*constHelper[i] - negHelper[i]
          else
            2*constParts[i-1] = 2*constHelper[i] + negHelper[i]
          endif
        else
          (constParts[i-1] = precision) /\ (negParts[i-1] = 0)
          /\
          (abs(constHelper[i]) > pow(2, corsLog2[i]-1)) /\ (abs(constHelper[i]) <= pow(2, corsLog2[i]))
        endif
      )
      /\
      (abs(constHelper[0]) > pow(2, corsLog2[0]-1))
      /\
      (abs(constHelper[0]) <= pow(2, corsLog2[0]))
      /\
      (corLog2 = sum(i in 0..m-1)(bits - corsLog2[i]));


predicate RDLModularAddition(array[int] of var 0..1: x, array[int] of var 0..1: y, array[int] of var 0..1: z, var int: corLog2) =
    let{
      array[0..rot-1] of var 0..1: xLow;
      array[0..n-rot-1] of var 0..1: xHigh;
      array[0..rot-1] of var 0..1: yLow;
      array[0..n-rot-1] of var 0..1: yHigh;
      array[0..rot-1] of var 0..1: zLow;
      array[0..n-rot-1] of var 0..1: zHigh;
      
      constraint forall(i in 0..rot-1)(xLow[i] = x[i]);
      constraint forall(i in 0..rot-1)(yLow[i] = y[i]);
      constraint forall(i in 0..rot-1)(zLow[i] = z[i]);
      
      constraint forall(i in 0..n-rot-1)(xHigh[i] = x[rot+i]);
      constraint forall(i in 0..n-rot-1)(yHigh[i] = y[rot+i]);
      constraint forall(i in 0..n-rot-1)(zHigh[i] = z[rot+i]);
      
      var int: corLog2Low;
      var int: corLog2High;
    }
    in
      DLModularAdditionV2(xLow, yLow, zLow, rot, corLog2Low)
      /\
      DLModularAdditionV2(xHigh, yHigh, zHigh, n-rot, corLog2High)
      /\
      corLog2 = corLog2Low + corLog2High;

      
predicate difRound(array[int] of var 0..1: leftInput, array[int] of var 0..1: rightInput, array[int] of var 0..1: leftOutput, array[int] of var 0..1: rightOutput, int: s, int: t, var int: pLog2) = 
    let {
      array[0..n-1] of var 0..1: modularAdditionRightInput;
      array[0..n-1] of var 0..1: xorLeftInput;
      array[0..n-1] of var 0..1: leftOutputBeforeRC;
      
      constraint forall(i in 0..n-1)(leftOutputBeforeRC[i] = xorLeftInput[(i-t+n) mod n]);
      constraint forall(i in 0..n-1)(rightInput[i] = modularAdditionRightInput[(i-s+n) mod n]);
    } 
    in
      rDifModularAddition(leftInput, modularAdditionRightInput, leftOutputBeforeRC, pLog2)
      /\
      forall(i in 0..n-1)(rightOutput[i]+xorLeftInput[i]+rightInput[i] in {0,2})
      /\
      forall(i in 0..n-1)(leftOutputBeforeRC[i]+leftOutput[i]+rc[i]+rc[(i-rot+n) mod n] in {0,2,4});

predicate linRound(array[int] of var 0..1: leftInput, array[int] of var 0..1: rightInput, array[int] of var 0..1: leftOutput, array[int] of var 0..1: rightOutput, int: s, int: t, var int: corLog2) =
    let {
      array[0..n-1] of var 0..1: modularAdditionRightInput;
      array[0..n-1] of var 0..1: modularAdditionOutput;
      
      constraint forall(i in 0..n-1)(rightInput[i]+rightOutput[i]+modularAdditionRightInput[(i-s+n) mod n] in {0,2});
      constraint forall(i in 0..n-1)(leftOutput[i]+modularAdditionOutput[i]+rightOutput[(i-t+n) mod n] in {0,2});
    }
    in
      linModularAddition(leftInput, modularAdditionRightInput, modularAdditionOutput, corLog2);

predicate DLRound(array[int] of var 0..1: leftInput, array[int] of var 0..1: rightInput, array[int] of var 0..1: leftOutput, array[int] of var 0..1: rightOutput, int: s, int: t, var int: cor) =
    let {
      array[0..n-1] of var 0..1: modularAdditionOutput;
      array[0..n-1] of var 0..1: modularAdditionRightInput;
      
      constraint forall(i in 0..n-1)(rightInput[i] = modularAdditionRightInput[(i-s+n) mod n]);
      constraint forall(i in 0..n-1)(leftOutput[i]+modularAdditionOutput[i]+rightOutput[(i-t+n) mod n] in {0,2});
    }          
    in
      RDLModularAddition(leftInput, modularAdditionRightInput, modularAdditionOutput, cor);
      

int: r0 = 2;
int: r1 = 5;
array[0..3, 0..1] of int: rots = array2d(0..3,0..1,[31,24,17,17,0,31,24,16]);
array[0..r0, 0..1, 0..n-1] of var 0..1: rDifVars;
array[0..r0-1] of var int: rDifProbs;
array[0..r1, 0..1, 0..n-1] of var 0..1: linVars;
array[0..r1-1] of var int: linCors;
array[0..r1, 0..1, 0..n-1] of var 0..1: rLinVars;
array[0..r1-1] of var int: rLinCors;

var int: RDLCTCor;

constraint forall(i in 0..r0-1)(difRound(rDifVars[i,0,..], rDifVars[i,1,..], rDifVars[i+1,0,..], rDifVars[i+1,1,..], rots[i mod 4, 0], rots[i mod 4, 1], rDifProbs[i]));
constraint DLRound(rDifVars[r0,0,..], rDifVars[r0,1,..], linVars[0,0,..], linVars[0,1,..], rots[r0 mod 4,0], rots[r0 mod 4, 1], RDLCTCor);
constraint forall(i in 0..r1-1)(linRound(linVars[i,0,..], linVars[i,1,..], linVars[i+1,0,..], linVars[i+1,1,..], rots[(i+r0+1) mod 4, 0], rots[(r0+1+i) mod 4, 1], linCors[i]));
constraint forall(i in 0..r1-1)(linRound(rLinVars[i,0,..], rLinVars[i,1,..], rLinVars[i+1,0,..], rLinVars[i+1,1,..], rots[(i+r0+1) mod 4, 0], rots[(r0+1+i) mod 4, 1], rLinCors[i]));

constraint forall(i in 0..r1, j in 0..n-1)((linVars[i,0,j] = rLinVars[i,0,(j-rot+n) mod n]) /\ (linVars[i,1,j] = rLinVars[i,1,(j-rot+n) mod n]));


var int: pU = sum(i in 0..r0-1)(rDifProbs[i]);
var int: cL = sum(i in 0..r1-1)(linCors[i]);
var int: rCL = sum(i in 0..r1-1)(rLinCors[i]);




var int: obj = pU + cL + rCL + RDLCTCor;



constraint sum(i in 0..n-1)(linVars[r1,0,i] + linVars[r1,1,i]) > 0;
constraint if rot == 0 then sum(i in 0..n-1)(rDifVars[0,0,i] + rDifVars[0,1,i]) > 0 endif;

solve minimize obj;

output ["\(rDifVars[i,0,..]) -- \(rDifVars[i,1,..]) -- D \(rDifProbs[i]) \n" | i in 0..r0-1];
output ["\(rDifVars[r0,0,..]) -- \(rDifVars[r0,1,..]) \n"];
output ["\(linVars[i,0,..]) -- \(linVars[i,1,..]) -- L \(linCors[i]),\(rLinCors[i]) \n" | i in 0..r1-1];
output ["\(linVars[r1,0,..]) -- \(linVars[r1,1,..]) \n"];
output ["\(pU) -- \(RDLCTCor) -- \(cL),\(rCL) \n"];

%[0..n-1] of var 0..1: leftInput;
%array[0..n-1] of var 0..1: rightInput;
%array[0..n-1] of var 0..1: leftOutput;
%array[0..n-1] of var 0..1: rightOutput;
%var int: pLog2;
%var int: corLog2;
%constraint linRound(leftInput, rightInput, leftOutput, rightOutput, 31, 24, corLog2);

      
%array[0..n-1] of var 0..1: xDif = array1d(0..n-1, [1,0,0,0,0,0,1,1,1,1,1,0,1,1,0,0,0,0,1,0,0,0,1,1,0,1,0,1,1,1,0,1]);
%array[0..n-1] of var 0..1: yDif = array1d(0..n-1, [0,0,0,0,0,0,1,1,0,1,1,1,0,1,1,1,1,1,1,0,0,1,0,0,0,0,1,0,0,0,0,1]);
%array[0..n-1] of var 0..1: zMask = array1d(0..n-1, [0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]);
%var int: DLCor;
%constraint DLModularAddition(xDif, yDif, zMask, DLCor);

%array[0..n-1] of var 0..1: xMask = array1d(0..n-1, [0,1,1,1,1]);
%array[0..n-1] of var 0..1: yMask = array1d(0..n-1, [0,0,0,1,1]);
%array[0..n-1] of var 0..1: zMask = array1d(0..n-1, [0,0,1,0,1]);
%var int: linCorLog2;
%constraint linModularAddition(xMask, yMask, zMask, linCorLog2);

%array[0..n-1] of var 0..1: xDif = array1d(0..n-1, [0,0,1,1,1]);
%array[0..n-1] of var 0..1: yDif = array1d(0..n-1, [0,1,1,0,0]);
%array[0..n-1] of var 0..1: zDif = array1d(0..n-1, [0,1,1,0,1]);
%var int: difProbLog2;
%constraint difModularAddition(xDif, yDif, zDif, difProbLog2);



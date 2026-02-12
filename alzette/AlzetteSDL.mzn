include "globals.mzn";
int: n = 32;
int: bits = 12;
int: precision = pow(2, bits);

predicate difModularAddition(array[int] of var 0..1: x, array[int] of var 0..1: y, array[int] of var 0..1: z, var int: pLog2) = 
    let {
      array[0..n-2] of var 0..1: neq;
      
      constraint forall(i in 0..n-2)
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
        forall(i in 1..n-1)
        (
          if neq[i-1] == 0 then 
            x[i]+y[i]+z[i]+x[i-1] in {0,2,4} 
          endif
        ) 
      ) 
      /\
      pLog2 = sum(i in 0..n-2)(neq[i]);


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
      
predicate DLModularAddition(array[int] of var 0..1: x, array[int] of var 0..1: y, array[int] of var 0..1: z, var int: cor) = 
    let {
      array[1..n] of var -precision..precision: negParts;
      array[1..n] of var -precision..precision: constParts;
      
      constraint negParts[n] = 0;
      constraint constParts[n] = precision;
    }
    in
      (
        forall(i in 1..n-1)
        (
          if z[i] == 1 then
            if x[i-1] == y[i-1] then
              2*negParts[i]=constParts[i+1]
              /\ 
              (
                if x[i-1] == 1 then
                  2*constParts[i] + constParts[i+1] = 2*negParts[i+1]
                else
                  2*constParts[i] = 2*negParts[i+1] + constParts[i+1]  
                endif
              )    
            else
              constParts[i] = negParts[i+1] /\ negParts[i] = 0
            endif
          else
            if x[i-1] == y[i-1] then
              2*negParts[i] = negParts[i+1]
              /\
              (
                if x[i-1] == 1 then
                  2*constParts[i] + negParts[i+1] = 2*constParts[i+1]
                else
                  2*constParts[i] = 2*constParts[i+1] + negParts[i+1]
                endif
              )
            else
              negParts[i] = 0 /\ constParts[i] = constParts[i+1]
            endif
          endif
        )
      )
      /\
      cor = negParts[1] + constParts[1];
      
      
predicate difRound(array[int] of var 0..1: leftInput, array[int] of var 0..1: rightInput, array[int] of var 0..1: leftOutput, array[int] of var 0..1: rightOutput, int: s, int: t, var int: pLog2) = 
    let {
      array[0..n-1] of var 0..1: modularAdditionRightInput;
      array[0..n-1] of var 0..1: xorLeftInput;
      
      constraint forall(i in 0..n-1)(leftOutput[i] = xorLeftInput[(i-t+n) mod n]);
      constraint forall(i in 0..n-1)(rightInput[i] = modularAdditionRightInput[(i-s+n) mod n]);
    } 
    in
      difModularAddition(leftInput, modularAdditionRightInput, leftOutput, pLog2)
      /\
      forall(i in 0..n-1)(rightOutput[i]+xorLeftInput[i]+rightInput[i] in {0,2});

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
      DLModularAddition(leftInput, modularAdditionRightInput, modularAdditionOutput, cor);
      

int: r0 = 6;
int: r1 = 6;
array[0..3, 0..1] of int: rots = array2d(0..3,0..1,[31,24,17,17,0,31,24,16]);
array[0..r0+r1+1,0..1,0..n-1] of var 0..1: vars;
array[0..r0+r1] of var int: vals;

constraint forall(i in 0..r0-1)(difRound(vars[i,0,..], vars[i,1,..], vars[i+1,0,..], vars[i+1,1,..], rots[i mod 4, 0], rots[i mod 4, 1], vals[i]));
constraint DLRound(vars[r0,0,..], vars[r0,1,..], vars[r0+1,0,..], vars[r0+1,1,..], rots[r0 mod 4,0], rots[r0 mod 4, 1], vals[r0]);
constraint forall(i in r0+1..r0+r1)(linRound(vars[i,0,..], vars[i,1,..], vars[i+1,0,..], vars[i+1,1,..], rots[i mod 4, 0], rots[i mod 4, 1], vals[i]));

var int: pU = sum(i in 0..r0-1)(vals[i]);
var int: cL = sum(i in r0+1..r0+r1)(2*vals[i]);

%int: wu = precision;
%int: wm = 1;
%int: wl = precision;

%var int: obj = wu*pU + wl*cL - wm*vals[r0];

%constraint 4*vals[r0] > precision;
var 1..bits: cM;
constraint abs(vals[r0]) > pow(2,cM-1);
constraint abs(vals[r0]) <= pow(2,cM);



var int: obj = pU + cL + bits - cM;



constraint (sum(i in 0..n-1)(vars[0,0,i]+vars[0,1,i])) > 0;
constraint (sum(i in 0..n-1)(vars[r0+r1+1,0,i]+vars[r0+r1+1,1,i])) > 0;

solve minimize obj;

output ["\(vars[i,0,..]) -- \(vars[i,1,..]) -- D \(vals[i]) \n" | i in 0..r0-1];
output ["\(vars[r0,0,..]) -- \(vars[r0,1,..]) \n"];
output ["\(vars[i,0,..]) -- \(vars[i,1,..]) -- L \(vals[i]) \n" | i in r0+1..r0+r1];
output ["\(vars[r0+r1+1,0,..]) -- \(vars[r0+r1+1,1,..]) \n"];
output ["\(pU) -- \(cM) -- \(cL) \n"];

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




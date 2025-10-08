include "globals.mzn";
int: n = 32;
int: bits = 20;
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
      array[0..n-1] of var 0..1: modularAdditionUpInput;
      array[0..n-1] of var 0..1: xorUpInput;
      
      constraint forall(i in 0..n-1)(modularAdditionUpInput[i] = leftInput[(i+s) mod n]);
      constraint forall(i in 0..n-1)(xorUpInput[i] = rightInput[(i-t+n) mod n]);
    } 
    in
      difModularAddition(modularAdditionUpInput, rightInput, leftOutput, pLog2)
      /\
      forall(i in 0..n-1)(leftOutput[i]+xorUpInput[i] + rightOutput[i] in {0,2});
      
predicate linRound(array[int] of var 0..1: leftInput, array[int] of var 0..1: rightInput, array[int] of var 0..1: leftOutput, array[int] of var 0..1: rightOutput, int: s, int: t, var int: corLog2) =
    let {
      array[0..n-1] of var 0..1: modularAdditionUpInput;
      array[0..n-1] of var 0..1: modularAdditionRightInput;
      array[0..n-1] of var 0..1: modularAdditionOutput;
      
      constraint forall(i in 0..n-1)(modularAdditionUpInput[i] = leftInput[(i+s) mod n]);
      constraint forall(i in 0..n-1)(modularAdditionOutput[i] + rightOutput[i] + leftOutput[i] in {0,2});
      constraint forall(i in 0..n-1)(rightInput[i] + modularAdditionRightInput[i] + rightOutput[(i+t) mod n] in {0,2});
    }
    in
      linModularAddition(modularAdditionUpInput, modularAdditionRightInput, modularAdditionOutput, corLog2);
      
predicate DLRound(array[int] of var 0..1: leftInput, array[int] of var 0..1: rightInput, array[int] of var 0..1: leftOutput, array[int] of var 0..1: rightOutput, int: s, int: t, var int: cor) =
    let {
      array[0..n-1] of var 0..1: modularAdditionUpInput;
      array[0..n-1] of var 0..1: modularAdditionOutput;
      
      constraint forall(i in 0..n-1)(modularAdditionUpInput[i] = leftInput[(i+s) mod n]);
      constraint forall(i in 0..n-1)(modularAdditionOutput[i] + rightOutput[i] + leftOutput[i] in {0,2});
    }          
    in
      DLModularAddition(modularAdditionUpInput, rightInput, modularAdditionOutput, cor);
      
int: r0 = 6;
int: r1 = 5;
int: alpha = 8;
int: beta = 3;
array[0..r0+r1+1,0..1,0..n-1] of var 0..1: vars;
array[0..r0+r1] of var int: vals;

constraint forall(i in 0..r0-1)(difRound(vars[i,0,..], vars[i,1,..], vars[i+1,0,..], vars[i+1,1,..], alpha, beta, vals[i]));
constraint DLRound(vars[r0,0,..], vars[r0,1,..], vars[r0+1,0,..], vars[r0+1,1,..], alpha, beta, vals[r0]);
constraint forall(i in r0+1..r0+r1)(linRound(vars[i,0,..], vars[i,1,..], vars[i+1,0,..], vars[i+1,1,..], alpha, beta, vals[i]));

var int: pU = sum(i in 0..r0-1)(vals[i]);
var int: cL = sum(i in r0+1..r0+r1)(2*vals[i]);

var 0..bits: cM;
constraint abs(vals[r0]) >= pow(2,cM);
constraint abs(vals[r0]) < pow(2,cM+1);



var int: obj = pU + cL + bits - cM;



constraint (sum(i in 0..n-1)(vars[0,0,i]+vars[0,1,i])) > 0;
constraint (sum(i in 0..n-1)(vars[r0+r1+1,0,i]+vars[r0+r1+1,1,i])) > 0;

solve minimize obj;

output ["\(vars[i,0,..]) -- \(vars[i,1,..]) -- D \(vals[i]) \n" | i in 0..r0-1];
output ["\(vars[r0,0,..]) -- \(vars[r0,1,..]) \n"];
output ["\(vars[i,0,..]) -- \(vars[i,1,..]) -- L \(vals[i]) \n" | i in r0+1..r0+r1];
output ["\(vars[r0+r1+1,0,..]) -- \(vars[r0+r1+1,1,..]) \n"];
output ["\(pU) -- \(cM) -- \(cL) \n"];
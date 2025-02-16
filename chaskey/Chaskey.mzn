include "globals.mzn";
int: n = 32;
int: bits = 20;
int: precision = pow(2, bits);



predicate DLModularAddition(array[int] of var 0..1: x, array[int] of var 0..1: y, array[int] of var 0..1: z, var int: corLog2) =
    let {
      array[0..n-1] of var -precision..precision: negParts;
      array[0..n-1] of var -precision..precision: constParts;
      array[1..n-1] of var -precision..precision: negHelper;
      array[1..n-1] of var -precision..precision: constHelper;

      constraint negParts[n-1] = 0;
      constraint constParts[n-1] = precision;
      
      var 1..bits: tmp;
    }
    in
      forall(i in 1..n-1)
      (
        if z[i] == 1 then
          negHelper[i] = constParts[i] /\ constHelper[i] = negParts[i]
        else
          negHelper[i] = negParts[i] /\ constHelper[i] = constParts[i]
        endif
      )
      /\
      forall(i in 1..n-1)
      (
        if x[i-1] == y[i-1] then
          2*negParts[i-1] = negHelper[i] 
          /\
          if x[i-1] == 1 then
            2*constParts[i-1] = 2*constHelper[i] - negHelper[i]
          else
            2*constParts[i-1] = 2*constHelper[i] + negHelper[i]
          endif
        else
          constParts[i-1] = constHelper[i] /\ negParts[i-1] = 0
        endif
      )
      /\
      abs(constParts[0] + negParts[0]) > pow(2, tmp-1)
      /\
      abs(constParts[0] + negParts[0]) <= pow(2, tmp)
      /\
      corLog2 = bits - tmp;

predicate difModularAddition(array[int] of var bool: x, array[int] of var bool: y, array[int] of var bool: z, var int: pLog2) = 
    let {
      array[0..n-2] of var bool: neq;
      
      constraint forall(i in 0..n-2)
      (
        neq[i] = (not (x[i] == y[i] /\ y[i] == z[i]))
      );
    }
    in
      (z[0] = x[0] xor y[0])
      /\
      (
        forall(i in 1..n-1)
        (
          if not neq[i-1] then 
            x[i-1] = (x[i] xor y[i] xor z[i]) 
          endif
        ) 
      ) 
      /\
      pLog2 = sum(i in 0..n-2)(bool2int(neq[i]));
      
predicate linModularAddition(array[int] of var bool: x, array[int] of var bool: y, array[int] of var bool: z, var int: cor2Log2) =
    let {
      array[1..n] of var bool: c;
      
      constraint c[n] = false;
    }
    in
      (
        forall(i in 1..n-1)
        (
          if not c[i+1] then
            x[i] = z[i] /\ y[i] = z[i] /\ c[i] = z[i]
          else
            (x[i] xor y[i] xor c[i] xor z[i]) = 1
          endif
        )
      )
      /\
      (
        if not c[1] then
          x[0] = z[0] /\ y[0] = z[0]
        endif
      )
      /\
      cor2Log2 = sum(i in 1..n)(bool2int(c[i]));
      
      
predicate difRoundComponent0(array[0..3, 0..n-1] of var bool: inputDif, array[0..3, 0..n-1] of var bool: outputDif, var int: pLog2, int: s, int: t) =
    let {
      var int: upModularAdditionPLog2;
      var int: downModularAdditionPLog2;
    } 
    in
      forall(i in 0..n-1)(inputDif[0,i] == outputDif[0, (i+s) mod n])
      /\
      difModularAddition(inputDif[0,..], inputDif[1,..], outputDif[1,..], upModularAdditionPLog2)
      /\
      difModularAddition(inputDif[2,..], inputDif[3,..], outputDif[2,..], downModularAdditionPLog2)
      /\
      forall(i in 0..n-1)(inputDif[3,i] == outputDif[3, (i+t) mod n])
      /\
      pLog2 = upModularAdditionPLog2 + downModularAdditionPLog2;

predicate difRoundComponent1(array[0..3, 0..n-1] of var bool: inputDif, array[0..3, 0..n-1] of var bool: outputDif) = 
    forall(i in 0..n-1)(outputDif[0,i] = inputDif[0,i] xor inputDif[1,i])
    /\
    forall(i in 0..n-1)(inputDif[1,i] == outputDif[2,i])
    /\
    forall(i in 0..n-1)(inputDif[2,i] == outputDif[1,(i+16) mod n])
    /\
    forall(i in 0..n-1)(outputDif[3,i] = inputDif[2,i] xor inputDif[3,i]);
    
predicate linRoundComponent0(array[0..3, 0..n-1] of var bool: inputDif, array[0..3, 0..n-1] of var bool: outputDif, var int: corLog2, int: s, int: t) = 
    let {
      var int: upModularAdditionCorLog2;
      var int: downModularAdditionCorLog2;
      array[0..n-1] of var bool: upModularAdditionUpInput;
      array[0..n-1] of var bool: downModularAdditionDownInput;
    }
    in
      forall(i in 0..n-1)(upModularAdditionUpInput[i] = inputDif[0,i] xor outputDif[0,(i+s) mod n])
      /\
      linModularAddition(inputDif[1,..], upModularAdditionUpInput, outputDif[1,..], upModularAdditionCorLog2)
      /\
      linModularAddition(inputDif[2,..], downModularAdditionDownInput, outputDif[2,..], downModularAdditionCorLog2)
      /\
      forall(i in 0..n-1)(downModularAdditionDownInput[i] = inputDif[3,i] xor outputDif[3,(i+t) mod n])
      /\
      corLog2 = upModularAdditionCorLog2 + downModularAdditionCorLog2;

predicate linRoundComponent1(array[0..3, 0..n-1] of var bool: inputDif, array[0..3, 0..n-1] of var bool: outputDif) = 
    forall(i in 0..n-1)
    (
      outputDif[0,i] == inputDif[0,i]
      /\
      (outputDif[2,i] = inputDif[0,i] xor inputDif[1,i])
      /\
      (outputDif[1,(i+16) mod n] = inputDif[2,i] xor inputDif[3,i])
      /\
      inputDif[3,i] == outputDif[3,i]
    );

predicate DLRoundComponent0(array[0..3, 0..n-1] of var bool: inputDif, array[0..3, 0..n-1] of var bool: outputMask, var int: upModularAdditionCor, var int: downModularAdditionCor, int: s, int: t) = 
    DLModularAddition(inputDif[0,..], inputDif[1,..], outputMask[1,..], upModularAdditionCor)
    /\
    DLModularAddition(inputDif[2,..], inputDif[3,..], outputMask[2,..], downModularAdditionCor);



predicate difRound(array[0..3, 0..n-1] of var bool: inputDif, array[0..3, 0..n-1] of var bool: outputDif, var int: pLog2) = 
    let {
      array[0..3, 0..n-1] of var bool: afterFirstQuarter;
      array[0..3, 0..n-1] of var bool: afterSecondQuarter;
      array[0..3, 0..n-1] of var bool: afterThirdQuarter;
      var int: leftPLog2;
      var int: rightPLog2;
      
    }
    in
      difRoundComponent0(inputDif, afterFirstQuarter, leftPLog2, 8, 5)
      /\
      difRoundComponent1(afterFirstQuarter, afterSecondQuarter)
      /\
      difRoundComponent0(afterSecondQuarter, afterThirdQuarter, rightPLog2, 13, 7)
      /\
      difRoundComponent1(afterThirdQuarter, outputDif)
      /\
      pLog2 = leftPLog2 + rightPLog2;
      
predicate linRound(array[0..3, 0..n-1] of var bool: inputMask, array[0..3, 0..n-1] of var bool: outputMask, var int: corLog2) = 
    let {
      array[0..3, 0..n-1] of var bool: afterFirstQuarter;
      array[0..3, 0..n-1] of var bool: afterSecondQuarter;
      array[0..3, 0..n-1] of var bool: afterThirdQuarter;
      var int: leftCorLog2;
      var int: rightCorLog2;
    }
    in
      linRoundComponent0(inputMask, afterFirstQuarter, leftCorLog2, 8, 5)
      /\
      linRoundComponent1(afterFirstQuarter, afterSecondQuarter)
      /\
      linRoundComponent0(afterSecondQuarter, afterThirdQuarter, rightCorLog2, 13, 7)
      /\
      linRoundComponent1(afterThirdQuarter, outputMask)
      /\
      corLog2 = leftCorLog2 + rightCorLog2;
      
predicate leftHalfDifRound(array[0..3, 0..n-1] of var bool: inputDif, array[0..3, 0..n-1] of var bool: outputDif, var int: pLog2) = 
    let {
      array[0..3, 0..n-1] of var bool: afterFirstQuarter;
    }
    in
      difRoundComponent0(inputDif, afterFirstQuarter, pLog2, 8, 5)
      /\
      difRoundComponent1(afterFirstQuarter, outputDif);

predicate rightHalfLinRound(array[0..3, 0..n-1] of var bool: inputMask, array[0..3, 0..n-1] of var bool: outputMask, var int: linCorLog2) =
    let {
      array[0..3, 0..n-1] of var bool: afterSecondQuarter;
      array[0..3, 0..n-1] of var bool: afterThirdQuarter;
    }
    in
      linRoundComponent1(inputMask, afterSecondQuarter)
      /\
      linRoundComponent0(afterSecondQuarter, afterThirdQuarter, linCorLog2, 13, 7)
      /\
      linRoundComponent1(afterThirdQuarter, outputMask);
   


int: r0 = 3;
int: r1 = 3;

array[0..r0, 0..127] of var bool: difVars;
array[0..r0-1] of var int: pLog2;



array[0..r1, 0..127] of var bool: linVars;
array[0..r1-1] of var int: corLog2;


var 0..bits: DLRoundDLCor0Log2;
var 0..bits: DLRoundDLCor1Log2;
var 0..2*bits: DLRoundDifProbLog2;

constraint forall(i in 0..r0-1)(difRound(array2d(0..3,0..n-1, difVars[i,..]), array2d(0..3,0..n-1, difVars[i+1,..]), pLog2[i]));
constraint forall(i in 0..r1-1)(linRound(array2d(0..3,0..n-1, linVars[i,..]), array2d(0..3,0..n-1, linVars[i+1,..]), corLog2[i]));

array[0..3, 0..n-1] of var bool: afterSecondQuarter;
array[0..3, 0..n-1] of var bool: afterThirdQuarter;
constraint leftHalfDifRound(array2d(0..3,0..n-1, difVars[r0,..]), afterSecondQuarter, DLRoundDifProbLog2);
constraint DLRoundComponent0(afterSecondQuarter, afterThirdQuarter, DLRoundDLCor0Log2, DLRoundDLCor1Log2, 13,7);
constraint linRoundComponent1(afterThirdQuarter, array2d(0..3,0..n-1, linVars[0,..]));

constraint sum(i in 0..127)(difVars[0,i]) > 0;
constraint sum(i in 0..127)(linVars[r1,i]) > 0;




var int: pU = sum(i in 0..r0-1)(pLog2[i]);
var int: cL = sum(i in 0..r1-1)(2*corLog2[i]);
var int: cM = DLRoundDLCor0Log2 + DLRoundDLCor1Log2;

var 0..100: obj = pU + cL + DLRoundDifProbLog2 + cM;

ann: search_ann;
search_ann = seq_search([
    int_search([linVars[r1,z] | z in 0..127], indomain_min, complete)
    ]);

solve minimize obj;

output ["\(bool2int(difVars[i,..])) -- D \(pLog2[i]) \n " | i in 0..r0-1];
output ["\(bool2int(difVars[r0,..])) -- D \(DLRoundDifProbLog2) \n"];
output ["\(bool2int(afterSecondQuarter)) -- DL \(DLRoundDLCor0Log2) , \(DLRoundDLCor1Log2) \n"];
output ["\(bool2int(afterThirdQuarter))\n"];
output ["\(bool2int(linVars[i,..])) -- L \(corLog2[i]) \n " | i in 0..r1-1];
output ["\(bool2int(linVars[r1,..])) \n"];
output ["\(pU),\(DLRoundDifProbLog2) -- \(DLRoundDLCor0Log2),\(DLRoundDLCor1Log2) -- \(cL) \n"];

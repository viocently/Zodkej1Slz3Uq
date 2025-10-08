include "globals.mzn";
int: n = 64;
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
      

predicate DLModularAddition(array[int] of var 0..1: x, array[int] of var 0..1: y, array[int] of var 0..1: z, var int: corLog2) =
    let {
      array[0..n-1] of var -precision..precision: negParts;
      array[0..n-1] of var -precision..precision: constParts;
      array[1..n-1] of var -precision..precision: negHelper;
      array[1..n-1] of var -precision..precision: constHelper;
      array[0..n-1] of var 1..bits: corsLog2;

      constraint negParts[n-1] = 0;
      constraint constParts[n-1] = precision;
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
          2*negParts[i-1] = negHelper[i] /\ corsLog2[i] = bits
          /\
          if x[i-1] == 1 then
            2*constParts[i-1] = 2*constHelper[i] - negHelper[i]
          else
            2*constParts[i-1] = 2*constHelper[i] + negHelper[i]
          endif
        else
          constParts[i-1] = precision /\ negParts[i-1] = 0
          /\
          abs(constHelper[i]) > pow(2, corsLog2[i]-1) /\ abs(constHelper[i]) <= pow(2, corsLog2[i])
        endif
      )
      /\
      abs(constParts[0] + negParts[0]) > pow(2, corsLog2[0]-1)
      /\
      abs(constParts[0] + negParts[0]) <= pow(2, corsLog2[0])
      /\
      corLog2 = sum(i in 0..n-1)(bits - corsLog2[i]);
      
      
predicate difRoundComponent0(array[0..3, 0..63] of var 0..1: inputDif, array[0..3, 0..63] of var 0..1: outputDif, var int: pLog2, int: s, int: t) =
    let {
      var int: upModularAdditionPLog2;
      var int: downModularAdditionPLog2;
    } 
    in
      forall(i in 0..63)(inputDif[0,i] == outputDif[0, (i+s) mod 64])
      /\
      difModularAddition(inputDif[0,..], inputDif[1,..], outputDif[1,..], upModularAdditionPLog2)
      /\
      difModularAddition(inputDif[2,..], inputDif[3,..], outputDif[2,..], downModularAdditionPLog2)
      /\
      forall(i in 0..63)(inputDif[3,i] == outputDif[3, (i+t) mod 64])
      /\
      pLog2 = upModularAdditionPLog2 + downModularAdditionPLog2;
      
predicate difRoundComponent1(array[0..3, 0..63] of var 0..1: inputDif, array[0..3, 0..63] of var 0..1: outputDif) = 
    forall(i in 0..63)(inputDif[0,i] + inputDif[1,i] + outputDif[0,i] in {0,2})
    /\
    forall(i in 0..63)(inputDif[1,i] == outputDif[2,(i+32) mod 64])
    /\
    forall(i in 0..63)(outputDif[1,i] == inputDif[2,i])
    /\
    forall(i in 0..63)(inputDif[2,i] + inputDif[3,i] + outputDif[3,i] in {0,2});
    
    
    
predicate linRoundComponent0(array[0..3, 0..63] of var 0..1: inputDif, array[0..3, 0..63] of var 0..1: outputDif, var int: corLog2, int: s, int: t) = 
    let {
      var int: upModularAdditionCorLog2;
      var int: downModularAdditionCorLog2;
      array[0..63] of var 0..1: upModularAdditionUpInput;
      array[0..63] of var 0..1: downModularAdditionDownInput;
    }
    in
      forall(i in 0..63)(inputDif[0,i] + outputDif[0,(i+s) mod 64] + upModularAdditionUpInput[i] in {0,2})
      /\
      linModularAddition(inputDif[1,..], upModularAdditionUpInput, outputDif[1,..], upModularAdditionCorLog2)
      /\
      linModularAddition(inputDif[2,..], downModularAdditionDownInput, outputDif[2,..], downModularAdditionCorLog2)
      /\
      forall(i in 0..63)(inputDif[3,i] + outputDif[3,(i+t) mod 64] + downModularAdditionDownInput[i] in {0,2})
      /\
      corLog2 = upModularAdditionCorLog2 + downModularAdditionCorLog2;
      
      
predicate linRoundComponent1(array[0..3, 0..63] of var 0..1: inputDif, array[0..3, 0..63] of var 0..1: outputDif) = 
    forall(i in 0..63)
    (
      outputDif[0,i] == inputDif[0,i]
      /\
      inputDif[0,i] + inputDif[1,i] + outputDif[2,(i+32) mod 64] in {0,2}
      /\
      inputDif[2,i] + inputDif[3,i] + outputDif[1,i] in {0,2}
      /\
      inputDif[3,i] == outputDif[3,i]
    );

predicate DLRoundComponent0(array[0..3, 0..63] of var 0..1: inputDif, array[0..3, 0..63] of var 0..1: outputMask, var int: upModularAdditionCor, var int: downModularAdditionCor, int: s, int: t) = 
    DLModularAddition(inputDif[0,..], inputDif[1,..], outputMask[1,..], upModularAdditionCor)
    /\
    DLModularAddition(inputDif[2,..], inputDif[3,..], outputMask[2,..], downModularAdditionCor);
      

predicate difRound(array[0..3, 0..63] of var 0..1: inputDif, array[0..3, 0..63] of var 0..1: outputDif, var int: pLog2) = 
    let {
      array[0..3, 0..63] of var 0..1: afterFirstQuarter;
      array[0..3, 0..63] of var 0..1: afterSecondQuarter;
      array[0..3, 0..63] of var 0..1: afterThirdQuarter;
      var int: leftPLog2;
      var int: rightPLog2;
    }
    in
      difRoundComponent0(inputDif, afterFirstQuarter, leftPLog2, 13, 16)
      /\
      difRoundComponent1(afterFirstQuarter, afterSecondQuarter)
      /\
      difRoundComponent0(afterSecondQuarter, afterThirdQuarter, rightPLog2, 17, 21)
      /\
      difRoundComponent1(afterThirdQuarter, outputDif)
      /\
      pLog2 = leftPLog2 + rightPLog2;
      
predicate linRound(array[0..3, 0..63] of var 0..1: inputMask, array[0..3, 0..63] of var 0..1: outputMask, var int: corLog2) = 
    let {
      array[0..3, 0..63] of var 0..1: afterFirstQuarter;
      array[0..3, 0..63] of var 0..1: afterSecondQuarter;
      array[0..3, 0..63] of var 0..1: afterThirdQuarter;
      var int: leftCorLog2;
      var int: rightCorLog2;
    }
    in
      linRoundComponent0(inputMask, afterFirstQuarter, leftCorLog2, 13, 16)
      /\
      linRoundComponent1(afterFirstQuarter, afterSecondQuarter)
      /\
      linRoundComponent0(afterSecondQuarter, afterThirdQuarter, rightCorLog2, 17, 21)
      /\
      linRoundComponent1(afterThirdQuarter, outputMask)
      /\
      corLog2 = leftCorLog2 + rightCorLog2;
      

predicate rightHalfLinRound(array[0..3, 0..63] of var 0..1: inputMask, array[0..3, 0..63] of var 0..1: outputMask, var int: linCorLog2) = 
    let {
      array[0..3, 0..63] of var 0..1: afterSecondQuarter;
      array[0..3, 0..63] of var 0..1: afterThirdQuarter;
    }      
    in
      linRoundComponent1(inputMask, afterSecondQuarter)
      /\
      linRoundComponent0(afterSecondQuarter, afterThirdQuarter, linCorLog2, 17, 21)
      /\
      linRoundComponent1(afterThirdQuarter, outputMask);


predicate leftHalfDifRound(array[0..3, 0..63] of var 0..1: inputDif, array[0..3, 0..63] of var 0..1: outputDif, var int: pLog2) = 
    let {
      array[0..3, 0..63] of var 0..1: afterFirstQuarter;
    }
    in
      difRoundComponent0(inputDif, afterFirstQuarter, pLog2, 13, 16)
      /\
      difRoundComponent1(afterFirstQuarter, outputDif);

      
int: r0 = 3;
int: r1 = 1;

array[0..r0, 0..255] of var 0..1: difVars;
array[0..r0-1] of var int: pLog2;

array[0..r1, 0..255] of var 0..1: linVars;
array[0..r1-1] of var int: corLog2;

var 0..bits: DLRoundDLCor0Log2;
var 0..bits: DLRoundDLCor1Log2;
var int: DLRoundDifProbLog2;


constraint forall(i in 0..r0-1)(difRound(array2d(0..3,0..63, difVars[i,..]), array2d(0..3,0..63, difVars[i+1,..]), pLog2[i]));
constraint forall(i in 0..r1-1)(linRound(array2d(0..3,0..63, linVars[i,..]), array2d(0..3,0..63, linVars[i+1,..]), corLog2[i]));

array[0..3, 0..63] of var 0..1: afterSecondQuarter;
array[0..3, 0..63] of var 0..1: afterThirdQuarter;
constraint leftHalfDifRound(array2d(0..3,0..63, difVars[r0,..]), afterSecondQuarter, DLRoundDifProbLog2);
constraint DLRoundComponent0(afterSecondQuarter, afterThirdQuarter, DLRoundDLCor0Log2, DLRoundDLCor1Log2, 17,21);
constraint linRoundComponent1(afterThirdQuarter, array2d(0..3,0..63, linVars[0,..]));

array[0..3, 0..n-1] of var 0..1: DLInputDif = array2d(0..3,0..n-1, afterSecondQuarter);
array[0..3, 0..n-1] of var 0..1: DLOutputMask = array2d(0..3,0..n-1, afterThirdQuarter);



constraint sum(i in 0..255)(difVars[0,i]) > 0;
constraint sum(i in 0..255)(linVars[r1,i]) > 0;

constraint forall(i in 0..63)(linVars[r1,i] == linVars[r1,i+64] /\ linVars[r1, i+64] == linVars[r1,i+128] /\ linVars[r1, i+128] == linVars[r1,i+192]);



var int: pU = sum(i in 0..r0-1)(pLog2[i]);
var int: cL = sum(i in 0..r1-1)(2*corLog2[i]);
var int: cM = DLRoundDLCor0Log2 + DLRoundDLCor1Log2;

var 0..100: obj = pU + cL + DLRoundDifProbLog2 + cM;


ann: search_ann;
search_ann = seq_search([
    int_search([linVars[r1,z] | z in 0..63], indomain_min, complete)
    ]);

solve minimize obj;

%output ["\(DLInputDif[0,..]) \n"];
%output ["\(DLInputDif[1,..]) \n"];
%output ["\(DLOutputMask[1,..]) \n"];
%output ["\n"];
%output ["\(DLInputDif[2,..]) \n"];
%output ["\(DLInputDif[3,..]) \n"];
%output ["\(DLOutputMask[2,..]) \n"];

output ["\(difVars[i,..]) -- D \(pLog2[i]) \n " | i in 0..r0-1];
output ["\(difVars[r0,..]) -- D \(DLRoundDifProbLog2) \n"];
output ["\(afterSecondQuarter) -- DL \(DLRoundDLCor0Log2) , \(DLRoundDLCor1Log2) \n"];
output ["\(afterThirdQuarter)"];
output ["\(linVars[i,..]) -- L \(corLog2[i]) \n " | i in 0..r1-1];
output ["\(linVars[r1,..]) \n"];
output ["\(pU),\(DLRoundDifProbLog2) -- \(DLRoundDLCor0Log2),\(DLRoundDLCor1Log2) -- \(cL) \n"];
      

    

    

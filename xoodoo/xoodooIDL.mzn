include "globals.mzn";

array[1..29,1..7] of int: DDT = array2d(1..29, 1..7, [0,0,0,0,0,0,0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]);

array[1..29,1..7] of int: LAT2 = array2d(1..29, 1..7, [0,0,0,0,0,0,0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]);

array[0..11,0..31] of int: RCs = array2d(0..11, 0..31, [0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]);

/*
--------------------------------- Differential part -------------------------------------
*/
predicate difRoundTheta(array[int,int,int] of var 0..1: a, array[int,int,int] of var 0..1: b) = 
    let {
      array[0..3,0..15] of var 0..1: colSum;
      
      constraint forall(x in 0..3, z in 0..15)
      (
        (colSum[x,z]+a[x,0,z]+a[x,1,z]+a[x,2,z]) in {0,2,4}
      );
    } 
    in
      forall(x in 0..3, y in 0..2, z in 0..15) 
      (
        (b[x,y,z]+a[x,y,z]+colSum[(x+3) mod 4,(z+27) mod 16]+colSum[(x+3) mod 4,(z+18) mod 16]) in {0,2,4}
      );
      
predicate difRoundRhoWest(array[int,int,int] of var 0..1: a, array[int,int,int] of var 0..1: b) =
    forall(x in 0..3, z in 0..15)
    (
      b[x,1,z] = a[(x+3) mod 4,1,z] /\ b[x,2,z] = a[x,2,(z+21) mod 16] /\ b[x,0,z] = a[x,0,z]
    );
    
predicate difRoundLota(array[int,int,int] of var 0..1: a, array[int,int,int] of var 0..1: b, array[int] of int: RC) =
    forall(x in 0..3, y in 0..2, z in 0..15)
    (
      if x == 0 /\ y == 0 then
        (b[0,0,z]+a[0,0,z]+RC[z]+RC[z+16]) in {0,2,4}
      else
        b[x,y,z] = a[x,y,z]
      endif
    );

predicate difRoundChi(array[int,int,int] of var 0..1: a, array[int,int,int] of var 0..1: b, array[int,int] of var int: pLog2) =
    forall(x in 0..3, z in 0..15)
    (
      table([a[x,0,z],a[x,1,z],a[x,2,z],b[x,0,z],b[x,1,z],b[x,2,z],pLog2[x,z]], DDT)
    );
    
predicate difRoundRhoEast(array[int,int,int] of var 0..1: a, array[int,int,int] of var 0..1: b) =
    forall(x in 0..3, z in 0..15)
    (
      b[x,1,z] = a[x,1,(z+31) mod 16] /\ b[x,2,z] = a[(x+2) mod 4, 2, (z+24) mod 16] /\ b[x,0,z] = a[x,0,z]
    );
    
/*
--------------------------------- DL differential part -------------------------------------
*/
predicate DLDifRoundTheta(array[int,int,int] of var -1..1: a, array[int,int,int] of var -1..1: b) = 
    let {
      array[0..3,0..15] of var -1..1: colSum;
      
      constraint forall(x in 0..3, z in 0..15)
      (
        if exists(y in 0..2)(a[x, y, z] == -1) then
          colSum[x,z] = -1
        else
          (colSum[x,z]+a[x,0,z]+a[x,1,z]+a[x,2,z]) in {0,2,4} /\ colSum[x,z] >= 0
        endif
      );
    } 
    in
      forall(x in 0..3, y in 0..2, z in 0..15) 
      (
        if a[x,y,z] == -1 \/ colSum[(x+3) mod 4,(z+27) mod 16] == -1 \/ colSum[(x+3) mod 4,(z+18) mod 16] == -1 then
          b[x,y,z] = -1
        else
          (b[x,y,z]+a[x,y,z]+colSum[(x+3) mod 4,(z+27) mod 16]+colSum[(x+3) mod 4,(z+18) mod 16]) in {0,2,4} /\ b[x,y,z] >= 0
        endif
      );
      
predicate DLDifRoundRhoWest(array[int,int,int] of var -1..1: a, array[int,int,int] of var -1..1: b) =
    forall(x in 0..3, z in 0..15)
    (
      (b[x,1,z] = a[(x+3) mod 4,1,z]) /\ (b[x,2,z] = a[x,2,(z+21) mod 16]) /\ b[x,0,z] = a[x,0,z]
    );
    
predicate DLDifRoundLota(array[int,int,int] of var -1..1: a, array[int,int,int] of var -1..1: b, array[int] of int: RC) =
    forall(x in 0..3, y in 0..2, z in 0..15)
    (
      if x == 0 /\ y == 0 then
        if a[0,0,z] == -1 then
          b[0,0,z] == -1
        else
          (b[0,0,z]+a[0,0,z]+RC[z]+RC[z+16]) in {0,2,4} /\ b[0,0,z] >= 0
        endif
      else
        b[x,y,z] = a[x,y,z]
      endif
    );
    
predicate DLDifRoundChi(array[int,int,int] of var -1..1: a, array[int,int,int] of var -1..1: b) = 
    forall(x in 0..3, z in 0..15)
    (
      if a[x,0,z] == 0 /\ a[x,1,z] == 0 /\ a[x,2,z] == 0 then
        b[x,0,z] = 0 /\ b[x,1,z] = 0 /\ b[x,2,z] == 0
      elseif a[x,0,z] == 1 /\ a[x,1,z] == 0 /\ a[x,2,z] == 0 then
        b[x,0,z] = 1 /\ b[x,1,z] = -1 /\ b[x,2,z] == -1
      elseif a[x,0,z] == 0 /\ a[x,1,z] == 1 /\ a[x,2,z] == 0 then
        b[x,0,z] = -1 /\ b[x,1,z] = 1 /\ b[x,2,z] == -1
      elseif a[x,0,z] == 0 /\ a[x,1,z] == 0 /\ a[x,2,z] == 1 then
        b[x,0,z] = -1 /\ b[x,1,z] = -1 /\ b[x,2,z] == 1
      else
        b[x,0,z] = -1 /\ b[x,1,z] = -1 /\ b[x,2,z] == -1
      endif
    );

predicate DLDifRoundRhoEast(array[int,int,int] of var -1..1: a, array[int,int,int] of var -1..1: b) =
    forall(x in 0..3, z in 0..15)
    (
      b[x,1,z] = a[x,1,(z+31) mod 16] /\ b[x,2,z] = a[(x+2) mod 4, 2, (z+24) mod 16] /\ b[x,0,z] = a[x,0,z]
    );
    


/*
--------------------------------- Linear part -------------------------------------
*/
predicate linRoundTheta(array[int,int,int] of var 0..1: a, array[int,int,int] of var 0..1: b) = 
    let {
      array[0..3,0..15] of var 0..1: colSum;
      
      constraint forall(x in 0..3, z in 0..15)
      (
        (colSum[x,z]+b[x,0,z]+b[x,1,z]+b[x,2,z]) in {0,2,4}
      );
    } 
    in
      forall(x in 0..3, y in 0..2, z in 0..15) 
      (
        (a[x,y,z]+b[x,y,z]+colSum[(x+1) mod 4,(z+5) mod 16]+colSum[(x+1) mod 4,(z+14) mod 16]) in {0,2,4}
      );
      
predicate linRoundRhoWest(array[int,int,int] of var 0..1: a, array[int,int,int] of var 0..1: b) =
    forall(x in 0..3, z in 0..15)
    (
      b[x,1,z] = a[(x+3) mod 4,1,z] /\ b[x,2,z] = a[x,2,(z+21) mod 16] /\ b[x,0,z] = a[x,0,z]
    );
    
predicate linRoundChi(array[int,int,int] of var 0..1: a, array[int,int,int] of var 0..1: b, array[int,int] of var int: c2Log2) = 
    forall(x in 0..3, z in 0..15)
    (
      table([a[x,0,z],a[x,1,z],a[x,2,z],b[x,0,z],b[x,1,z],b[x,2,z],c2Log2[x,z]], LAT2)
    );
    
    
predicate linRoundRhoEast(array[int,int,int] of var 0..1: a, array[int,int,int] of var 0..1: b) =
    forall(x in 0..3, z in 0..15)
    (
      b[x,1,z] = a[x,1,(z+31) mod 16] /\ b[x,2,z] = a[(x+2) mod 4, 2, (z+24) mod 16] /\ b[x,0,z] = a[x,0,z]
    );
    
/*
--------------------------------- DL linear part -------------------------------------
*/
predicate DLLinRoundTheta(array[int,int,int] of var -1..1: a, array[int,int,int] of var -1..1: b) = 
    let {
      array[0..3,0..15] of var -1..1: colSum;
      
      constraint forall(x in 0..3, z in 0..15)
      (
        if exists(y in 0..2)(b[x, y, z] == -1) then
          colSum[x,z] = -1
        else
          (colSum[x,z]+b[x,0,z]+b[x,1,z]+b[x,2,z]) in {0,2,4} /\ colSum[x,z] >= 0
        endif
      );
    } 
    in
      forall(x in 0..3, y in 0..2, z in 0..15) 
      (
        if b[x,y,z] == -1 \/ colSum[(x+1) mod 4,(z+5) mod 16] == -1 \/ colSum[(x+1) mod 4,(z+14) mod 16] == -1 then
          a[x,y,z] = -1
        else
          (a[x,y,z]+b[x,y,z]+colSum[(x+1) mod 4,(z+5) mod 16]+colSum[(x+1) mod 4,(z+14) mod 16]) in {0,2,4} /\ a[x,y,z] >= 0
        endif
      );

predicate DLLinRoundRhoWest(array[int,int,int] of var -1..1: a, array[int,int,int] of var -1..1: b) =
    forall(x in 0..3, z in 0..15)
    (
      b[x,1,z] = a[(x+3) mod 4,1,z] /\ b[x,2,z] = a[x,2,(z+21) mod 16] /\ b[x,0,z] = a[x,0,z]
    );
    
predicate DLLinRoundChi(array[int,int,int] of var -1..1: a, array[int,int,int] of var -1..1: b) = 
    forall(x in 0..3, z in 0..15)
    (
      if b[x,0,z] == 0 /\ b[x,1,z] == 0 /\ b[x,2,z] == 0 then
        a[x,0,z] = 0 /\ a[x,1,z] = 0 /\ a[x,2,z] == 0
      elseif b[x,0,z] == 1 /\ b[x,1,z] == 0 /\ b[x,2,z] == 0 then
        a[x,0,z] = 1 /\ a[x,1,z] = -1 /\ a[x,2,z] == -1
      elseif b[x,0,z] == 0 /\ b[x,1,z] == 1 /\ b[x,2,z] == 0 then
        a[x,0,z] = -1 /\ a[x,1,z] = 1 /\ a[x,2,z] == -1
      elseif b[x,0,z] == 0 /\ b[x,1,z] == 0 /\ b[x,2,z] == 1 then
        a[x,0,z] = -1 /\ a[x,1,z] = -1 /\ a[x,2,z] == 1
      else
        a[x,0,z] = -1 /\ a[x,1,z] = -1 /\ a[x,2,z] == -1
      endif
    );
    
predicate DLLinRoundRhoEast(array[int,int,int] of var -1..1: a, array[int,int,int] of var -1..1: b) =
    forall(x in 0..3, z in 0..15)
    (
      b[x,1,z] = a[x,1,(z+31) mod 16] /\ b[x,2,z] = a[(x+2) mod 4, 2, (z+24) mod 16] /\ b[x,0,z] = a[x,0,z]
    );
    

% differential part (2 rounds)


array[0..3,0..2,0..15] of var 0..1: a03;
array[0..3,0..2,0..15] of var 0..1: a04;
array[0..3,0..2,0..15] of var 0..1: a10;
array[0..3,0..2,0..15] of var 0..1: a11;
array[0..3,0..2,0..15] of var 0..1: a12;
array[0..3,0..2,0..15] of var 0..1: a13;
array[0..3,0..2,0..15] of var 0..1: a14;
array[0..3,0..2,0..15] of var 0..1: a20;
array[0..3,0..2,0..15] of var 0..1: a21;
array[0..3,0..2,0..15] of var 0..1: a22;
array[0..3,0..2,0..15] of var 0..1: tmpa23;
array[0..3,0..15] of var 0..3: p0Log2;
array[0..3,0..15] of var 0..3: p1Log2;

constraint difRoundChi(a03, a04, p0Log2);
constraint difRoundRhoEast(a04,a10);

constraint difRoundTheta(a10,a11);
constraint difRoundRhoWest(a11,a12);
constraint difRoundLota(a12,a13,row(RCs,7));
constraint difRoundChi(a13,a14, p1Log2);
constraint difRoundRhoEast(a14,a20);

constraint difRoundTheta(a20,a21);
constraint difRoundRhoWest(a21,a22);
constraint difRoundLota(a22,tmpa23,row(RCs,8));


var int: pU = 2* sum(x in 0..3, z in 0..15)(p0Log2[x,z]+p1Log2[x,z]);






array[0..3,0..2,0..15] of var -1..1: a23;
array[0..3,0..2,0..15] of var -1..1: a24;
array[0..3,0..2,0..15] of var -1..1: a30;
array[0..3,0..2,0..15] of var -1..1: a31;
array[0..3,0..2,0..15] of var -1..1: a32;
array[0..3,0..2,0..15] of var -1..1: a33;
array[0..3,0..2,0..15] of var -1..1: a34;
array[0..3,0..2,0..15] of var -1..1: a40;
array[0..3,0..2,0..15] of var -1..1: a41;
array[0..3,0..2,0..15] of var -1..1: a42;
array[0..3,0..2,0..15] of var -1..1: a43;



constraint forall(x in 0..3, y in 0..2, z in 0..15)(a23[x,y,z] = tmpa23[x,y,z]);


constraint DLDifRoundChi(a23,a24);
constraint DLDifRoundRhoEast(a24,a30);

constraint DLDifRoundTheta(a30,a31);
constraint DLDifRoundRhoWest(a31,a32);
constraint DLDifRoundLota(a32,a33,row(RCs,9));
constraint DLDifRoundChi(a33,a34);
constraint DLDifRoundRhoEast(a34,a40);

constraint DLDifRoundTheta(a40,a41);
constraint DLDifRoundRhoWest(a41,a42);
constraint DLDifRoundLota(a42,a43,row(RCs,10));





array[0..3,0..2,0..15] of var -1..1: b22;
array[0..3,0..2,0..15] of var -1..1: b23;
array[0..3,0..2,0..15] of var -1..1: b30;
array[0..3,0..2,0..15] of var -1..1: b31;
array[0..3,0..2,0..15] of var -1..1: b32;
array[0..3,0..2,0..15] of var -1..1: b33;
array[0..3,0..2,0..15] of var -1..1: b40;
array[0..3,0..2,0..15] of var -1..1: b41;
array[0..3,0..2,0..15] of var -1..1: b42;
array[0..3,0..2,0..15] of var -1..1: b43;




constraint DLLinRoundChi(b22,b23);
constraint DLLinRoundRhoEast(b23,b30);

constraint DLLinRoundTheta(b30,b31);
constraint DLLinRoundRhoWest(b31,b32);
constraint DLLinRoundChi(b32,b33);
constraint DLLinRoundRhoEast(b33,b40);

constraint DLLinRoundTheta(b40,b41);
constraint DLLinRoundRhoWest(b41,b42);
constraint DLLinRoundChi(b42,b43);

var int: CM = sum(x in 0..3, y in 0..2, z in 0..15)(bool2int( b22[x,y,z] == -1 /\ a23[x,y,z] == -1 ))
            + sum(x in 0..3, y in 0..2, z in 0..15)(bool2int( b32[x,y,z] == -1 /\ a33[x,y,z] == -1 ))
            + sum(x in 0..3, y in 0..2, z in 0..15)(bool2int( b42[x,y,z] == -1 /\ a43[x,y,z] == -1 ));

array[0..3,0..15] of var 0..6: c02Log2;
array[0..3,0..2,0..15] of var 0..1: tmpb43;
array[0..3,0..2,0..15] of var 0..1: b50;
array[0..3,0..2,0..15] of var 0..1: b51;
array[0..3,0..2,0..15] of var 0..1: b52;
array[0..3,0..2,0..15] of var 0..1: b53;

constraint linRoundRhoEast(tmpb43,b50);

constraint linRoundTheta(b50,b51);
constraint linRoundRhoWest(b51,b52);
constraint linRoundChi(b52,b53,c02Log2);




constraint forall(x in 0..3, y in 0..2, z in 0..15)(b43[x,y,z] = tmpb43[x,y,z]);

var int: CL = 2*sum(x in 0..3, z in 0..15)(c02Log2[x,z]);

          
     
int: wu = 1;
int: wm = 2;   
int: wl = 1;


constraint (sum(x in 0..3, y in 0..2, z in 0..15)(b53[x,y,z])) > 0;

var int: obj = wu*pU + wm*CM + wl*CL;

ann: search_ann;
search_ann = seq_search([
    int_search([tmpb43[0,0,z] | z in 0..15], indomain_min, complete),
    int_search([tmpb43[0,1,z] | z in 0..15], indomain_min, complete),
    int_search([tmpb43[0,2,z] | z in 0..15], indomain_min, complete),
    int_search([tmpb43[1,0,z] | z in 0..15], indomain_min, complete),
    int_search([tmpb43[1,1,z] | z in 0..15], indomain_min, complete),
    int_search([tmpb43[1,2,z] | z in 0..15], indomain_min, complete),
    int_search([tmpb43[2,0,z] | z in 0..15], indomain_min, complete),
    int_search([tmpb43[2,1,z] | z in 0..15], indomain_min, complete),
    int_search([tmpb43[2,2,z] | z in 0..15], indomain_min, complete),
    int_search([tmpb43[3,0,z] | z in 0..15], indomain_min, complete),
    int_search([tmpb43[3,1,z] | z in 0..15], indomain_min, complete),
    int_search([tmpb43[3,2,z] | z in 0..15], indomain_min, complete)]);



solve minimize obj;


output ["Differential Propagation: \n \(a03) \n \(a04) \n \(a10) \n \(a11) \n \(a12) \n \(a13) \n \(a14) \n \(a23) \n"];
output ["Linear Propagation: \n \(b43) \n \(b53) \n"];
output ["\(pU) -- \(CM) -- \(CL)"];


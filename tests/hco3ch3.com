***,hco3ch3 CAME vibrational spectrum
memory,1000,m

orient,mass
geometry={
O1
C1,O1,B1
O2,C1,B2,O1,A1
H1,O2,B3,C1,A2,O1,0.0
O3,C1,B4,O1,A3,H1,180.0
C2,O3,B5,C1,A4,O1,0.0
H2,C2,B6,O3,A5,C1,180.0
H3,C2,B7,O3,A6,C1,D1
H4,C2,B7,O3,A6,C1,-D1
}
 
B1=                  1.20826386 ANG
B2=                  1.34632938 ANG
A1=                125.37105052 DEGREE
B3=                  0.96620055 ANG
A2=                104.93520844 DEGREE
B4=                  1.33375801 ANG
A3=                126.65441654 DEGREE
B5=                  1.43708653 ANG
A4=                113.41206325 DEGREE
B6=                  1.08375252 ANG
A5=                105.23462880 DEGREE
B7=                  1.08668081 ANG
A6=                110.40633463 DEGREE
D1=                 60.53271261 DEGREE
 
mass,iso
basis=vtz-f12
{hf
start,atden}
ccsd(t)-f12a
optg
{freq,symm=auto, print=0}

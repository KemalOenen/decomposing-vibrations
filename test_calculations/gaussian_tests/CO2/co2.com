%chk=c2h2.chk
#p opt freq hf/3-21g

CO2 Calculation

0 1
X1
C1 X1 1.
O1 C1 r1 X1 90.
O2 C1 r1 X1 90. O1 180.

r1 1.1558

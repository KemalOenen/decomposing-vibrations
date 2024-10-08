%chk=c2h2.chk
#p opt freq hf/6-311g(d,p)

Cs CHOOH

0 1
C1
O1 C1 r1
O2 C1 r2 O1 a1
H1 C1 r3 O1 a2 O2 180.
H2 O2 r4 C1 a3 O1 180.

r1   1.1918
r2   1.3543
r3   1.0820
r4   0.9642
a1 122.5617
a2 123.6452
a3 114.7730

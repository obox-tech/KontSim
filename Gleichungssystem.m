syms ddx ddth1 ddth2 mm m1 m2 l1 l2 F g th2

[ddx,ddth1,ddth2] =solve([ddx*(mm+m1+m2)+ddth1*(m1/2+m2)*l1+ddth2*m2*l2/2==F,+ddx*(m1/2+m2)*l1+ddth1*(m1/3+m2)*l1*l1+ddth2*l1*l2*m2/2==0 ,+ddx*m2*l2/2+ddth1*m2/2*l1*l2+ddth2*m2/3*l2*l2==m2/2*l2*g*th2],[ddx , ddth1, ddth2])

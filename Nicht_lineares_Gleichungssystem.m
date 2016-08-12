syms ddx ddth1 ddth2 mm m1 m2 l1 l2 F g th2 dth1 dth2 th1

[ddx,ddth1,ddth2] =solve([ddx*(mm+m1+m2)+ddth1*(m1/2+m2)*l1*cos(th1)+ddth2*m2*l2/2*cos(th2)==F+(m1/2+m2)*l1*sin(th1)*dth1*dth1+m2/2*l2*cos(th2),
    ddx*(m1/2+m2)*l1*cos(th1)+ddth1*(m1/3+m2)*l1*l1+ddth2*cos(th2-th1)*l1*l2*m2/2==m2/2*l1*l2*sin(th2-th1)*dth2*dth2*(m1/2+m2)*g*l1*sin(th1) ,
    ddx*m2*l2/2*cos(th2)+ddth1*m2/2*l1*l2*cos(th2-th1)+ddth2*m2/3*l2*l2==m2/2*l2*(g*sin(th2)-l1*sin(th2-th1)*dth2*dth2) ],[ddx , ddth1, ddth2])
pretty([ddx,ddth1,ddth2])
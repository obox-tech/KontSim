function xpunkt=RechteSeite(t,x)

%---global definierte variablen uebernommen---
global m1 m2 l1 l2 mm xc g tend


%---Regelungsparameter aus Angabe uebernommen---
kin=1 ; kth1=32.8641 ; kth1d=-0.385336 ; kth2=-51.80481 ;
kth2d=-8.42461 ; kx=-2.769165 ; kxd=-3.334 ;

%---Definition der Kraft---
F=-kx*x(1)+kxd*x(2)+kin*x(3)+kth1*x(4)+kth2*x(6)+kth1d*x(5)+kth2d*x(7) ;


%---Definition der Ableitung vom Zustandsvektor x---
xpunkt=[-x(2);
    F+(m1/2+m2)*l1*sin(x(4))*x(5)^2+m2/2*l2*sin(x(6))*x(7)^2;
    x(1);
    x(5);
    m2/2*l1*l2*sin(x(6)-x(4))*x(7)^2*(m1/2+m2)*g*l1*sin(x(4));
    x(7);
    m2/2*l2*(g*sin(x(6))-l1*sin(x(6)-x(4))*x(7)^2)] ;
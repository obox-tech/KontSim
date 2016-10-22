function ma=Massenmatrix(t,x)

%---global definierte variablen uebernommen---
global m1 m2 l1 l2 mm xc g tend


%---definition der einzelnen Eintraege in der Massenmatrix---
%---mit x(i), das i-te Element im Zustandsvektor verwenden---
ma(1,1)=-1;
ma(2,2)= mm+m1+m2;
ma(2,5)=(m1/2+m2)*l1*cos(x(4));
ma(2,7)=m2/2*l2*cos(x(6));
ma(3,3)=1;
ma(4,4)=1;
ma(5,2)=(m1/2+m2)*l1*cos(x(4));
ma(5,5)=(m1/3+m2)*l1^2;
ma(5,7)=m2/2*l1*l2*cos(x(6)-x(4));
ma(6,6)=1;
ma(7,2)=m2/2*l2*cos(x(6));
ma(7,5)=m2/2*l1*l2*cos(x(6)-x(4));
ma(7,7)=m2/3*l2^2;
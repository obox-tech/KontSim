global m1 m2 l1 l2 mm xc g tend fx
%---globale Variablen definiert---
m1=0.01 ; m2=0.01 ; l1=0.5 ; l2=0.7 ; mm=0.2 ; g=9.81 ; xc=0.2 ; tend=5.00 ; fx=0 ;


p=2*(m1*m2+4*m1*mm+3*m2*mm+m1*m1);


A=[0    -1     0      0       0                    0                   0;
   0     0     0      0       0               3*g*m1*m2/p              0;
   1     0     0      0       0                    0                   0;
   0     0     0      0       1                    0                   0;
   0     0     0      0       0       -9*g*m2*(m1+2*mm)/l1/p           0;
   0     0     0      0       0                    0                   1;
   0     0     0      0       0  3*g*(m1*(m1+4*m2+4*mm)+12*m2*mm)/l2/p 0];

B=[     0;
   (8*m1+6*m2)/p;
        0;
        0;
  -3*(4*m1+2*m2)/l1/p;
        0;
     6*m1/l2/p];
 
C=[1 0 0 0 0 0 0;
   0 0 0 1 0 0 0;
   0 0 0 0 0 1 0];

D=[ 0;
    0;
    0];


states={'e' 'dx' 'in' 'th1' 'dth1' 'th2' 'dth2'};
inputs={'F'};
outputs={'e' 'th1' 'th2'};

sys_ss=ss(A,B,C,D,'statename', states, 'inputname' , inputs, 'outputname' , outputs);

%co=ctrb(sys_ss);
%controllability=rank(co)


Q=eye(7);
r=1;

k=lqr(A,B,Q,r)

Ac = [(A-B*k)];
Bc = [B];
Cc = [C];
Dc = [D];


sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:5;
r =0.2*ones(size(t));
[y,t,x]=lsim(sys_cl,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),t,y(:,3),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
set(get(AX(3),'Ylabel'),'String','2nd pendulum angle (radians)')
title('Step Response with LQR Control')
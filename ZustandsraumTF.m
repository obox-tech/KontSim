syms ddx ddth1 ddth2 th2 X TH1 TH2
s=tf('s');
m1=0.01 ; m2=0.01 ; l1=0.5 ; l2=0.7 ; mm=0.2 ; g=9.81;


A=mm+m1+m2;
B=m1/3+m2;
C=m2/3*l2^2;
D=(m1/2+m2)*l1;
E=m2/2*l2;
G=m2/2*l1*l2;
H=m2/2*l2*g;

X_F=(2*(6*g*m1 + 18*g*m2 - 4*l2*m1*s^2 - 12*l2*m2*s^2 + 9*l1^2*l2*m2*s^2))/(s^2*(6*l2*l1^2*m1^2*s^2 - 9*g*l1^2*m1^2 + 24*l2*l1^2*m1*m2*s^2 - 36*g*l1^2*m1*m2 + 6*l2*l1^2*m2^2*s^2 - 36*g*l1^2*m2^2 + 18*l2*mm*l1^2*m2*s^2 - 8*l2*m1^2*s^2 + 12*g*m1^2 - 26*l2*m1*m2*s^2 + 48*g*m1*m2 - 8*l2*mm*m1*s^2 + 12*g*mm*m1 - 6*l2*m2^2*s^2 + 36*g*m2^2 - 24*l2*mm*m2*s^2 + 36*g*mm*m2));
TH1_F=-(6*l1*(3*g*m1 + 6*g*m2 - 2*l2*m1*s^2 - l2*m2*s^2))/(s^2*(6*l2*l1^2*m1^2*s^2 - 9*g*l1^2*m1^2 + 24*l2*l1^2*m1*m2*s^2 - 36*g*l1^2*m1*m2 + 6*l2*l1^2*m2^2*s^2 - 36*g*l1^2*m2^2 + 18*l2*mm*l1^2*m2*s^2 - 8*l2*m1^2*s^2 + 12*g*m1^2 - 26*l2*m1*m2*s^2 + 48*g*m1*m2 - 8*l2*mm*m1*s^2 + 12*g*mm*m1 - 6*l2*m2^2*s^2 + 36*g*m2^2 - 24*l2*mm*m2*s^2 + 36*g*mm*m2));
TH2_F=(6*(2*m1 + 6*m2 - 3*l1^2*m1 - 6*l1^2*m2))/(6*l2*l1^2*m1^2*s^2 - 9*g*l1^2*m1^2 + 24*l2*l1^2*m1*m2*s^2 - 36*g*l1^2*m1*m2 + 6*l2*l1^2*m2^2*s^2 - 36*g*l1^2*m2^2 + 18*l2*mm*l1^2*m2*s^2 - 8*l2*m1^2*s^2 + 12*g*m1^2 - 26*l2*m1*m2*s^2 + 48*g*m1*m2 - 8*l2*mm*m1*s^2 + 12*g*mm*m1 - 6*l2*m2^2*s^2 + 36*g*m2^2 - 24*l2*mm*m2*s^2 + 36*g*mm*m2);


sys_tf = [X_F ; TH1_F ; TH2_F];

inputs={'F'};
outputs={'e' 'th1' 'th2'};
states={'e' 'dx' 'in' 'th1' 'dth1' 'th2' 'dth2'};

set(sys_tf,'InputName',inputs)
set(sys_tf,'OutputName',outputs)

sys_tf

sys_ss=ss(sys_tf);

sys_ss

%set(sys_ss,'statename', states)
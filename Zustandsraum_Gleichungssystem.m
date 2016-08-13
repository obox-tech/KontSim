syms ddx ddth1 ddth2 mm m1 m2 l1 l2 F g th2 X TH1 TH2 X_F TH1_F TH2_F s

A=mm+m1+m2;
B=m1/3+m2;
C=m2/3*l2^2;
D=(m1/2+m2)*l1;
E=m2/2*l2;
G=m2/2*l1*l2;
H=m2/2*l2*g;

[X_F,TH1_F,TH2_F] =solve([A*s^2*X_F+D*s^2*TH1_F+E*s^2*TH2_F==1,D*s^2*X_F+B*s^2*TH1_F+G*s^2*TH2_F==0,E*s^2*X_F+G*s^2*TH1_F+C*s^2*TH2_F==H*TH2_F],[X_F,TH1_F,TH2_F])
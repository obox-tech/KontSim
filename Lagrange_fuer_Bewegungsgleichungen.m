syms l1 l2 phi_1 phi_2 phi_p1 phi_p2 phi_pp1 phi_pp2
syms a a_p a_pp mm m1 m2 g I_1 I_2 F

frg=3;                                     %Anzahl der Freiheitsgrade
n=3;                                       %Anzahl der Koerper

q=[a ; phi_1 ; phi_2];                     %Minimalkoordinaten
q_p=[a_p ; phi_p1 ; phi_p2];               %zeitliche Ableitungen
q_pp=[a_pp ; phi_pp1 ; phi_pp2];

%---- Drehmatrix Stab 1
T_IK1 = [cos(phi_1) sin(phi_1) 0;
        -sin(phi_1) cos(phi_1) 0;
              0          0     1];
%---- Drehmatrix Stab 2
T_IK2 = [cos(phi_2) sin(phi_2) 0;
        -sin(phi_2) cos(phi_2) 0;
              0          0     1];

%---- Ortsvektoren
I_r_Sm = [a;0;0];
I_r_S1 = [a+l1/2*sin(phi_1) ; l1/2*cos(phi_1) ; 0];
I_r_Q2 = [a+l1*sin(phi_1) ; l1*cos(phi_1) ; 0];
K1_r_Q1S1 = [0; l1/2; 0];
K2_r_Q2S2 = [0; l2/2; 0];
I_r_S2 = I_r_Q2 + T_IK2 * K2_r_Q2S2;

%---- Traegheitstensoren in den koerperfesten Koordinatensystemen
K1_I_S1 = diag([0 0 I_1]);
K2_I_S2 = diag([0 0 I_2]);

%---- Winkelgeschwindigkeitsvektoren der Staebe
K_om1 = [0 ; 0 ; -phi_p1];
K_om2 = [0 ; 0 ; -phi_p2];

%---- JACOBI-Matrizen der Translation
J_Tm = jacobian(I_r_Sm, q);
J_T1 = jacobian(I_r_S1, q);
J_T2 = jacobian(I_r_S2, q);

%---- JACOBI-Matrizen der Rotation
J_R1 = jacobian(K_om1, q_p);
J_R2 = jacobian(K_om2, q_p);

%---- Geschwindigkeitsvektoren
I_v_Sm = J_Tm*q_p ;
I_v_S1 = J_T1*q_p ; 
I_v_S2 = J_T2*q_p ;

%---- kinetische Energie
T = 1/2*(mm*(I_v_Sm.'*I_v_Sm)+m1*(I_v_S1.'*I_v_S1)+m2*(I_v_S2.'*I_v_S2) ...%Translation
    +K_om1.'*K1_I_S1*K_om1+K_om2.'*K2_I_S2*K_om2);           %Rotation
T = simplify(T);                                             %Vereinfachung

%---- potentielle Energie
V=-(m1*I_r_S1.'+m2*I_r_S2.')*[0 ; -g ; 0];

%---- Ableitungen fuer LAGRANGEsche Gleichung 2. Art
dTdv = simplify(jacobian(T,q_p).');         %mit transponieren zu Spaltenvektor gemacht
dTdq = simplify(jacobian(T,q).');
dVdq = simplify(jacobian(V,q).');

%---- Elemente der Bewegungsgleichung M(q)*q_pp + f(q,q_p) = 0
disp('System-Massenmatrix M')
M = simplify(jacobian(dTdv,q_p))
disp('System-Vektorfunktion f')
f = simplify(jacobian(dTdv,q)*q_p+dVdq-dTdq-[F;0;0])

%==========================================================================
%---- Linearisierung um die Gleichgewichtslage:
%     phi_1 = 0, phi_2 = 0, a = 0

disp(' ')
disp('Elemente der linearisierten Bewegungsgleichung')
disp('System-Massenmatrix M0')
M0 = subs(M,{phi_1, phi_2, a},{0, 0, 0})
f0 = subs(f,{a, phi_1, phi_2, a_p, ...
    phi_p1, phi_p2},{0, 0, 0, 0, 0, 0});
disp('Auslenkungs-proportionaler Anteil')
Q = subs(jacobian(f,q),{a, phi_1, phi_2, a_p, ...
    phi_p1, phi_p2},{0, 0, 0, 0, 0, 0})
disp('Steifigkeitsmatrix K')
K = 1/2*(Q+Q.')
disp('Matrix der nichtkonservativen Kräfte')
N = 1/2*(Q-Q.')
disp('gesschw.-proportionaler Anteil')
P = subs(jacobian(f,q_p),{a, phi_1, phi_2, a_p, ...
    phi_p1, phi_p2},{0, 0, 0, 0, 0, 0})
disp('Daempfungsmatrix')
D = 1/2*(P+P.')
disp('gyroskopischer Anteil')
G = 1/2*(P-P.')

%==========================================================================
%----Erstellen und Simulieren der Zustandsraumdarstellung
% ohne Integral
syms x th1 th2 x_p th1_p th2_p
syms x_pp th1_pp th2_pp

y = [q.',q_p.'].';
y_p = [q_p.',x_pp , th1_pp , th2_pp].';

A = [zeros(3),eye(3);
    -M0^(-1)*Q, -M0^(-1)*P];
A = double(subs(A,{mm, m1, m2, l1, l2, g, I_1, I_2}, ...
    {0.2, 0.01, 0.01, 0.5, 0.7, 9.81, 2.0833e-04, 4.0833e-04}))

B = [zeros(3,1);M0^(-1)*[1;0;0]];
B = double(subs(B,{mm, m1, m2, l1, l2, g, I_1, I_2}, ...
    {0.2, 0.01, 0.01, 0.5, 0.7, 9.81, 2.0833e-04, 4.0833e-04}))

C = [1 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1 0 0 0]

D = [0; 0; 0]


Q=eye(6);
r=1;

k = lqr(A,B,Q,r)

Ac = [(A-B*k)];
Bc = [B];
Cc = [C];
Dc = [D];

states = {'x' 'th1' 'th2' 'x_p' 'th1_p' 'th2_p'};
inputs = {'F'};
outputs = {'x' 'th1' 'th2'};

sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:5;
F =0.2*ones(size(t));
[y,t,x]=lsim(sys_cl,F,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
hold on
line(t,y(:,3),'parent',AX(2),'color','g')
hold off
set(get(AX(1),'Ylabel'),'String','cart position (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angles (radians)')
title('Step Response with LQR Control')


syms l1 l2 phi_1 phi_2 phi_p1 phi_p2 phi_pp1 phi_pp2
syms y y_p y_pp mm m1 m2 g I_1 I_2

frg=3;                                     %Anzahl der Freiheitsgrade
n=3;                                       %Anzahl der Koerper

q=[y ; phi_1 ; phi_2];                     %Minimalkoordinaten
q_p=[y_p ; phi_p1 ; phi_p2];               %zeitliche Ableitungen
q_pp=[y_pp ; phi_pp1 ; phi_pp2];

%---- Drehmatrix Stab 1
T_IK1 = [cos(phi_1) sin(phi_1) 0;
        -sin(phi_1) cos(phi_1) 0;
              0          0     1];
%---- Drehmatrix Stab 2
T_IK2 = [cos(phi_2) sin(phi_2) 0;
        -sin(phi_2) cos(phi_2) 0;
              0          0     1];

%---- Ortsvektoren
I_r_Sm = [y;0;0];
I_r_S1 = [y+l1/2*sin(phi_1) ; l1/2*cos(phi_1) ; 0];
I_r_Q2 = [y+l1*sin(phi_1) ; l1*cos(phi_1) ; 0];
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
T= 1/2*(mm*(I_v_Sm.'*I_v_Sm)+m1*(I_v_S1.'*I_v_S1)+m2*(I_v_S2.'*I_v_S2) ...%Translation
    +K_om1.'*K1_I_S1*K_om1+K_om2.'*K2_I_S2*K_om2);           %Rotation
T=simplify(T);                                             %Vereinfachung

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
f = simplify(jacobian(dTdv,q)*q_p+dVdq-dTdq)


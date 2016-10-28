%---- Ermitteln der Bewegungsgleichungen
%     definieren der Systemvariablen
syms l1 l2 th1 th2 th1_p th2_p th1_pp th2_pp
syms x x_p x_pp mm m1 m2 g I_1 I_2 F xc

frg=3;                                     %Anzahl der Freiheitsgrade
n=3;                                       %Anzahl der Koerper

q=[x ; th1 ; th2];                     %Minimalkoordinaten
q_p=[x_p ; th1_p ; th2_p];               %zeitliche Ableitungen
q_pp=[x_pp ; th1_pp ; th2_pp];

%---- Drehmatrix Stab 1
T_IK1 = [cos(th1) sin(th1) 0;
        -sin(th1) cos(th1) 0;
              0          0     1];
%---- Drehmatrix Stab 2
T_IK2 = [cos(th2) sin(th2) 0;
        -sin(th2) cos(th2) 0;
              0          0     1];

%---- Ortsvektoren
I_r_Sm = [x;0;0];
I_r_S1 = [x+l1/2*sin(th1) ; l1/2*cos(th1) ; 0];
I_r_Q2 = [x+l1*sin(th1) ; l1*cos(th1) ; 0];
K1_r_Q1S1 = [0; l1/2; 0];
K2_r_Q2S2 = [0; l2/2; 0];
I_r_S2 = I_r_Q2 + T_IK2 * K2_r_Q2S2;

%---- Traegheitstensoren in den koerperfesten Koordinatensystemen
K1_I_S1 = diag([0 0 I_1]);
K2_I_S2 = diag([0 0 I_2]);

%---- Winkelgeschwindigkeitsvektoren der Staebe
K_om1 = [0 ; 0 ; -th1_p];
K_om2 = [0 ; 0 ; -th2_p];

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
%     th1 = 0, th2 = 0, x = 0

disp(' ')
disp('Elemente der linearisierten Bewegungsgleichung')
disp('System-Massenmatrix M0')
M0 = subs(M,{th1, th2, x},{0, 0, 0})
f0 = subs(f,{x, th1, th2, x_p, ...
    th1_p, th2_p},{0, 0, 0, 0, 0, 0});
disp('Auslenkungs-proportionaler Anteil')
Q = subs(jacobian(f,q),{x, th1, th2, x_p, ...
    th1_p, th2_p},{0, 0, 0, 0, 0, 0})
disp('Steifigkeitsmatrix K')
K = 1/2*(Q+Q.')
disp('Matrix der nichtkonservativen Kraefte')
N = 1/2*(Q-Q.')
disp('gesschw.-proportionaler Anteil')
P = subs(jacobian(f,q_p),{x, th1, th2, x_p, ...
    th1_p, th2_p},{0, 0, 0, 0, 0, 0})
disp('Daempfungsmatrix')
D = 1/2*(P+P.')
disp('gyroskopischer Anteil')
G = 1/2*(P-P.')

%==========================================================================
%----Erstellen und Simulieren der Zustandsraumdarstellung

A = [zeros(3),eye(3);
    -M0^(-1)*Q, -M0^(-1)*P];
A = double(subs(A,{mm, m1, m2, l1, l2, g, I_1, I_2}, ...
    {0.2, 0.01, 0.01, 0.5, 0.7, 9.81, 2.0833e-04, 4.0833e-04}));
A(7,7) = 0;
A(7,1) = -1

B = [zeros(3,1);M0^(-1)*[1;0;0]];
B = double(subs(B,{mm, m1, m2, l1, l2, g, I_1, I_2}, ...
    {0.2, 0.01, 0.01, 0.5, 0.7, 9.81, 2.0833e-04, 4.0833e-04}));
B(7,1) = 0
Bxc = [0; 0; 0; 0; 0; 0; 1]

C = [1 0 0 0 0 0 0;
   0 1 0 0 0 0 0;
   0 0 1 0 0 0 0]

D = [0; 0; 0]


Q=eye(7);
r=1;

%----lqr Regelungsentwurf
k = lqr(A,B,Q,r)

%----neue Zustandsraumsystemmatrizen nach Parameterruekfuehrung
Ac = [(A-B*k)];
Bc = [Bxc];
Cc = [C];
Dc = [D];

states = {'x' 'th1' 'th2' 'x_p' 'th1_p' 'th2_p' 'in'};
inputs = {'F'};
outputs = {'x' 'th1' 'th2'};

sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

%----definieren des Simulationszeitraums
t = 0:0.01:8;

%----definition des konstanten 0.2m offsets als Input
u =0.2*ones(size(t));

%----Simulation des erstellten Systems ueber gegebene Zeit mit bekanntem
%Input
[y,t,x]=lsim(sys_cl,u,t);


%----Drei einzelne Diagramme in einem Fenster
figure(1);
ax(1) = subplot(3,1,1);
    plot(ax(1),t,y(:,1),'b');
    title(ax(1),'cart position');     %Titel, Beschriftungen, Kommentare,
    ylim([-0.1,0.25]);                %andere Farben, andere skalierungen,
    grid on                           %da kann man sich noch frei austoben.
ax(2) = subplot(3,1,2);               %relativ einfach verstaendliche 
    plot(ax(2),t,y(:,2),'r');         %loesung. Ws nicht Laufzeit optimiert
    title(ax(2),'angle th 1');
    grid on
ax(3) = subplot(3,1,3);
    plot(ax(3),t,y(:,3),'g');
    title(ax(3),'angle th 2');
    grid on

%----Plotten der Ausgangsgroessen
% [AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
% hold on
% line(t,y(:,3),'parent',AX(2),'color','g')
% hold off
% set(get(AX(1),'Ylabel'),'String','cart position (m)')
% set(get(AX(2),'Ylabel'),'String','pendulum angles (radians)')
% title('Step Response with LQR Control')

%----Berechnung der Eigenwerte
Eigenwerte = eig(Ac)
function Geregeltes_Doppelpendel(th1ic,th2ic)

% if nargin == 2                          % Start mit vorgegebenen Anfangswerten?
%    x0 = [0; 0; 0; th1ic ; 0; th2ic ; 0 ] ;   % Anfangswerte  [phi1 ; phi2 ; omega1 ; omega2]
% else
%    x0 = [0; 0; 0; 0 ; 0; 0 ; 0 ] ;            % Default-Anfangswerte
% end                                          % Hab ich irgendwo im Internet gefunden, so koennen wir die funktion spaeter aufrufen mit verschiedenen Startwinkeln.

%---globale Variablen angekuendigt---
global m1 m2 l1 l2 mm xc g tend
%---globale Variablen definiert---
m1=0.01 ; m2=0.01 ; l1=0.5 ; l2=0.7 ; mm=0.2 ; g=9.81 ; xc=0.2 ; tend=5.00 ;


%---Zustandsvektor x[e;dx;in;th1;dth1;th2;dth2]
x0=[0; 0; 0; 0; 0; 0; 0];
tspan = [0 ; tend];

%---mit 'options' die Massenmatrix auf der linken Seite
%   zu beruecksichtigen gegeben, 'MaxStep' mit 0.01 damit die Genauigkeit
%   gross genug ist---
options = odeset('Mass' , @Massenmatrix , 'MaxStep' , 0.01) ;
[t x]  = ode45 (@RechteSeite , tspan , x0 , options) ; %---integrieren der DGL mit ode45---


%---auslesen der Werte aus dem Zustandsvektor x---
y=xc-[1 0 0 0 0 0 0]*x';
th1=[0 0 0 1 0 0 0]*x';
th2=[0 0 0 0 0 1 0]*x';

%---zeichnen der Werte ueber der Zeit---
subplot(3,1,1) ; plot (t , y) , grid on , title  ('y(t) (Schlitten)')
subplot(3,1,2) ; plot (t , th1) , grid on , title  ('th1(t) (unteres Pendel)')
subplot(3,1,3) ; plot (t , th2) , grid on , title  ('th2(t) (oberes Pendel)')

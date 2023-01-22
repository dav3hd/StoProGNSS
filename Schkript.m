%% Stochastische Prozesse Beleg 4 - ﻿Zeitreihenanalyse von GNSS-Daten mit power-law noise
close all; clear; clc; format long

dat = readmatrix("aufgabe5.dat");
u = dat(:,3);
t = dat(:,2);

clearvars dat

%% Aufgabe 1: Sicherstellen, dass u(t) eine gerade Anzahl an Werten enthält

if mod(length(u),2)==0
else
u = u(1:end-1); % Kürzen der Zeitreihe um einen Wert
end

figure(1)
hold on
plot(t,u,'b.-')
title("GNSS Up-Reihe")
xlabel("Zeit [t]")
ylabel("u[t]")
hold off

%% Aufgabe 2: Bestimmen der Kenngrößen

I = length(u); % Länge der Zeitreihe
dt = t(2)-t(1);
T = I*dt; % Abtastweite im Frequenzbereich
fg = pi/dt; % Grenzfrequenz oder Nyquistfrequenz

%% Aufgabe 3: Erste Parameterschätzung

L = u;

% Naeherungsvektor der Unbekannten:
u0_ = 0;
u1_ = 1;
a1_ = 1;
a2_ = 1;
b1_ = 1;
b2_ = 1;

X0 = [u0_, u1_, a1_, b1_, a2_, b2_];

A = [];
    A(:,1) = ones(length(u),1);
    A(:,2) = t;
    A(:,3) = cos(2*pi.*t/365.25);
    A(:,4) = sin(2*pi.*t/365.25);
    A(:,5) = cos(4*pi.*t/365.25);
    A(:,6) = sin(4*pi.*t/365.25);


Li = ones(1000,1)*9999;
Ldach = ones(1000,1)*-9999;
while round(Li,6) ~= round(Ldach,6)

L0 = [];
L0(:,1) = u0_ + u1_.*t + a1_*cos(2*pi.*t/365.25) + a2_*sin(2*pi.*t/365.25) + b1_*cos(4*pi.*t/365.25) + b2_*sin(4*pi.*t/365.25);

% gekuerzter Beobachtungsvektor:
l = L - L0;

% P = Einheitsmatrix

% Ausgleichung
% ------------

Xdach = A\u;

% Ausgeglichene Beobachtungen
xdach = (A'*A)\A'*l;
v = A*xdach - l;
Ldach = L + v;

% Hauptrechenprobe
Xdach = [ X0 ]' + xdach; % Ausgeglichene Parameter


u0_ = Xdach(1);
u1_ = Xdach(2);
a1_ = Xdach(3);
a2_ = Xdach(4);
b1_ = Xdach(5);
b2_ = Xdach(6);

Li = [];
Li(:,1) = u0_ + u1_.*t + a1_*cos(2*pi.*t/365.25) + a2_*sin(2*pi.*t/365.25) + b1_*cos(4*pi.*t/365.25) + b2_*sin(4*pi.*t/365.25);

if round(Li,6) == round(Ldach,6)
    fprintf('Hauptrechenprobe war erfolgreich.\n')
else
    fprintf('Hauptrechenprobe ist fehlgeschlagen.\n')
end

end

figure
hold on
plot(t,u,'b.-')
plot(t,Li,'r-',LineWidth=3)
hold off


%% Aufgabe 4: Darstellen der Zeitreihe der Residuen

e = v; % Berechnen der Residuen (Verbesserungen aus Ausgleichung)

figure;
plot(1:I,e); % Diagramm (t, e(t))
xlabel('t'); ylabel('e(t)');
title('Diagramm der Residuen');


%% Aufgabe 5: Schätzen der Energiedichte P(f)

clearvars -except e u t I dt T fg Xdach A Li

E = fft(e);
E = E.*conj(E)/I;

P(1) = E(1);
P(2:I/2) = 2*E(2:I/2);
P(I/2+1) = E(I/2+1);

domega = 2*pi/T;
n = 1:round(fg/domega);
f(:,1) = n .* domega; 
f_deg = f * 180/pi;

figure;
plot(f_deg,P(2:end)); % Diagramm (f, P(f))
xlabel('f'); ylabel('P(f)');
title('Diagramm der Energiedichte');


%% Aufgabe 6 - Bestimmung Spektralindex

P_t = log( P(2:end) );
f_t = log( f );

Pp = polyfit(f_t,P_t,1);

fprintf("Kappa: %.4f\n", Pp(1))
fprintf("Kappa0: %.4f\n", Pp(2))



figure
loglog(f_t,P_t,'.')



%% Aufgabe 7 - Kovarianzmatrix des power-law noise

U = eye(I);
h1 = 1;
for i = 2:I
    h = (h1/(i-1)) * (i - (Pp(1)/2) - 2);
    for j = 1:(I-i+1)
        U(j,j+i-1) = h;
    end
    h1 = h;
end
c0 = 1;
CPL = c0 * U' * U;


figure
imagesc(CPL)

%% Aufgabe 8 - Wiederholung Parameterschätzung

C_PL = CPL;
s0 = std(u-A*Xdach);
C = eye(I)*s0^2 + C_PL;
iC = C \ eye(length(C));
x_cov = (A' * iC * A) \ A' * iC * u;
x_cov = [x_cov(1:2); x_cov(3:4); x_cov(5:6)];

%*365.25/(2*pi)

disp('Die Werte der geschätzten Parameter und ihre Genauigkeiten sowie die Standardabweichung der Gewichtseinheit sind:')
disp(x_cov)

C_x = (A'* iC *A) \ eye(width(A));
C_x = C_x./(C_x(1,1)*C_x(2,2)*C_x(3,3)*C_x(4,4)*C_x(5,5)*C_x(6,6))^(1/6);

figure
imagesc(C_x)
title('Korrelationsmatrix der geschätzten Parameter')
xlabel('Zeilenindex')
ylabel('Spaltenindex')
colorbar

u_fit = A*x_cov;

figure
plot(t,u)
hold on
plot(t,u_fit,'r-')
title('Zeitreihe u(t) und Trajektorie des geschätzten Modells')
xlabel('Zeit, t (d)')
ylabel('Höhe, u(t) (mm)')
legend('u(t)','Modell')


disp('Vergleich der Parameter unter 3. und 8.:')
disp('Parameter 3. 8.')
disp([Xdach x_cov])

figure
hold on
grid on
plot(t,u_fit-Li,'r-')
title("Differenz zwischen Zeitreihe aus Aufgabe 3 und Aufgabe 8")
hold off


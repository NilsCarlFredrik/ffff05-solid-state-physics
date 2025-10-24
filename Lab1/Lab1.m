%% Metallstav

Im = [0.516 0.1 0.996 1.494 2.015 2.490 2.963 3.488 3.888];
Um = [0.303 0.061 0.597 0.931 1.33 1.73 2.27 3.00 3.75];
Tm = [296 298 298 305 312 320 343 375 407]; % Oeversatt termisk spaenning fraan tabell
R = Um./Im;

figure(1)
plot(Tm,R, 'o-')
grid on
xlabel('T') %[K]
ylabel('R') %[Ohm]

%% Halvledare rum 1
format long
data = importdata('kisel.txt');
Th = data(1:end,1);
Uh = data(1:end,2);

lnU = log(1./Uh);
T = 1./Th;

figure(2)
plot(T,lnU)
grid on
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([3*10^-3, 15*10^(-3)]);
ylim([0 9])
xlabel('1/T') % [1/K]
ylabel('ln(1/U)') % [1/m?V]

%% Linjeanpassning, ber?kning av Ea
% Anpassar linjer f?r de tv? omr?dena f?r att finna T1, dvs vid det T som 

Tion = T(T>8*10^(-3) & T<10*10^(-3));                         % range of linear sequence for low T, ionization region
% Tion = T(T>5.5*10^(-3) & T<7*10^(-3));                           
Uion = lnU(strfind(T',Tion'):strfind(T',Tion')+length(Tion)-1); % Corresponding range for U values in relation to Tlag
kmion = polyfit(Tion,Uion,1);
ionLine = @(T) kmion(1).*T + kmion(2);
% lowLineT = arrayfun(lowLine,Tlow);

Tsat = T(T>3*10^(-3) & T<4.5*10^(-3));                          % range of saturation sequence
Usat = lnU(strfind(T',Tsat'):strfind(T',Tsat')+length(Tsat)-1); % Corresponding range for U values in relation to Tlag
kmsat = polyfit(Tsat,Usat,1);
satLine = @(T) kmsat(1).*T + kmsat(2);

%Ber?kning av aktiveringsenergi och dopkoncentration
k = 1.38064852*10^(-23);            %Boltzmann
Eaj = -2*k*kmion(1);                % Aktiveringsenergi [J]
Eaev = Eaj/(1.6022*10^(-19));       % Aktiveringsenergi [Ev]

mh = 0.69;
T1 = 1/0.005591;    % Temperature at line intersection
Nv = 2.5*10^(25)*(mh*T1/300)^(3/2); % Effektiv tillstaandstaethet
Na = Nv*exp(-Eaj/(k*T1));           % Dopkoncentration

%% Plots linjeanpassning halvledare

figure(3)
hold off
plot(Tion,ionLine(Tion))
hold on
plot(Tion,Uion)
hold off
grid on


figure(4)
hold off
plot(Tsat,satLine(Tsat))
hold on
plot(Tsat,Usat)
hold off
grid on

% Plottar m?tv?rden tsm med linjeapproximationer
figure(5)
plot(T,lnU)
hold on
plot(T,satLine(T))
plot(T,ionLine(T))
hold off
grid on
xlabel('1/T')
ylabel('ln(1/U)')
xlim([3*10^-3, 14*10^(-3)]);
ylim([6.5 9.5])

%% Halvledare rum 2 (Halleffekt)

h = 1.5*10^(-6);
e = 1.602176634*10^(-19);
Ih = 1*10^(-3);

% Maetning 1
Is1 = [2.981 2.808 2.608 2.398 2.192 2.006 1.799 13603 1.408 1.200 0.994 0.798 0.603 0.395 0.220 0]; % [A]
Vp1 = -[0.207 0.255 0.308 0.377 0.462 0.545 0.666 0.802 0.933 1.077 1.232 1.380 1.524 1.676 1.810 1.962]; % [mV]
B1 = 0.1.*[6.2 6.31 5.9 5.7 5.3 5 4.8 4.1 3.7 3.1 2.7 2.1 1.6 1.0 0.5 0]; % [Te-1]

% Maetning 2
Is2 = [2.862 2.214 1.645 0.912 0.639 0.123 0]; % [A]
Vp2 = -[3.733 3.516 3.196 2.677 2.465 2.101 2.011]; % [mV
B2 = -0.1.*[6.3 5.3 4.2 2.5 0.5 0.2 0]; % [Te-1]

%linjeanpassning 
kmHall1 = polyfit(B1,Vp1,1);
n1 = Ih/(h*e*kmHall1(1));
% p1 = Ih.*B1./(Vp1.*h.*e); %laddningsbaerare FEL
hallLine = @(B) kmHall1(1).*B + kmHall1(2);

kmHall2 = polyfit(B2,Vp2,1);
% p2 = Ih.*B2./(Vp2.*h.*e);
n2 = Ih/(h*e*kmHall2(1));

%% Plot Hall
figure(6)
hold off
plot(B1,Vp1, 'x')
hold on
plot(B2,Vp2, 'x')
% title('Test 1')
xlabel('B')
ylabel('U')
legend('1','2')
grid on
hold off

%%

figure(7)
hold off
plot(B1,Vp1, 'o')
hold on
plot(B1,hallLine(B1))
% title('Test 1')
xlabel('B')
ylabel('U')
legend('M?tdata','Linjeanpassning')
grid on
hold off



















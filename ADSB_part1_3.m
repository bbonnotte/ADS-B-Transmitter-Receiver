%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 Benjamin BONNOTTE & Abdoul GUISSET - G2                 %
%                     Projet Communication Numerique                      %
%                                 --------                                %
%                SIMULATION D'UN EMETTEUR / RECEPTEUR ADS-B               %
%                              ADSB_part1_3                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all; close all; clc;

%% Definitions des variables

Fe = 20*10^(6);   %Frequence d'echantillonage
Te = 1/Fe;  % Periode d'echantillonage
Ds = 10*10^(6);   % Debit symbole
Ts = 1/Ds;  % Temps symbole
Fse = 20;   % Facteur de sur-echantillonage
Fp = 1090*10^(6);   %Frequence porteuse des signeaux ADS-B
Nfft = 512; % Nombre de points fft
t=0:Te:50*Ts-Te;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=2; % Nombre de symboles dans la constellation
Npaq = 112; % Nombres de bits par paquets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eb_No = 0:1:10; % SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tp = 8; % Duree du signal de synchro sp (en microsecondes)

%% Emetteur
% Bits aleatoires
sb = randi([0 1],1,Npaq); % Generation d'une sequence aleatoire de 1000 bits

% Association bits -> Symboles
ss = pammod(sb,M);

% Sur-echantillonage afin d'etre coherent avec le filtre (a Fe)
ssup = upsample(ss,Fse); % entre chaque on met "Fse" zero

% Filtrage
g1 = -0.5*ones(1,Fse/2);
g2 = 0.5*ones(1,Fse/2);
g= (1/sqrt(Fse)).*[g1 g2]; % Fonction porte

% Question 11 : Affichage de sl(t) pour les 25 premiers bits
sl= 0.5 + conv(ssup,g); % Convolution avec le filtre

ga=fliplr(g); %Filtre de reception 

% Matrices receptrices
variance = zeros(1, length(Eb_No));
TEB_real = zeros(1, length(Eb_No));
TEB_real_error = zeros(1, length(Eb_No));
EbNo_l = ones(1, length(Eb_No));

%% Synchonisation - Desynchronisation
delta_t = randi(100); % Delai de propagation temporelle
delta_f = randi([-1*10.^3 1*10.^3]); % Delai de propagation frequentielle
desyncro_t = zeros(1,delta_t); % Desynchro temporelle
desyncro_f = zeros(1,delta_f); % Desynchro frequentielle
sp = [1 1 0 0 0 0 0 0]; % Preambule
correlation = zeros(delta_t,delta_f);

sl_det = [sp, sl]; % Signal sl deterministe pour t[0,Tp]
t=(-delta_t+1:length(sl_det)); % t prenant en compte le decalage
sl_desynchro = [desyncro_t sl_det].*exp(-i*2*pi*delta_f*Te.*t); % Signal sl(t) desynchronise

% SNR de 1 a 10 -> TEB pratique
for i=1:1:10
    EbNo_l(i) = 10.^(i/10);
    variance = 1/(2*EbNo_l(i)); % Bruit
    sigma = sqrt(variance);
    yl = sl + (sigma * randn(1,length(sl))); % Signal bruite

    yl_desynchro = sl_desynchro; % Signal bruite desynchronise
    
     % Re-Synchronisation signal par estimations des parametres
   correlation_part1 = sum(yl_desynchro); % Numerateur
   correlation_part2 = sqrt(sum(abs(sp).^2)).*sqrt(sum(abs(yl_desynchro).^2)); % Denominateur
   correlation = correlation_part1/correlation_part2;

       % Esimation des parametres pour la synchronisation    
    [max_f indice_f]=max(abs(correlation));
    [max_t indice_t]=max(max_f);

    delta_t_estime= delta_t(indice_t) % Delai de propagation temporelle estime
    delta_f_estime= delta_f(indice_f) % Delai de propagation frequentielle estime
    
    rl= conv(yl_desynchro, g);
    
    rln = downsample(rl, Fse);
    An = rln;
    An(An >= 0) = 1;
    An(An < 0) = -1;
    b = An;
    b(b == -1) = 0;
    b(b == 1) = 1;
    b = b(2:end-3);
    TEB_real(i) = tauxDerreurBinaire(sb, b);  % Fonction qui caclul le TEB
end

%% Afficahge du TEB pour trames desynchronisees
Pb=0.5*erfc(sqrt(Eb_No)); % TEB theorique

%% Afficahge du TEB 
figure(1),
semilogy(10*log10(Eb_No),Pb,'r') ;
hold on
semilogy(10*log10(Eb_No),TEB_real,'g');
title('TEB en fonction du SNR en db');
legend('TEB theorique','TEB pratique');
hold off
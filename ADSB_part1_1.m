%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 Benjamin BONNOTTE & Abdoul GUISSET - G2                 %
%                     Projet Communication Numerique                      %
%                                 --------                                %
%                SIMULATION D'UN EMETTEUR / RECEPTEUR ADS-B               %
%                              ADSB_part1_1                               %
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
% t=0:Te:50*Ts-Te;
f = (0:Nfft-1)/Nfft-0.5; % Axe en frequence 
M=2; % Nombre de symboles dans la constellation

n = 1000; % Nombre de bits du message

M=2; % Nombre de symboles dans la constellation
Npaq = 112; % Nombres de bits par paquets

%% Emetteur
% Bits aleatoires
sb = randi([0 1],1,Npaq); % Generation d'une sequence aleatoire de 1000 bits

% figure(1), plot(sb);
% title('Emission de 112 bits');

% Association bits -> Symboles
ss = pammod(sb,M);

% Sur-echantillonage afin d'etre coherent avec le filtre (a Fe)
ssup = upsample(ss,Fse); % entre chaque on met "Fse" zero

% Filtraage

g1 = -0.5*ones(1,Fse/2); % Filtre 1
g2 = 0.5*ones(1,Fse/2); % Filtre 2

g= (1/sqrt(Fse)).*[g1 g2]; % Filtre

% Question 11 : Affichage de sl(t) pour les 25 premiers bits
sl= 0.5 + conv(ssup,g); % Convolution avec le filtre
% figure(2), plot(sl(1:25*Fse)); % 25 premiers bits
% xlabel('time');
% ylabel('sl(t)');
% title('Allure temporelle du signal sl(t) pour les 25 premirs bits');

% Question 12 : Diagramme de l'oeil de sl(t)
% eyediagram(sl(1:100*Fse),2*Fse,2*Ts);

% Question 13 : DSP de sl(t)
% figure(4),
% subplot(2,1,1);
% pwelch(sl,ones(1,Nfft),0,Nfft,Fe,'centered')
% xlabel('Normalized frequency')
% ylabel('Magnitude')
% title('DSP pratique de sl(t)')
% 
% % DSP theorique 
theoric_DSP=(1/(4*Ts))*(abs(fft(g, Nfft)).^2+abs(fft(g, Nfft)).^2); %DSP theorique
% 
% subplot(2,1,2);
% plot(f,fftshift(theoric_DSP))
% xlabel('Normalized frequency')
% ylabel('Magnitude')
% title('DSP théorique de sl(t)')
experimental_DSP=fftshift(abs(fft(sl, Nfft)).^2); %DSP pratique

figure(4),
semilogy(f,fftshift(theoric_DSP),'g');
hold on
semilogy(f,experimental_DSP,'r') ;
title('DSP');
legend('DSP Theorique','DSP Pratique');
hold off
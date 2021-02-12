%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SIG - Initiation Traitment Signal : TP2 Extraction ECG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all; 
clc;

%% 1/
% Chargement des donn�es
load(..);
load(..);

%% 2/
% Initialisation : les deux fichiers ont exactement la meme duree, donc pas
% besoin de creer des variables separees.
fe = 1000;                      % Frequence d'echantillonnage en Hz
N = ..;          % Taille du vecteur
t = linspace(0, N/fe, N);       % Echelle temporelle

% Affichage des series temporelles
figure(1); 
subplot(2, 1, 1); plot(.., ..);
xlabel('Temps (s)'); ylabel('Amplitude (mV)'); title('Signal ECG au repos')
subplot(2, 1, 2); plot(.., ..);
xlabel('Temps (s)'); ylabel('Amplitude (mV)'); title('Signal ECG pendant effort')

%% 3/
NFFT = 2 .^ nextpow2(N);                   % Next power of 2 from length of x
FFT_ECG_Rep = fft(.., ..);        % TF repos 
FFT_ECG_Eff = fft(.., ..);        % TF effort
f = fe/2 * linspace(0, 1, NFFT/2);         % Semie-echelle fr�quentielle

% Affichage des spectres (module uniquement, pas la phase) 
figure(2);
subplot(2, 1, 1); plot(f, 2*abs(FFT_ECG_Rep(1:NFFT/2)));
xlabel('Frequency (Hz)'); ylabel('|FFTECGRep(f)|'); 
title('Single-Sided Amplitude Spectrum of XECGRep');
subplot(2, 1, 2); plot(f, 2*abs(FFT_ECG_Eff(1:NFFT/2)));
xlabel('Frequency (Hz)'); ylabel('|FFTECGEff(f)|'); 
title('Single-Sided Amplitude Spectrum of XECGEff');

%% 4/ Extraction de la s�rie temporelle RR
%% Etape 1 : Seuillage
seuil_Rep = ..;
seuil_Eff = ..;

figure(3);
% A vous de jouer


%% Etape 2 :
% Vecteur logique repos
for i = ..  
    if ..
        Diag_ECG_Rep(1, i) = 1;
    else
        Diag_ECG_Rep(1, i) = 0;
    end
end

% Vecteur logique effort
for i = ..  
    if ..
        Diag_ECG_Eff(1, i) = 1;
    else
        Diag_ECG_Eff(1, i) = 0;
    end
end

figure(4);
% Affichage vecteur logique

%% Etapes 3 :
% Repos
Diff_Diag_ECG_Rep_P = find(diff(Diag_ECG_Rep) == 1); Diff_Diag_ECG_Rep_P(end) = [];
Diff_Diag_ECG_Rep_N = find(diff(Diag_ECG_Rep) == -1);
A_Rep = [];
B_Rep = [];
for k = 1 : length(Diff_Diag_ECG_Rep_P)
    [A_Rep(k), B_Rep(k)] = max(X_ECG_Rep(Diff_Diag_ECG_Rep_P(1, k) : Diff_Diag_ECG_Rep_N(1, k)));
    B_Rep(k) = B_Rep(k) + Diff_Diag_ECG_Rep_P(1, k);
end
Peak_ECG_Rep = zeros(size(X_ECG_Rep));
Peak_ECG_Rep(B_Rep) = X_ECG_Rep(B_Rep);

% Effort
Diff_Diag_ECG_Eff_P = find(diff(Diag_ECG_Eff) == 1); Diff_Diag_ECG_Eff_P(end) = [];
Diff_Diag_ECG_Eff_N = find(diff(Diag_ECG_Eff) == -1);
A_Eff = [];
B_Eff = [];
for k = 1 : length(Diff_Diag_ECG_Eff_P)
    [A_Eff(k), B_Eff(k)] = max(X_ECG_Eff(Diff_Diag_ECG_Eff_P(1, k) : Diff_Diag_ECG_Eff_N(1, k)));
    B_Eff(k) = B_Eff(k) + Diff_Diag_ECG_Eff_P(1, k);
end
Peak_ECG_Eff = zeros(size(X_ECG_Eff));
Peak_ECG_Eff(B_Eff) = X_ECG_Eff(B_Eff);

figure(1);  % Reprise de la figure 1
subplot(2, 1, 1); hold on; plot(t(B_Rep), Peak_ECG_Rep(B_Rep), 'r*');
subplot(2, 1, 2); hold on; plot(t(B_Eff), Peak_ECG_Eff(B_Eff), 'r*');

%% Etape 4 :
RR_ECG_Rep = diff(t(..));
RR_ECG_Eff = diff(t(..));

figure(5);
% Affichage series temporelles RR

%% 5/
fc1 = ..;    % Fr�quence de coupure haute
fc2 = ..;   % Fr�quence de coupure basse
fe_RR = 1;    %Freq d'echantillonnage theorique de la s�rie RR

%Filtarge des series RR entre 0.04 Hz ET 0.4 Hz
RR_ECG_Rep_filtre = bandf(RR_ECG_Rep, (fc1*2)/fe_RR, (fc2*2)/fe_RR);
RR_ECG_Eff_filtre = bandf(RR_ECG_Eff, (fc1*2)/fe_RR, (fc2*2)/fe_RR);

% TF des signaux
nb_points_fft = ..;
f = fe_RR/2 * linspace(0, 1, nb_points_fft/2);
TF_RR_ECG_Rep_filtre = fft(.., ..);
TF_RR_ECG_Eff_filtre = fft(.., ..);

figure(6);
% Affichage demi-spectre series RR

%% 6/
PSD_RR_ECG_Rep_filtre = ..;
PSD_RR_ECG_Eff_filtre = ..;

figure(7);
% Affichage PSD des series RR filtrees

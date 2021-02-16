% Course: MAE 295 GNSS Signal Processing and SDR Design 
% Date: December 14th, 2020
% Author: Alex Nguyen
% Description: Updated Project, Save: C/N0, sigIQ^2, and PRN Sequence Data

clc; clear; close all

%% Load Saved Data:
load('C_N0dBVec.mat')       % C/N0 Data & sigIQ Squared
load('PRNdata.mat')         % Each SV PRN Sequences

%% Setup Parameters
% GPS L1 C/A Values:
fL1 = 1575.42e6;            % GPS L1 C/A Frequency [Hz] 
Rc = 1.02300325e6;          % Chip Rate [chips/sec]
Nc = 1023;                  % Number of Chips per Subaccum
Tc = 1e-3/Nc;               % Chip Sampling Period [s]
fc = 1/Tc;                  % Chip Frequency [Hz]
fsampIQ = 2500/1e-3;        % IQ sampling frequency [Hz]
teml = 0.8*Tc;              % Early Minus Late [Chips]

% Subaccum:
Tsub = 1e-3;                % Subaccum Sampling Interval
SubAccum = 1;               % Subaccum for L1 C/A (e.g. M = 1)
N = floor(fsampIQ*Tsub);    % # of Subaccums per Subaccum

% Load Data:
Tfull = 60;                                    % Time Interval of Data to Load [s]
tDataV = Tfull/Tsub;                           % Full Time Data [s] 
Ylen = floor(Tfull*fsampIQ);
fid = fopen('mystery_data_file2.bin','r','l'); % Open File
Y = fread(fid, [2, Ylen], 'int16')';           % Read Binary Data File
Y = Y(:,1) + 1i*Y(:,2);
fclose(fid);                                   % Close File 

%% PRN Identifier
% LSFR Parameter Setup:
n = 10; 
m = 2^n - 1;                    % Max-Length Period
a0Vec = ones(1, 10);            % Initial State
ciVec1 = [3, 10]';              % Characterisitic Polynomial
ciVec2 = [2, 3, 6, 8, 9, 10]';  
          
% Max-Length LSFR Sequences:
f1_D = generateLfsrSequence(n, ciVec1, a0Vec); 
f2_D = generateLfsrSequence(n, ciVec2, a0Vec); 

% SV Code Delay Chips:
SV = [5 6 7 8 17 18 139 140 141 251 252 254 255 256 257 258 469 470 471 ...
      472 473 474 509 512 513 514 515 516 859 860 861 862 863 950 947 ...
      948 950]';
  
% Preallocation:
Gpm_os = zeros(length(SV), N);
  
for sv = 1:length(SV)
      % GPS PRN Gold Code:
      G = mod(circshift(f2_D, SV(sv)) + f1_D, 2);
      Gpm = 2*G - 1;
      
      % Oversampled GPS PRN Gold Sequence:
      Gpm_os(sv, :) = oversampleSpreadingCode(Gpm, Rc/fsampIQ, N, m);
end

% Save PRN Data:
save('PRNdata.mat', 'Gpm_os', 'f1_D', 'f2_D');

%% Signal Acquisition
% Doppler Frequency and Code Start Time Mesh:
fLow = -6000; fHigh = 6000; df = 100; f = fLow:df:fHigh;
tLow = 0; tHigh = Tsub - 1/fsampIQ; t = linspace(tLow, tHigh, fsampIQ*Tsub)';

% SV Code Delay Chips:
SV = [5 6 7 8 17 18 139 140 141 251 252 254 255 256 257 258 469 470 471 ...
      472 473 474 509 512 513 514 515 516 859 860 861 862 863 950 947 ...
      948 950]';

% Preallocation:
Sktilde_sq = cell(length(SV), 1);
C_N0dB = zeros(length(SV),1); 

for sv = 1:length(SV)
    % Initialize:
    Sk_sq = 0;
    
    % GPS PRN Gold Code:
    G = mod(circshift(f2_D, SV(sv)) + f1_D, 2);
    Gpm = 2*G - 1;
    
    % Oversampled GPS PRN Gold Sequence:
    Gpm_os = oversampleSpreadingCode(Gpm, Rc/fsampIQ, N, m);
    
    for M = 1:SubAccum
        % Initialize:
        k0 =  1; kend = N;
        Zk = zeros(length(f), length(t));
        
        for ff = 1:length(f)
            % Complex Mixing of Complex Sequence:
            xk = Y(k0:kend).*exp(-1i*2*pi*f(ff)*t);
            
            % Discrete Fourier Transform:
            Xr = fft(xk);       % Mixed Complex Signal
            Cr = fft(Gpm_os);   % PRN Sequence
            
            % Circular Correlation Sequence
            Zr = Xr.*conj(Cr);  
            zk = ifft(Zr);
            
            if sv == 29
                % Save Peak Sk Value:
                Sk_sq29 = zk;
                Sk0_sq29 = max(max(zk));
            end

            % |Sk_tilde|.^2 Per Subaccum:
            Zk(ff, :) = zk.*conj(zk);  
            
        end
        % Sum Previous Subaccums:
        Sk_sq = Sk_sq + Zk;
        
        % Update Complex Signal Indices:
        k0 = k0 + N; kend = kend + N;
 
    end
    % Average GPS L1 C/A Signal Value:
    Sk_sq = Sk_sq/SubAccum;
    Sktilde_sq{sv} = Sk_sq;  % Save 
       
    % GPS L1 C/A Signal Max Value:
    [Pk_t, loc_t] = max(max(Sk_sq));   
    [Pk_f, loc_f] = max(max(Sk_sq'));
    
    % Carrier-to-Noise Ratio [dB-Hz]:
    C_N0 = (Pk_t - 2*sigIQ_sq)/2/sigIQ_sq/Tsub;  
    C_N0dB(sv) = 10*log10(C_N0);       

    % Determine Sigma_IQ:
    if sv == 37
        % No Signal Present @ sv = 37:
        sigIQ_sq = mean(mean(Sk_sq))/2;
    end
end

%Save Data:
save('C_N0dBVec.mat', 'C_N0dB', 'sigIQ_sq')


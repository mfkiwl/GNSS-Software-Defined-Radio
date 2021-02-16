% Course: MAE 295 GNSS Signal Processing and SDR Design 
% Date: December 14th, 2020
% Author: Alex Nguyen
% Description: Updated Project Part 1, Top-Level Acquisition and Tracking
% Script for PRN 32

clc; clear; close all

%% Load Saved Data:
load('C_N0dBVec.mat')       % C/N0 Data & sigIQ Squared
load('PRNdata.mat')         % Each SV PRN Sequences

%% Setup Parameters
% Speed of Light [m/s]:
c = 299792458;

% GPS L1 C/A Values:
fL1 = 1575.42e6;            % GPS L1 C/A Frequency [Hz] 
Rc = 1.02300325e6;          % Chip Rate [chips/sec]
Nc = 1023;                  % Number of Chips per Subaccum
Tc = 1e-3/Nc;               % Chip Sampling Period [s]
fc = 1/Tc;                  % Chip Frequency [Hz]
fsampIQ = 2500/1e-3;        % IQ sampling frequency [Hz]
teml = 0.8*Tc;              % Early Minus Late [Chips]

% Subaccum:
Tsub = 1e-3;                            % Subaccum Sampling Interval
N = floor(fsampIQ*Tsub);                % # of Subaccums per Subaccum
t = linspace(0, Tsub - 1/fsampIQ, N)';  % Subaccum Time Interval [s]

% Load Data:
Tfull = 60;                                    % Time Interval of Data to Load [s]
tDataV = Tfull/Tsub;                           % Full Time Data [s] 
Ylen = floor(Tfull*fsampIQ);
fid = fopen('mystery_data_file2.bin','r','l'); % Open File
Y = fread(fid, [2, Ylen], 'int16')';           % Read Binary Data File
Y = Y(:,1) + 1i*Y(:,2);
fclose(fid);                                   % Close File                        

%% (a) Acquire a PRN Sequence
% Initialize:
sv = 32;                               % PRN Signal to Acquire                             
SubAccum = 1;                          % Subaccum # for GPS L1 C/A 

% Perform Acquisition:
[ts_hat, fDk, C_N0dB, sig, Sk_sq] = acquireSV(sv, SubAccum); 

%% (b) Initialize Beat Carrier Phase Estimate
% Carrier Phase:
th_beat = 0;                     % (Arbitary) Initial Beat Carrier Phase                         

%% (c) Initiailze Moving-Window Average
% Initial In-Phase and Quadrature Peak Value:
Ip0_sq = real(Sk_sq); Qp0_sq = imag(Sk_sq);
SkAvg = Sk_sq;

%% (d) Initialize Phase Tracking Loop Filter 
% Initialize:
BL_target = 10;  % Target Bandwidth [Hz]
loopOrder = 3;   % Loop Order
BL_codeT = 0.2;  % Code Tracking Loop [Hz]  

% Phase Tracking Loop Filter:
[Ad, Bd, Cd, Dd, BL_act] = configureLoopFilter(BL_target, Tsub, loopOrder);

% Initial State Matrix:
M = [Ad(1, 1) - 1, Ad(1, 2); ...
       Cd(1)     ,   Cd(2)];
 
b = [0; 2*pi*fDk];

% Initialize States:
xk = M\b;

%% (f) Tracking Signal
% Preallocate:
ts_hatV = zeros(tDataV, 1);
th_beatV = ts_hatV;
fDV = ts_hatV;
IVec = ts_hatV; QVec = ts_hatV;
ek_p = ts_hatV; ek_d = ts_hatV;
C_N0dBV = ts_hatV;
vTotalkV = ts_hatV;

% Initialize:
jk = 1; jkend = N;       % Input Accumulation Indicies
s.BL_target = BL_codeT;  % Dll Target Bandwidth
s.IsqQsqAvg = 0;         % Average Sk Value

% Gold Code SV Sequence for Subaccum:
Gpm = Gpm_os(sv, :)';  

for ii = 1:tDataV
% Carrier Phase Estimate:
th_hat = 2*pi*fDk*t + th_beat;   

% Correlator Block:
[Spk, Sek, Slk] = ...
    CorrelatorBlk(Y(jk:jkend), Gpm, fsampIQ, ts_hat, th_hat, teml);

% Update Moving-Window Average:
SkAvg = Spk.*conj(Spk);
C_N0 = (SkAvg - 2*sigIQ_sq)/2/sigIQ_sq/Tsub; 
C_N0dB = 10*log10(C_N0);

% Initialize updatePll.m Structure:
s.Ipk = real(Spk); s.Qpk = imag(Spk);
s.xk = xk;
s.Ad = Ad; s.Bd = Bd;
s.Cd = Cd; s.Dd = Dd;

% Phase Tracking Loop:
[xkp1, vk, ek_pll] = updatePll(s);

% Initialize updateDll.m Structure:
s.IsqQsqAvg = (s.IsqQsqAvg*ii + SkAvg)/ii;
s.sigmaIQ = sqrt(sigIQ_sq);
s.Ipk = real(Spk); s.Qpk = imag(Spk); 
s.Iek = real(Sek); s.Qek = imag(Sek);
s.Ilk = real(Slk); s.Qlk = imag(Slk);
s.vpk = vk/2/pi/fL1;
s.Tc = Tc;

% Carrier-Aided Tracking Loop:
[vTotalk, ek_dll] = updateDll(s);

% Save:
fDV(ii) = fDk;
ts_hatV(ii) = ts_hat;
th_beatV(ii) = th_beat;
IVec(ii) = real(Spk); 
QVec(ii) = imag(Spk); 
ek_d(ii) = ek_dll; ek_p(ii) = ek_pll;
C_N0dBV(ii) = C_N0dB;
vTotalkV(ii) = vTotalk;

% Update:
xk = xkp1;                                % Loop Filter State
fDk = vk/2/pi;                            % Approximate Doppler   
ts_hat = ts_hat + (Nc*Tc)/(1 + vTotalk);  % Code Start Time
th_beat = th_beat + vk*Tsub;              % Beat Carrier Phase Estimate  
jk = jk + N; jkend = jkend + N;           % Complex Subaccum Indices
end

%% Estimate of j0 Index
j0 = ts_hatV(1)*fsampIQ;

% Print Results:
fprintf('Part (ii): \n')
fprintf('Sample # After First Complete C/A Code Interval: '); disp(vpa(j0));

%% Plot Results
t = linspace(0, Tfull, length(ts_hatV));

% (a) Code Error [chips]
figure(1); 
plot(t, ek_d/Tc, 'Color', '#D95319', 'linewidth', 1)
xlabel('Time [s]');
ylabel('$e_k$', 'interpreter', 'latex');
legend('Arctangent (AT)')
title(sprintf('(a) SV # %d: Code Error [Chips]', sv)) 

% (b) Phase Error [deg]
figure(2);
plot(t, ek_p*360/2/pi, 'Color', '#A2142F', 'linewidth', 1)
xlabel('Time [s]');
ylabel('$e_k$', 'interpreter', 'latex');
legend('Arctangent (AT) ')
title(sprintf('(b) SV # %d: Phase Tracking Loop Error [deg]', sv)) 

% (c) Estimated Doppler Frequency [Hz]
figure(3);
plot(t, fDV, 'Color', '#77AC30', 'linewidth', 1.5);
xlabel('Time [s]');
ylabel('$fD_k$', 'interpreter', 'latex');
title(sprintf('(c) SV # %d: Estimated Doppler Frequency [Hz]', sv));

% (d) Estimated Carrier-to-Noise Ratio [dB-Hz]
figure(4);
plot(t, real(C_N0dBV), '.' , 'Color', '#4DBEEE', 'linewidth', 1)
xlabel('Time [s]');
ylabel('$C/N0$', 'interpreter', 'latex');
title(sprintf('(d) SV # %d: Estimated Carrier-to-Noise Ratio [dB-Hz]', sv)) 

% (e) Tracked In-Phase and Quadrature Components
figure(5);
plot(t, IVec, 'Color', '#EDB120', 'linewidth', 1.5); hold on; 
plot(t, QVec, 'Color', '#7E2F8E', 'linewidth', 1.5);
xlabel('Time [s]');
legend('$\tilde{I}_k$', '$\tilde{Q}_k$', 'interpreter', 'latex', 'location', 'best'); 
title(sprintf('(e) SV # %d: Tracked In-Phase and Quadrature Subaccum', sv));

% (f) Code Phase Estimate [m]
figure(6)
plot(t, ts_hatV*c, 'Color', '#77AC30', 'linewidth', 1.5);
xlabel('Time [s]');
title(sprintf('(f) SV # %d: Code Phase Estimate [m]', sv));

% (g) Phase Portrait of Sk AFTER Lock
iiPostLock = 2500;  % Sample index achieving PLL Lock

figure(7)
plot(IVec(iiPostLock:end), QVec(iiPostLock:end), '.');
xlabel('$\tilde{I}_k$', 'interpreter', 'latex');
ylabel('$\tilde{Q}_k$', 'interpreter', 'latex');
title(sprintf('(g) SV # %d: Sk Complex Plane after Lock', sv));

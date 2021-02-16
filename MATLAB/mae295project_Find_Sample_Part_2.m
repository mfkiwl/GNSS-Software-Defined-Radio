% Course: MAE 295 GNSS Signal Processing and SDR Design 
% Date: December 14th, 2020
% Author: Alex Nguyen
% Description: Updated Project Top-Level Acquisition and Tracking Script
% for PRN 10, 11, 14, 31, and 32

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
sv = 10;                               % PRN Signal to Acquire                             
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
tsV = zeros(tDataV, 1);
fDV = tsV;
IVec = tsV;
QVec = tsV;
C_N0dBV = tsV;
ek_p = tsV;
ek_d = tsV;
vTotalkV = tsV;

% Initialize Input Signal Indicies:
jk = 1; jkend = N;

for ii = 1:tDataV
% Carrier Phase:
th_hat = 2*pi*fDk*t + th_beat;   

% Perform Correlation Over Subaccum:
Gpm = Gpm_os(sv, :)';  % PRN Code

[Spk, Sek, Slk] = CorrelatorBlk(Y(jk:jkend), Gpm, fsampIQ, ts_hat, th_hat, teml);

% Update Moving-Window Average:
SkAvg = abs(Spk).^2;
C_N0 = (SkAvg - 2*sigIQ_sq)/2/sigIQ_sq/Tsub; 
C_N0dBV(ii) = 10*log10(C_N0);

% Initialize updatePll.m Structure:
s.Ipk = real(Spk); 
s.Qpk = imag(Spk);
s.xk = xk;
s.Ad = Ad;
s.Bd = Bd;
s.Cd = Cd;
s.Dd = Dd;

% Phase Tracking Loop:
[xkp1, vk, ek_pll] = updatePll(s);

% Initialize updateDll.m Structure:
s.BL_target = BL_codeT;
s.IsqQsqAvg = SkAvg;
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
tsV(ii) = ts_hat;
IVec(ii) = real(Spk); QVec(ii) = imag(Spk); 
ek_d(ii) = ek_dll; ek_p(ii) = ek_pll;
vTotalkV(ii) = vTotalk;

% Update:
xk = xkp1;
ts_hat = ts_hat + (Nc*Tc)/(1 + vTotalk);  % Code Start Time
fDk = vk/2/pi;                            % Approximate Doppler   
th_beat = th_beat + vk*Tsub;              % Beat Carrier Phase Constant 
jk = jk + N; jkend = jkend + N;           % Input Signal Indices
end

% Incase Zero Array Values:
tsV(tsV == 0) = [];
fDV(fDV == 0) = [];
th_hat(th_hat == 0) = [];
IVec(IVec == 0) = []; QVec(QVec == 0) = [];
SkAvg(SkAvg == 0) = [];
C_N0dBV(C_N0dBV == 0) = [];

%% Find Sample Rx
% tRxRaw (Samples):
samp = 135000600; 

% Calculate True Samples:
tr_samp = samp/fsampIQ;
ind = tsV(tsV >= tr_samp);
tr_samp_true = vpa(ind(1)*fsampIQ);

% Print Results:
fprintf('tRxRaw (Samples) = %d or ', ceil(tr_samp_true)); disp(tr_samp_true);

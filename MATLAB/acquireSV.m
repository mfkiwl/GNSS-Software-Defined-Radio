function [ts0, fD0, C_N0dB, sig, Pk_t] = acquireSV(sv, SubAccum) 
% acquireSV: Acquires the SV signal then outputs the inital Doppler, 
%            intial code start time, Carrier-to-Noise Ratio, determine if a 
%            signal is present, and Initial peak signal power used in moving
%            window average.
%
% Inputs
%       sv: PRN #
% SubAccum: # of Coherent Sumsubaccums 
%
% Outputs
%      ts0: Initial Code Start Time [s]
%      fD0: Initial Doppler Estimate [Hz]
%   C_N0dB: Carrier-to-Noise Ratio [dB]
%      sig: (1) Signal Present or (0) No Signal Present
%
%-------------------------------------------------------------------------%
%% Load Data (Previously Computed):
load('PRNdata.mat')
load('C_N0dBVec.mat')

%% Initialize Parameters:
% GPS L1 C/A Values:
Rc = 1.02300325e6;          % Chip Rate [chips/sec]
Tc = 1e-3/1023;             % Chip Sampling Period [s]
fsampIQ = 2500/1e-3;        % IQ sampling frequency [Hz]
Tfull = 5;                  % Time Interval of Data to Load [s]
Ylen = floor(Tfull*fsampIQ);

% Subaccum:
Tsub = 1e-3;                % Subaccum Sampling Interval
N = floor(fsampIQ*Tsub);    % # of Subaccums per Subaccum

% Load Data:
fid = fopen('mystery_data_file2.bin','r','l');  % Open File
Y = fread(fid, [2, Ylen], 'int16')';           % Read Binary Data File
Y = Y(:,1) + 1i*Y(:,2);
fclose(fid);                                   % Close File

%% Signal Acquisition:
% Doppler Frequency Time Mesh:
fLow = -6000; fHigh = 6000; df = 150; f = fLow:df:fHigh;

% Code Start Time Mesh:
tLow = 0; tHigh = Tsub - 1/fsampIQ; t = linspace(tLow, tHigh, N)';

% SV Code Delay Chips:
SV = [5 6 7 8 17 18 139 140 141 251 252 254 255 256 257 258 469 470 471 ...
      472 473 474 509 512 513 514 515 516 859 860 861 862 863 950 947 ...
      948 950]';
  
% Carrier-to-Noise Ratio for NO Signal Present:
C_N0dB_avg = mean(C_N0dB(C_N0dB > 40)); 
clear C_N0dB

% Initialize:
Sk_sq = 0;
Gpm = Gpm_os(sv, :)';

for M = 1:SubAccum
    % Initialize:
    k0 = 1; kend = N;
    Zk = zeros(length(f), length(t));
    
    for ff = 1:length(f)
        % Complex Mixing of Complex Sequence:
        xk = Y(k0:kend).*exp(-1i*2*pi*f(ff)*t);
        
        % Discrete Fourier Transform:
        Xr = fft(xk);       % Mixed Complex Signal
        Cr = fft(Gpm);      % PRN Sequence
        
        % Circular Correlation Sequence
        Zr = Xr.*conj(Cr);
        zk = ifft(Zr);
        
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

% GPS L1 C/A Signal Max Value:
[Pk_t, loc_t] = max(max(Sk_sq));
[~, loc_f] = max(max(Sk_sq'));

% Carrier-to-Noise Ratio [dB-Hz]:
C_N0 = (Pk_t - 2*sigIQ_sq)/2/sigIQ_sq/Tsub;
C_N0dB = 10*log10(C_N0);

% Determine if GPS L1 C/A Signal is Present:
if C_N0dB > C_N0dB_avg
    % Approximate Values:
    ts0 = t(loc_t);     % Code Start Time Offset [s]
    fD0 = f(loc_f);     % Apparent Doppler Frequency [Hz]
    
    % Determine if Signal Present:
    sig = 1;
   
else
    sig = 0;
end

end
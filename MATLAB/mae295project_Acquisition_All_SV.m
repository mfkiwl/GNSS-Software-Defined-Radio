% Course: MAE 295 GNSS Signal Processing and SDR Design 
% Date: December 14th, 2020
% Author: Alex Nguyen
% Description: Updated Project, Acquisition of All Signals 

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


%% Signal Acquisition
% Doppler Frequency and Code Start Time Mesh:
fLow = -6000; fHigh = 6000; df = 150; f = fLow:df:fHigh;
tLow = 0; tHigh = Tsub - 1/fsampIQ; t = linspace(tLow, tHigh, fsampIQ*Tsub)';

% SV Code Delay Chips:
SV = [5 6 7 8 17 18 139 140 141 251 252 254 255 256 257 258 469 470 471 ...
      472 473 474 509 512 513 514 515 516 859 860 861 862 863 950 947 ...
      948 950]';
  
% Carrier-to-Noise Ratio for NO Signal Present:
C_N0dB_avg = mean(C_N0dB(C_N0dB > 40)); 
clear C_N0dB

% Preallocation:
Sktilde_sq = cell(length(SV), 1);

for sv = 1:length(SV)
    % Initialize:
    Sk_sq = 0;

    % Oversampled GPS PRN Gold Sequence:
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
            Cr = fft(Gpm);   % PRN Sequence
            
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
    Sktilde_sq{sv} = Sk_sq;  % Save 
       
    % GPS L1 C/A Signal Max Value:
    [Pk_t, loc_t] = max(max(Sk_sq));   
    [Pk_f, loc_f] = max(max(Sk_sq'));
    
    % Carrier-to-Noise Ratio [dB-Hz]:
    C_N0 = (Pk_t - 2*sigIQ_sq)/2/sigIQ_sq/Tsub;  
    C_N0dB = 10*log10(C_N0);           
    
    % Determine if GPS L1 C/A Signal is Present
    if C_N0dB > C_N0dB_avg
        % Plot Acquired Signals:
        figure(5 + sv);
        surf(t, f, Sk_sq, 'FaceColor', 'interp')
        title(sprintf('SV # %d', sv));
        xlabel('Code Start Time [s]'); ylabel('Frequency [Hz]');
        
        % Approximate Values:
        ts = t(loc_t);     % Code Start Time Offset [s]
        fd = f(loc_f);     % Apparent Doppler Frequency [Hz]
        
        % Print Results:
        fprintf(['SV # %d Signal Acquired: fD = %d Hz, ts = %4.4e sec, and '...
            'C/N0 = %4.4f dB-Hz \n'], sv, fd, ts, C_N0dB); 
    end       
end

%%  Plot ALL GPS L1 C/A Signals:
subp = 0;     % (0) No Plot, (1) Plot

if subp == 1
% Initialize:
ii = 1; jj = 1; sv = 1;

while sv <= 37
    % Reset Counter: 
    jj = 1;
    
    while jj <= 9        
        % Plot Signal Acquisition:
        figure(ii)
        subplot(3, 3, jj)
        surf(t, f, Sktilde_sq{sv}, 'FaceColor', 'interp');
        title(sprintf('SV # %d', sv));
        
        % Update Counter:
        jj = jj + 1;
        sv = 1 + sv;
        
        % Last SV:
        if sv == 38 
            break
        end        
    end
    
    % Update Figure Parameters:
    ii = ii + 1;
    
    % Subplot Title:
    if sv ~= 38
        sgtitle(sprintf('Signal Acquistion for PRN # %d to %d', sv - 9, sv - 1));
    elseif sv == 38
        sgtitle(sprintf('Signal Acquistion for PRN # %d ', sv - 1));
    end
end
end
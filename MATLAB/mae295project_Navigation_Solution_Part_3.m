% Course: MAE 295 GNSS Signal Processing and SDR Design 
% Date: December 17th, 2020
% Author: Alex Nguyen
% Description: Project Navigation Solution.

clc; clear; close all;

%% Load Saved Data:
load('C_N0dBVec.mat')       % C/N0 Data & sigIQ Squared
load('PRNdata.mat')         % Each SV PRN Sequences

%% Setup 
% Initialize Constants:
N = 5;                       % Measurement #
c = 299792458;               % Speed of Light [m/s]
fsampIQ = 2.5e6;             % Sampling Frequency [Hz]
Ts = 1/fsampIQ;              % Uniform Sample Spacing [sec/sample] 
off = 442681.598;            % Constant Reciever Time Offset

% Load Data:
X = load('RxandSVClk.txt');
PRN = X(:, 1);               % PRN #
tr_samp = X(:, 2)*Ts;        % Reciever Time [s]
tr  = tr_samp + off;         % Reciever Time with Offset [s]
ts = X(:, 3);                % SV Time [s]

% Compute Pseudoranges:
z = c*(tr - ts);

%% Pseduorange Variance:
Nc = 1023;          % Number of Chips per Subaccum
Tc = 1e-3/Nc;       % Chip Sampling Period [s]
teml = 0.8*Tc;      % Early Minus Late [Chips]
d = teml/Tc;
Bdll = 0.2; 
C_N0 = 10.^(C_N0dB./10);
CN0 = [C_N0(10); C_N0(11); C_N0(14); C_N0(31); C_N0(32)];
sig_delt = c*sqrt((d*Bdll*Tc.^2)./(2.*CN0));  % Standard Deviation

% Pseudorange Measurment Noise: 
w = [sig_delt(1), sig_delt(2), sig_delt(3), sig_delt(4), sig_delt(5)]';

%% Initialize Pseudorange Equation Values
% Augmented Pseudorange Model:
delR = @(x, sv_pos) norm(x(1:3) - sv_pos); 

% PRN 10:
rSvECEF10 = [9121305.7077300, -23601900.3375000, 7954520.8506900]';
dtIono10 = 2.273761e-08;
dtTropo10 = 1.622817e-08;
dtSV10 = -0.000100688843623000;
z_10 = @(x) delR(x, rSvECEF10) + c*(x(4)/c - dtSV10) + c*(dtTropo10 + dtIono10);

% PRN 11:
rSvECEF11 = [-22226444.3115000, -7263521.9296100, 12021436.2202000]';
dtIono11 = 1.945468e-08;
dtTropo11 = 1.321500e-08;
dtSV11 = -0.000664972769672000;
z_11 = @(x) delR(x, rSvECEF11) + c*(x(4)/c - dtSV11) + c*(dtTropo11 + dtIono11);

% PRN 14:
rSvECEF14 = [-7443729.8892700, -15525794.7538000, 20526967.7486000]';
dtIono14 = 1.439808e-08;
dtTropo14 = 9.258317e-09;
dtSV14 = -0.000036453122726400;
z_14 = @(x) delR(x, rSvECEF14) + c*(x(4)/c - dtSV14) + c*(dtTropo14 + dtIono14);

% PRN 31:
rSvECEF31 = [-6485206.0043400, -24828766.7942000, 6135430.7428900]';
dtIono31 = 1.545366e-08;
dtTropo31 = 1.003285e-08;
dtSV31 = 0.000246467440852000;
z_31 = @(x) delR(x, rSvECEF31) + c*(x(4)/c - dtSV31) + c*(dtTropo31 + dtIono31);

% PRN 32
rSvECEF32 = [4189914.7336700, -14908405.1919000, 21590991.0269000]';
dtIono32 = 1.854755e-08;
dtTropo32 = 1.245384e-08;
dtSV32 = -0.000125803551001000;
z_32 = @(x) delR(x, rSvECEF32) + c*(x(4)/c - dtSV32) + c*(dtTropo32 + dtIono32);

% Measurment Equations:
h = @(x) [z_10(x); z_11(x); z_14(x); z_31(x); z_32(x)] + w*randn;

%% Nonlinear Least Squares Estimator 
% Initialize:
xhat0 = zeros(4, 1);          % Initial Reciever Guess
gamma = 1e-3;                 % Desired Criteria < 1 mm
delXMag = 1;                  % Initial Trajectory Error
count = 1;                    % counter

% Preallocation:
xhat_est = xhat0;

while delXMag > gamma
    % Redefine:
    xhat = xhat0;
    
    % Jacobian Matrix:
    H = @(x) [(x(1) - rSvECEF10(1))/delR(x, rSvECEF10), (x(2) - rSvECEF10(2))/delR(x, rSvECEF10), (x(3) - rSvECEF10(3))/delR(x, rSvECEF10), 1 ; ...
              (x(1) - rSvECEF11(1))/delR(x, rSvECEF11), (x(2) - rSvECEF11(2))/delR(x, rSvECEF11), (x(3) - rSvECEF11(3))/delR(x, rSvECEF11), 1 ; ...
              (x(1) - rSvECEF14(1))/delR(x, rSvECEF14), (x(2) - rSvECEF14(2))/delR(x, rSvECEF14), (x(3) - rSvECEF14(3))/delR(x, rSvECEF14), 1 ; ...
              (x(1) - rSvECEF31(1))/delR(x, rSvECEF31), (x(2) - rSvECEF31(2))/delR(x, rSvECEF31), (x(3) - rSvECEF31(3))/delR(x, rSvECEF31), 1 ; ...
              (x(1) - rSvECEF32(1))/delR(x, rSvECEF32), (x(2) - rSvECEF32(2))/delR(x, rSvECEF32), (x(3) - rSvECEF32(3))/delR(x, rSvECEF32), 1]; ...
        
    % Evaluate Jacobian:
    Hk = H(xhat);

    % Solve Least Squares
    delXi = z - h(xhat);
    delXV = (Hk'*Hk)\Hk'*delXi;

    % Adjust x(4) Estimate w.r.t. Time:
    delXV = [delXV(1:3); delXV(4)/c];
    xhat = [xhat(1:3); xhat(4)/c];
    
    % Save:
    xhat_est = [xhat_est xhat];

    % Update:
    xhat0 = xhat + delXV;
    delXMag = norm(delXV(1:3));
%     delXMag = norm(delXV);
    count = count + 1;
end

% Print Results:
fprintf('Navigation Solution: \n');

% Final Estimated State:
x_est = xhat_est(:, end);
fprintf('Estimated States (xr, yr, zr, cDeltatr) = [%4.4f m, %4.4f m, %4.4f m, %4.4e sec] \n', x_est);

% Longitude and Latitude:
p = ecef2lla(x_est(1:3)');
fprintf('Latitude = %4.4f, Longitude = %4.4f, and Altitude = %4.4f \n\n', p(1), p(2), p(3))

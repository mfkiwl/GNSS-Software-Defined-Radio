function [Ad, Bd, Cd, Dd, BL_act] = configureLoopFilter(BL_target, Tsub, loopOrder)
% configureLoopFilter : Configure a discrete-time loop filter for a feedback
% tracking loop.
%
% INPUTS
%
% BL_target ----- Target loop noise bandwidth of the closed-loop system, in Hz.
%
% Tsub ---------- Subaccumulation interval, in seconds. This is also the loop
% update (discretization) interval.
%
% loopOrder ----- The order of the closed-loop system. Possible choices are
% 1, 2, or 3.
%
% OUTPUTS
%
% Ad,Bd,Cd,Dd --- Discrete-time state-space model of the loop filter.
%
% BL_act -------- The actual loop noise bandwidth (in Hz) of the closed-loop
% tracking loop as determined by taking into account the
% discretized loop filter, the implicit integration of the
% carrier phase estimate, and the length of the accumulation
% interval.
%
%--------------------------------------------------------------------------

if loopOrder == 1                % First-Order Case
% Define Parameter K:
K = 4*BL_target;

% Type 1 Transfer Function:
Fs = K*tf(1, 1);

elseif loopOrder == 2            % Second-Order Case
% Define Parameters K & a:
K = 8*BL_target/3;
a = K/2;

% Type 2 Transfer Function:
Fs = K*tf([1 a], [1 0]);

elseif loopOrder == 3            % Third-Order Case
% Define Parameters K, a, & b:
a = 1.2*BL_target;
b = a.^2/2;
K = 2*a;

% Type 3 Transfer Function:
Fs = K*tf([1 a b], [1 0 0]);

else                             % Wrong Loop-Order Choice
    error('Choose the Loop Order to be Either 1, 2, or 3!');
end

% Convert Loop Filter to a DT State-Space Model:
Fz = c2d(Fs, Tsub, 'zoh');       % DT Transfer Function
[Ad, Bd, Cd, Dd] = ssdata(Fz);

% Linearized DT Costas Open-Loop System:
NCO = Tsub*tf(1, [1 -1], Tsub);    % Number-Controlled Oscillator (NCO)
PD = 1/2*tf([1 1], [1 0], Tsub); % DT Averager
sysOpenLoop = PD*Fz*NCO;    

% Linearized DT Costas Closed-Loop System:
H_z = sysOpenLoop/(1 + sysOpenLoop);

% Calculate Actual Loop Noise Bandwidth:
walias = pi/Tsub;
wvec = (0:10000)'*(walias/10000);
[magvec, ~] = bode(H_z, wvec);
magvec = magvec(:);
BL_act = sum(magvec.^2)*(wvec(2) - wvec(1))/(2*pi); % Approx. Bn 

end



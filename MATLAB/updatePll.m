function [xkp1, vk, ek_pll] = updatePll(s)
% updatePll : Perform a single update step of a phase tracking loop with an
% arctangent phase detector.
%
% INPUTS
%
% s ------------- A structure with the following fields:
%
% Ipk -------- The in-phase prompt accumulation over the interval from
% tkm1 to tk.
%
% Qpk -------- The quadrature prompt accumulation over the interval from
% tkm1 to tk.
%
% xk --------- The phase tracking loop filter’s state at time tk. The
% dimension of xk is N-1, where N is the order of the loop’s
% closed-loop transfer function.
%
% Ad,Bd,Cd,Dd -- The loop filter’s state-space model.
%
% OUTPUTS
%
% xkp1 -------- The loop filter’s state at time tkp1. The dimension of xkp1
% is N-1, where N is the order of the loop’s closed-loop
% transfer function.
%
% vk ---------- The Doppler frequency shift that will be used to drive the
% receiver’s carrier-tracking numerically controlled
% oscillator during the time interval from tk to tkp1, in
% rad/sec.
%
%+------------------------------------------------------------------------------+

% Phase Detector (Discriminator):
ek_pll = atan(s.Qpk./s.Ipk);           % Arc-Tangent (AT)

% Next Time Step Calculations (@ tkp1):
xkp1 = s.Ad*s.xk + s.Bd*ek_pll;        % Loop Filter States  
vk = s.Cd*s.xk + s.Dd*ek_pll;          % Doppler Frequency Shift 

end





function [Spk, Sek, Slk] = CorrelatorBlk(Y, Gpm, fsampIQ, ts_hat, th_hat, teml)
% Inputs 
%       Y: Complex Signal
%     Gpm: Oversampled PRN # 
% fsampIQ: Signal Sampling Frequency
%  ts_hat: Code Start Time Estimate 
%  th_hat: Beat Carrier Phase Estimate
%    teml: Early-Minus-Late Time 
%
% Outputs
%     Spk: Prompt Signal
%     Sek: Early Signal
%     Slk: Late Signal
%-------------------------------------------------------------------------%

% Circular Shift PRN:
Tprom   = circshift(Gpm, round(ts_hat*fsampIQ));
TemlAdv = circshift(Gpm, floor((ts_hat - teml/2)*fsampIQ));
TemlDel = circshift(Gpm, floor((ts_hat + teml/2)*fsampIQ));

% Wipe Real Signal with Doppler and Code Start Estimates:
Sk = Y.*exp(-1i*th_hat);

% Summate Subaccums:
Spk = sum(Sk.*Tprom);    % Prompt
Sek = sum(Sk.*TemlAdv);  % Advance
Slk = sum(Sk.*TemlDel);  % Delay

end
function [codeOS] = oversampleSpreadingCode(code,delChip,Ns,Np)
%[codeOS] = oversampleSpreadingCode(code,delChip,Ns,Np)
%
% Oversample the +/-1-valued input pseudorandom spreading code.
%
% INPUTS
%
% code ------------ N-by-1 vector containing the +/-1-valued input
%                   pseudorandom spreading code.  The spreading code is
%                   assumed to be periodic with period Np.
%
% delChip --------- Sampling interval, measured in spreading code chip
%                   intervals.  To ensure good auto- and cross-correlation
%                   properties, the sampling rate should not be an integer
%                   multiple of the chipping rate.  Thus, n*delChip should not
%                   equal 1 for any integer n.
% 
% Ns -------------- Number of samples required for the output oversampled
%                   code.
% 
% Np -------------- Assumed number of chips in one spreading code period.
%
%
% OUTPUTS
%
% codeOS ---------- Ns-by-1 vector containing the +/-1-valued oversampled
%                   spreading code.

if((sum(code == 1) + sum(code == -1))~= length(code))
  error('Spreading code is assumed to be +/-1-valued');
end

if(abs(mod(1/delChip,1)) < 1e-6)
  error(['To ensure good auto- and cross-correlation '...
         'properties, the sampling rate should not be an integer '...
         'multiple of the chipping rate.  Thus, n*delChip should not '...
         'equal 1 for any integer n.']);
end

codeOS = zeros(Ns,1);
for ii=1:Ns
  k = floor((ii-1)*delChip);
  k = mod(k, Np);
  codeOS(ii) = code(k+1);
end
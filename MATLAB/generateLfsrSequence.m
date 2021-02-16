function [lfsrSeq] = generateLfsrSequence(n,ciVec,a0Vec)
%
% Generate a 1/0-valued linear feedback shift register (LFSR) sequence.
%
% INPUTS
%
% n ------ Number of stages in the linear feedback shift register.
%
% ciVec -- Nc-by-1 vector whose elements give the indices of the Nc nonzero
% connection elements. For example, if the characteristic polynomial
% of an LFSR is f(D) = 1 + D^2 + D^3, then ciVec = [2,3]’ or [3,2]’.
%
% a0Vec -- n-by-1 1/0-valued initial state of the LFSR, where a0Vec = [a(-1),
% a(-2), ..., a(-n)]’. In defining the initial LFSR state, a
% Fibonacci LFSR implementation is assumed.
%
% OUTPUTS
%
% lfsrSeq -- m-by-1 vector whose elements are the 1/0-valued LFSR sequence
% corresponding to n, ciVec, and a0Vec, where m = 2^n - 1. If the
% sequence is a maximal-length sequence, then there is no
% repetition in the m sequence elements.
%
%+------------------------------------------------------------------------------+
 
m = 2^n - 1;             % LSFR Sequence Period
Nc = length(ciVec);      % Characterisitc Polynomial Length
c(ciVec(:)) = 1;         % Connection Vector 
a = a0Vec;               % LSFR State (Assuming Fibonacci Implementation)
lfsrSeq = zeros(m, 1);   % LFSR Sequence Preallocation
b = zeros(Nc - 1, 0);    % Logic Gate Result Preallocation

for k = 1:m-1
    % Logic Gate Operation(s)
    b(1) = xor(a(ciVec(1)), a(ciVec(2)));

    if Nc > 2
        for i = 1:Nc - 2
            
            b(i + 1) = xor(a(ciVec(i + 2)), b(i));

        end
    end
    
    % Shifting LSFR Sequence
    j = 1:n-1;
    a(n + 1 - j) = a(n - j);
    a(1) = b(Nc - 1);
    
    % Save LSFR Sequence Output
    lfsrSeq(k + 1) = a(end);
end

lfsrSeq(1) = a0Vec(end); % Initial LFSR Sequence Output

end

function r = delta2ratio(deltaD)
%DELTA2RATIO Convert delta2H (per mil VSMOW) to normalized isotope ratio.
%
%   r = R / R_VSMOW
%
%   Therefore:
%       r = 1 + deltaD/1000


r = 1 + deltaD ./ 1000;

end

function D = effective5_latent2_bounds(B)
%EFFECTIVE5_LATENT2_BOUNDS
% Conservative rectangular L1/L2 bounds derived from the 32 corners of
% the five-state admissible hyperrectangle.
%
% The inverse also rejects candidate L values whose representative
% E vector falls outside the original E-state bounds.

lo=[B.parameterDefinitions.low].';
hi=[B.parameterDefinitions.high].';

C=dec2bin(0:31)-'0';
Lc=NaN(32,2);

for i=1:32
    E=lo+(hi-lo).*C(i,:).';
    Lc(i,:)=effective5_to_latent2(E,B).';
end

D=struct( ...
    'id',{'L1','L2'}, ...
    'name',{'overallBiochemicalLevel','homologSlopeContrast'}, ...
    'low',{min(Lc(:,1)),min(Lc(:,2))}, ...
    'high',{max(Lc(:,1)),max(Lc(:,2))});

end

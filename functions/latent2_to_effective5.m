
function E = latent2_to_effective5(L,B)
%LATENT2_TO_EFFECTIVE5
% Representative five-state vector for L1/L2 with all three near-null
% coordinates fixed to zero at the reference.

L=L(:);

if numel(L)~=2
    error('L must contain [L1 L2].');
end

z=B.V2*L;
E=B.referenceEffective+B.parameterRange.*z;

end

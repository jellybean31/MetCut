
function [L,z,zResidual] = effective5_to_latent2(E,B)
%EFFECTIVE5_TO_LATENT2 Project five effective states onto L1/L2.

E=E(:);

if numel(E)~=5
    error('E must contain E01-E05.');
end

z=(E-B.referenceEffective)./B.parameterRange;
L=B.V2.'*z;

zProjected=B.V2*L;
zResidual=z-zProjected;

end

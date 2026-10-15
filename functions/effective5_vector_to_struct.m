
function S = effective5_vector_to_struct(E)
%EFFECTIVE5_VECTOR_TO_STRUCT Convert [E01..E05] vector to bridge struct.

E=E(:);
if numel(E)~=5
    error('Expected five effective biochemical states E01-E05.');
end

S=struct;
S.acetylStarter_permil      = E(1);
S.retainedPrecursor_permil = E(2);
S.plastidNADPH_permil      = E(3);
S.plastidNADH_permil       = E(4);
S.erNADPH_permil           = E(5);
end

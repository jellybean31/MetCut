function F = effective7_structural4_inverse_forward( ...
    M,intracellularNode,timeDay,P,Z)
%EFFECTIVE7_STRUCTURAL4_INVERSE_FORWARD
% Forward map for architecture-v2 inverse calibration:
%
%   Z -> Q -> Erep -> wax deltaD
%
% Inputs
% ------
% M                 validated model struct
% intracellularNode intracellular state at one time node
% timeDay           scalar model time
% P                 hybrid structural4 prior from
%                   effective7_structural4_hybrid_inverse_prior
% Z                 4xN standardized inverse coordinates
%
% Output
% ------
% F struct with Q, Erep, and homolog deltaD matrix.

if size(Z,1)~=4
    error('Z must be 4 x N.');
end

[Q,~]=effective7_structural4_hybrid_transform(P,Z,'z_to_q');

n=size(Z,2);
E=NaN(7,n);
wax=NaN(4,n);

for i=1:n
    E(:,i)=P.referenceEffective7 + P.rawStructuralBasis*Q(:,i);
    wax(:,i)=effective7_node_wax( ...
        M,intracellularNode,timeDay,E(:,i));
end

F=struct;
F.version='2026-08-29_effective7_structural4_inverse_forward_v1';
F.Z=Z;
F.Q=Q;
F.Effective7Representative=E;
F.waxDeltaD_permil=wax;

end

function [Q,Z] = effective7_structural4_hybrid_transform(P,input,mode)
%EFFECTIVE7_STRUCTURAL4_HYBRID_TRANSFORM
% Transform between standardized inverse coordinates Z and raw
% structural coordinates Q using a hybrid prior.
%
% mode='z_to_q':
%   Q = mean + factor * Z
%
% mode='q_to_z':
%   Z = inverseFactor * (Q - mean)
%
% Columns are samples:
%   Z is 4 x N
%   Q is 4 x N

mode=lower(strtrim(string(mode)));

switch mode
    case "z_to_q"
        Z=input;
        if size(Z,1)~=4
            error('Z must be 4 x N.');
        end
        Q=P.mean + P.factor*Z;

    case "q_to_z"
        Q=input;
        if size(Q,1)~=4
            error('Q must be 4 x N.');
        end
        Z=P.inverseFactor*(Q-P.mean);

    otherwise
        error('mode must be ''z_to_q'' or ''q_to_z''.');
end

end

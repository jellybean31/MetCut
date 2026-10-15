function y = effective7_node_wax( ...
    Mreference,intracellularNode,timeDay,E7)
%EFFECTIVE7_NODE_WAX Exact architecture-v2 C27/C29/C31/C33 wax vector.
%
% E7 order:
%   1 E01 acetylStarter_permil
%   2 E02 retainedPrecursor_permil
%   3 E03 plastidNADPH_permil
%   4 E04 plastidNADH_permil
%   5 E05 erNADPH_permil
%   6 E06 stage6Delta_permil
%   7 E07 stage7Delta_permil
%
% The evaluator reuses the validated effective5 biochemical
% bridge for E01-E05, then applies the validated exact-D ER stage-state
% adapter. 
%
% With E06=E07=0, this must reproduce effective5_node_wax.

E7=E7(:);
if numel(E7)~=7
    error('Expected seven effective biochemical coordinates E01-E07.');
end
if any(~isfinite(E7))
    error('Effective biochemical coordinates must be finite.');
end

cfg=Mreference.config.biochemistry;
effectiveStates=effective5_vector_to_struct(E7(1:5));

if isfield(Mreference.biochemistry,'reactionAlphaValues')
    reactionAlpha=Mreference.biochemistry.reactionAlphaValues;
else
    reactionAlpha=[];
end

Bio=wax_biochemistry_from_intracellular_effective5( ...
    timeDay,intracellularNode,cfg,effectiveStates,reactionAlpha);

stageTraits=wax_er_stage_state_traits();
stageTraits.stage6Delta_permil=E7(6);
stageTraits.stage7Delta_permil=E7(7);

W2=apply_wax_er_stage_states_exactD( ...
    Bio.wax.deltaDMatrix,Bio.alkaneCarbonNumber,stageTraits);

y=W2.deltaDMatrix(:,1);

end

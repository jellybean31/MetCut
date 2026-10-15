
function M = attach_wax_architecture_v2_to_model(M,stageTraits)
%ATTACH_WAX_ARCHITECTURE_V2_TO_MODEL
% Apply architecture-v2 hook to an existing model output.
%
% Preserves:
%   M.biochemistry.wax
%
% Adds:
%   M.biochemistry.waxStageExpanded
%   M.biochemistry.architecture

if nargin<2 || isempty(stageTraits)
    stageTraits=wax_er_stage_state_traits();
end

if ~isfield(M,'biochemistry')
    error('M.biochemistry is required.');
end

M.biochemistry=apply_wax_architecture_v2_to_biochemistry( ...
    M.biochemistry,stageTraits);

end

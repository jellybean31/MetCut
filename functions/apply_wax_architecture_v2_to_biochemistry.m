
function B = apply_wax_architecture_v2_to_biochemistry(B,stageTraits)
%APPLY_WAX_ARCHITECTURE_V2_TO_BIOCHEMISTRY
% Production-safe architecture-v2 hook.
%
% Preserves architecture-v1 fields, especially:
%   B.wax
%
% Adds:
%   B.waxStageExpanded
%   B.architecture


if nargin<2 || isempty(stageTraits)
    stageTraits=wax_er_stage_state_traits();
end

if ~isstruct(B)
    error('B must be a biochemical output structure.');
end

if ~isfield(B,'wax') || ~isfield(B.wax,'deltaDMatrix')
    error('B.wax.deltaDMatrix is required.');
end

if isfield(B,'alkaneCarbonNumber')
    C=B.alkaneCarbonNumber;
elseif isfield(B.wax,'alkaneCarbonNumber')
    C=B.wax.alkaneCarbonNumber;
else
    error(['Could not resolve alkane carbon numbers from either ' ...
        'B.alkaneCarbonNumber or B.wax.alkaneCarbonNumber.']);
end

W=apply_wax_er_stage_states_exactD( ...
    B.wax.deltaDMatrix,C,stageTraits);

B.waxStageExpanded=W;

B.architecture=struct;
B.architecture.name='wax-biochemistry-v2';
B.architecture.baseArchitecture='wax-biochemistry-v1';
B.architecture.stageStateImplementation='effectiveERStageAdapter';
B.architecture.stageTraits=stageTraits;
B.architecture.zeroStateNestsV1=true;

end

function [M,B] = rerun_full_model_isotope_forcing_with_wax_timing( ...
    Mreference,deltaDPrecipInterval_permil,stageTraits,accumTraits)
%RERUN_FULL_MODEL_ISOTOPE_FORCING_WITH_WAX_TIMING
%
%
% INPUTS
% ------
% Mreference
%   integrated reference model
%
% deltaDPrecipInterval_permil
%   nTime-1 precipitation isotope forcing
%
% stageTraits
%   architecture-v2 ER-stage traits
%   default: zero stage6/stage7 states
%
% accumTraits
%   wax accumulation timing traits
%   default: wax_accumulation_default_traits()
%
% OUTPUTS
% -------
% M
%   replayed model with:
%       M.biochemistry.waxStageExpanded
%       M.biochemistry.waxAccumulation
%
% B
%   original fractionation budget from the replay.
%


if nargin<3 || isempty(stageTraits)
    stageTraits=wax_er_stage_state_traits();
    stageTraits.enabled=true;
    stageTraits.stage6Delta_permil=0;
    stageTraits.stage7Delta_permil=0;
end

if nargin<4 || isempty(accumTraits)
    accumTraits=wax_accumulation_default_traits();
end

[M,B]=rerun_full_model_isotope_forcing( ...
    Mreference,deltaDPrecipInterval_permil);

M=attach_wax_architecture_v2_to_model(M,stageTraits);

M=attach_wax_accumulation_to_model( ...
    M,accumTraits,'waxStageExpanded');

M.temporalTransferExtension=struct;
M.temporalTransferExtension.version= ...
    '2026-08-30_H6_architecture_v2_wax_timing_replay_v1';
M.temporalTransferExtension.architectureSource='waxStageExpanded';
M.temporalTransferExtension.stageTraits=stageTraits;
M.temporalTransferExtension.accumulationTraits=accumTraits;
M.temporalTransferExtension.hydraulicsReused=true;

end

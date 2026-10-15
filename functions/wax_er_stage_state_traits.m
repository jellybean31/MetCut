
function T = wax_er_stage_state_traits()
%WAX_ER_STAGE_STATE_TRAITS Default architecture-v2 ER stage states.
T=struct;
T.enabled=true;
T.stage6Delta_permil=0;
T.stage7Delta_permil=0;
T.hydrogenAtomsPerStage=4;
T.stage6MinimumCarbonNumber=29;
T.stage7MinimumCarbonNumber=31;
T.stateMode='constantAlphaRelativeToBaselineProduct';
T.Rstd_VSMOW=155.76e-6;
T.description=[ ...
    'Effective ER stage states applied in exact D-atom space. ' ...
    'Zero-state identity is mandatory.'];
end

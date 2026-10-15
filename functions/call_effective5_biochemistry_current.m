
function Bio = call_effective5_biochemistry_current(Mreference,Mup,E)
%CALL_EFFECTIVE5_BIOCHEMISTRY_CURRENT
% Exact adapter for current project signature:
%
% Bio = wax_biochemistry_from_intracellular_effective5( ...
%     timeDays,intracellular,cfg,effectiveStates,reactionAlpha)

if ~isfield(Mup,'intracellular')
    error('Candidate upstream model does not contain intracellular.');
end

intracellular=Mup.intracellular;

if isfield(intracellular,'timeDays')
    timeDays=intracellular.timeDays;
elseif isfield(Mreference,'biochemistry') && ...
        isfield(Mreference.biochemistry,'timeDays')
    timeDays=Mreference.biochemistry.timeDays;
else
    error('Could not resolve timeDays.');
end

if isfield(Mreference,'config') && ...
        isfield(Mreference.config,'biochemistry')
    cfg=Mreference.config.biochemistry;
else
    cfg=struct;
end

effectiveStates=effective5_vector_to_struct(E);

if isfield(Mreference,'biochemistry') && ...
        isfield(Mreference.biochemistry,'reactionAlphaValues') && ...
        ~isempty(Mreference.biochemistry.reactionAlphaValues)

    reactionAlpha=Mreference.biochemistry.reactionAlphaValues;
else
    % The bridge itself has a documented midpoint fallback when empty.
    reactionAlpha=[];
end

Bio=wax_biochemistry_from_intracellular_effective5( ...
    timeDays,intracellular,cfg,effectiveStates,reactionAlpha);

end

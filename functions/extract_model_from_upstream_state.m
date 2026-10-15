
function M = extract_model_from_upstream_state(U)
%EXTRACT_MODEL_FROM_UPSTREAM_STATE Resolve model carrying intracellular state.

if isstruct(U) && isfield(U,'intracellular')
    M=U;
    return
end

paths={ ...
    {'model'}, ...
    {'M'}, ...
    {'state'}, ...
    {'upstreamModel'}};

for k=1:numel(paths)
    p=paths{k}{1};
    if isstruct(U) && isfield(U,p) && isstruct(U.(p)) && ...
            isfield(U.(p),'intracellular')
        M=U.(p);
        return
    end
end

error(['Unable to locate intracellular state inside Ustate. ' ...
    'Inspect fieldnames(Ustate) from the preflight candidate test.']);
end

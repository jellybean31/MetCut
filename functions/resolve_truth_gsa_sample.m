
function sample = resolve_truth_gsa_sample(Rpilot)
%RESOLVE_TRUTH_GSA_SAMPLE Reconstruct retained synthetic-truth GSA sample.

sample=struct;
sample.I01=2.03690139;
sample.I02=5.38306961;

M=Rpilot.truth.model;

if isfield(M,'biochemistry') && ...
        isfield(M.biochemistry,'sourceParameterValues')

    v=M.biochemistry.sourceParameterValues(:).';

    for j=1:min(7,numel(v))
        sample.(sprintf('B%02d',j))=v(j);
    end
end

fprintf(['resolve_truth_gsa_sample: using parameter ' ...
    'truth sample reconstructed from pilot.\n']);
end

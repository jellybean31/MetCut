
function T = synthetic_recovery_make_truth_sample(seed)
%SYNTHETIC_RECOVERY_MAKE_TRUTH_SAMPLE
% Deterministic interior truth point for identifiability tests.

if nargin<1 || isempty(seed), seed=20260828; end

D=synthetic_recovery_parameter_definitions();


rng(seed,'twister');

vals=NaN(1,numel(D.integrated));
names=cell(1,numel(D.integrated));

u=[0.61 0.43 0.35 0.72 0.58 0.39 0.66 0.52 0.74];

for j=1:numel(D.integrated)
    q=D.integrated(j);

    switch lower(q.distribution)
        case 'uniform'
            vals(j)=q.low+(q.high-q.low)*u(j);
        case 'loguniform'
            vals(j)=exp(log(q.low)+(log(q.high)-log(q.low))*u(j));
        otherwise
            error('Unsupported distribution.');
    end

    names{j}=char(matlab.lang.makeValidName( ...
        string(q.id)+"_"+string(q.name)));
end

T=array2table(vals,'VariableNames',names);
T=addvars(T,1,'Before',1,'NewVariableNames','RunID');

end

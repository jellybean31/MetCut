function y = effective5_predict_final_homologs(Truth,x)
%EFFECTIVE5_PREDICT_FINAL_HOMOLOGS Predict final C27/C29/C31/C33 deltaD.
x=x(:).';
if numel(x)~=5, error('x must contain five effective biochemical states.'); end
I=Truth.model.intracellular;
cfg=Truth.model.config.biochemistry;
E=struct('acetylStarter_permil',x(1), ...
    'retainedPrecursor_permil',x(2), ...
    'plastidNADPH_permil',x(3), ...
    'plastidNADH_permil',x(4), ...
    'erNADPH_permil',x(5));
P=wax_integrated_biochemical_prior_definitions();
a=0.5.*(P.core.reaction.lowerAlpha(:).'+P.core.reaction.upperAlpha(:).');
a(7)=1;
Bio=wax_biochemistry_from_intracellular_effective5(Truth.model.timeDays,I,cfg,E,a);
chains=Bio.alkaneCarbonNumber(:);
homologs=Truth.config.biochemistry.homologs(:);
y=NaN(numel(homologs),1);
for k=1:numel(homologs)
    q=find(chains==homologs(k),1);
    if isempty(q), error('Homolog C%d not found.',homologs(k)); end
    y(k)=Bio.wax.deltaDMatrix(q,end);
end
end

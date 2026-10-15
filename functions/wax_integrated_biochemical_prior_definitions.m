function P = wax_integrated_biochemical_prior_definitions()
%WAX_INTEGRATED_BIOCHEMICAL_PRIOR_DEFINITIONS
% Development pre-priors for integrated leaf-wax biochemical uncertainty.




P = struct;
P.status = 'integrated biochemical prior-development bounds; not final physiological priors';
P.sourceReferenceWater_permil = -100;
P.sourceOffsetDefinition = ['Delta_source,ref = deltaD_preExchange_source(ref) - ', ...
    'deltaD_water(ref); propagated away from the reference as constant alpha'];
P.exchangeDefinition = 'f_exchange is the exact-D atom-fraction mixing fraction with local water';
P.sourcePropagationDefinition = ['alpha_source = R_source(ref)/R_water(ref); ', ...
    'R_source(local) = alpha_source * R_water(local)'];
P.reactionPriorDefinitions = wax_reaction_prior_definitions();

P.sourceParameterNames = { ...
    'delta_acetylStarter_vs_water', ...
    'delta_pNADPH_preExchange_vs_water', ...
    'delta_pNADH_preExchange_vs_water', ...
    'delta_erNADPH_preExchange_vs_water', ...
    'f_exchange_pNADPH', ...
    'f_exchange_erNADPH', ...
    'f_postMalonate_plastid', ...
    'f_postMalonate_er'};

P.sourceParameterType = { ...
    'isotopeOffset','isotopeOffset','isotopeOffset','isotopeOffset', ...
    'exchangeFraction','exchangeFraction','exchangeFraction','exchangeFraction'};

P.core.sourceLower = [-150 -200 -150 -175   0 0   0 0];
P.core.sourceUpper = [ -50  -50  -25  -25 0.5 0.5 0.5 0.5];

P.broad.sourceLower = [-200 -300 -250 -250   0 0   0 0];
P.broad.sourceUpper = [   0   50   50   50   1 1   1 1];

P.outer.sourceLower = [-300 -500 -400 -400   0 0   0 0];
P.outer.sourceUpper = [ 100  200  200  200   1 1   1 1];

P.core.reaction = P.reactionPriorDefinitions.core;
P.broad.reaction = P.reactionPriorDefinitions.broad;
P.outer.reaction = P.reactionPriorDefinitions.outer;

P.core.hierarchicalResidualFraction = 0.15;
P.broad.hierarchicalResidualFraction = 0.15;
P.outer.hierarchicalResidualFraction = 0.15;

% Three alternative reaction prior structures.
P.reactionStructures = {'independent','grouped','hierarchical'};
P.reactionStructureDescription.independent = [ ...
    'Six reaction alphas sampled independently within the tier bounds.'];
P.reactionStructureDescription.grouped = [ ...
    'One shared plastid-hydride alpha, one shared ER-hydride alpha, ', ...
    'and one shared proton alpha; CER3 fixed at one.'];
P.reactionStructureDescription.hierarchical = [ ...
    'Latent plastid-hydride, ER-hydride, and proton group means with ', ...
    'smaller reaction-specific bounded deviations.'];

% Currently supported upstream uncertainty choices.
P.plastidNADHWaterExchange.fixedFraction = 0;
P.plastidNADHWaterExchange.reason = [ ...
    'Retained at zero in this development analysis because the present ', ...
    'literature review specifically supports NADPH-water exchange more ', ...
    'strongly than an independently parameterized plastid NADH exchange.'];
P.cer3Source = 'erNADPH';
P.cer3AlphaFixedAtOne = true;
P.pathwayFractionsSampled = false;
P.pathwayFractionsReason = [ ...
    'Effective pre-exchange pool offsets are sampled instead of separately ', ...
    'sampling non-identifiable pathway fractions and endmember contrasts.'];

% Build source-bound tables and validate nesting.
tiers = {'core','broad','outer'};
for i=1:numel(tiers)
    nm = tiers{i};
    lo = P.(nm).sourceLower;
    hi = P.(nm).sourceUpper;
    if numel(lo)~=numel(P.sourceParameterNames) || numel(hi)~=numel(P.sourceParameterNames)
        error('Source-prior bounds do not match source parameter names.');
    end
    if any(~isfinite(lo)) || any(~isfinite(hi)) || any(lo>hi)
        error('Invalid source-prior bounds in %s tier.',nm);
    end
    fIdx = strcmp(P.sourceParameterType,'exchangeFraction');
    if any(lo(fIdx)<0) || any(hi(fIdx)>1)
        error('Exchange-fraction bounds must lie in [0,1].');
    end
    P.(nm).sourceTable = table(P.sourceParameterNames(:),P.sourceParameterType(:), ...
        lo(:),hi(:),'VariableNames',{'Parameter','Type','Lower','Upper'});
end

if any(P.core.sourceLower < P.broad.sourceLower) || ...
        any(P.core.sourceUpper > P.broad.sourceUpper)
    error('Core source bounds must lie inside broad bounds.');
end
if any(P.broad.sourceLower < P.outer.sourceLower) || ...
        any(P.broad.sourceUpper > P.outer.sourceUpper)
    error('Broad source bounds must lie inside outer bounds.');
end
end

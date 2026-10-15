function P = wax_reaction_prior_definitions()
%WAX_REACTION_PRIOR_DEFINITIONS
% Evidence-tiered PRE-PRIORS for reaction-specific H incorporation alpha.
%
% alpha_inc = R_incorporated / R_source
% epsilon_inc = 1000*(alpha_inc - 1) per mil
%
% Tiers
%   reference : null regression value alpha = 1
%   core      : conservative working effective range
%   broad     : wider working range for prior robustness
%   outer     : mechanistic stress-test envelope; NOT a physiological prior
%
% Hydride outer lower bound 1/2.7 is retained as the diagnostic limit
% corresponding to complete expression of the FabG-like intrinsic KIE.
%
% Proton terms are effective water->product parameters and may include
% proton-donor equilibration and kinetic expression.

P = struct;
P.definition = 'alpha_inc = R_incorporated / R_source';
P.parameterization = 'epsilon_inc = 1000*(alpha_inc - 1)';
P.status = 'pre-prior-development ranges; not final physiological priors';

names = { ...
    'alpha_pKAR', ...
    'alpha_pENR_hydride', ...
    'alpha_pENR_proton', ...
    'alpha_erKCR', ...
    'alpha_erECR_hydride', ...
    'alpha_erECR_proton', ...
    'alpha_CER3'};

types = { ...
    'hydride','hydride','proton', ...
    'hydride','hydride','proton','cer3'};

coreLo = [0.80 0.80 0.85 0.80 0.80 0.85 1.00];
coreHi = [1.00 1.00 1.05 1.00 1.00 1.05 1.00];

broadLo = [0.70 0.70 0.80 0.70 0.70 0.80 1.00];
broadHi = [1.00 1.00 1.10 1.00 1.00 1.10 1.00];

outerHydrideLo = 1/2.7;
outerLo = [outerHydrideLo outerHydrideLo 0.50 ...
           outerHydrideLo outerHydrideLo 0.50 1.00];
outerHi = [1.00 1.00 1.10 1.00 1.00 1.10 1.00];

reference = ones(1,numel(names));

P.parameterNames = names;
P.parameterType = types;
P.referenceAlpha = reference;
P.core = makeTier(coreLo,coreHi,names,types,'core');
P.broad = makeTier(broadLo,broadHi,names,types,'broad');
P.outer = makeTier(outerLo,outerHi,names,types,'outer mechanistic');

P.dynamicModulation.epsilonHalfRange_permil = [50 50 50 50 50 50 0];
P.dynamicModulation.description = [ ...
    'Future within-plant/time-varying modulation envelope around a sampled ', ...
    'baseline reaction epsilon. Development placeholder, not final prior.'];

P.cer3FixedAtOne = true;
P.independenceAssumption = [ ...
    'Current Monte Carlo pre-priors sample reaction baselines independently. ', ...
    'Later hierarchical/correlated priors may be warranted.'];
P.evidenceInterpretation = { ...
    'pKAR: strongest hydride-KIE analogue; intrinsic KIE is outer bound only'; ...
    'pENR hydride: separate reaction retained; no direct plant product alpha'; ...
    'pENR proton: effective water-to-product term; not a pure enzyme KIE'; ...
    'ER KCR: hydride-transfer analogue, but no direct plant product alpha'; ...
    'ER ECR hydride: distinct enzyme architecture; independent broad prior'; ...
    'ER ECR proton: effective water-to-product term'; ...
    'CER3: fixed at unity until stronger isotope evidence exists'};
end

function T = makeTier(lo,hi,names,types,label)
if any(lo<=0) || any(hi<=0) || any(lo>hi)
    error('Invalid alpha bounds in %s tier.',label);
end
T = struct;
T.label = label;
T.lowerAlpha = lo;
T.upperAlpha = hi;
T.lowerEpsilon_permil = 1000.*(lo-1);
T.upperEpsilon_permil = 1000.*(hi-1);
T.parameterTable = table(names(:),types(:),lo(:),hi(:), ...
    T.lowerEpsilon_permil(:),T.upperEpsilon_permil(:), ...
    'VariableNames',{'ReactionParameter','Type','LowerAlpha','UpperAlpha', ...
    'LowerEpsilon_permil','UpperEpsilon_permil'});
end

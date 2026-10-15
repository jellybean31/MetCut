function T = wax_reaction_incorporation_traits_null()
%WAX_REACTION_INCORPORATION_TRAITS_NULL
% Null reaction-incorporation isotope-effect configuration.
%
% alpha is defined as:
%
%   alpha_inc = R_H,incorporated / R_H,source
%
% All values are 1 here, so this configuration must exactly reproduce the
% post-malonate source-only model after reaction classes are resolved.

T = struct;
T.plastid.KAR.alphaHydride = 1;
T.plastid.ENR.alphaHydride = 1;
T.plastid.ENR.alphaProton = 1;
T.er.KCR.alphaHydride = 1;
T.er.ECR.alphaHydride = 1;
T.er.ECR.alphaProton = 1;
T.cer3.alphaAldehydeHydride = 1;

T.definition = 'alpha_inc = R_incorporated / R_source';
T.isNullRegressionConfiguration = true;
end

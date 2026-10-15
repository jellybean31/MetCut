function P = rooting_arch_mapping_payload()
%ROOTING_ARCH_MAPPING_PAYLOAD_CANDIDATE
%
% Candidate deterministic constructor for the numerical payload
% consumed by evaluate_rooting_arch_mapping.


P=struct;

P.predictorNames=[ ...
"intercept"; ...
"z_logPET"; ...
"z_logMAPoverPET"; ...
"broadleavedTrees"; ...
"needleleavedTrees"; ...
"shrubs"; ...
"semiShrubs"; ...
"grasses"; ...
"forbs"];

P.coefficients=[ ...
3.4182016958106329      1.5426497415679388
-0.48744778771303471   -0.2507382877658707
0.5234049723685541      0.16827777938041252
-0.75400782621402462   -0.22703225936574409
-0.41427373127496842    0.011051846840582913
0.24592809435179327     0.16401490850119863
0.49053616591236671     0.30499917683684347
0.20319018039710179     0.032797480597302586
0.48227057459808387     0.21472015597702063];

P.PET_logMean=6.8633590023910553;
P.PET_logSD=0.43972269727668833;
P.aridity_logMean=-0.26708054840958101;
P.aridity_logSD=0.86795300314577317;
P.responseParameterization="ALR_f1_over_f3__f2_over_f3";
P.outputFamily="logistic_normal_simplex";

end

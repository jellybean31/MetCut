function L = leaf_parallel_hydraulic_network( ...
    psiS, ...
    E, ...
    dWSdt, ...
    KXA, ...
    KAS, ...
    KSE, ...
    KAE)
%LEAF_PARALLEL_HYDRAULIC_NETWORK
%
% Algebraic two-domain leaf hydraulic network.
%
%
% Topology:
%
%                         KAE
%               +-------------------->
%               |
%               |
%   xylem ----> A
%      KXA      |
%               |
%               +----> S ------------>
%                    KAS       KSE
%
%                               evaporation site
%
%
% Hydraulic states:
%
%   psiS
%       finite symplastic/cellular hydraulic state
%
%   psiA
%       massless apoplastic hydraulic node
%
%   psiE
%       massless liquid-side evaporation-site hydraulic node
%
%
% Water fluxes:
%
%   JXA
%       xylem -> apoplast
%
%   JAE
%       apoplast -> evaporation site
%
%   JAS
%       apoplast -> symplast
%
%   JSE
%       symplast -> evaporation site
%
%
% Symplastic storage:
%
%   dWS/dt =
%       JAS - JSE
%
%
% Evaporation:
%
%   E =
%       JAE + JSE
%
%
% Apoplast:
%
%   JXA =
%       JAE + JAS



%% ========================================================================
% Validate inputs
% ========================================================================

inputValues = [ ...
    psiS, ...
    E, ...
    dWSdt, ...
    KXA, ...
    KAS, ...
    KSE, ...
    KAE];


if any(~isfinite(inputValues))

    error( ...
        'All leaf hydraulic-network inputs must be finite.');

end


if KXA <= 0

    error( ...
        'KXA must be positive.');

end


if KAS <= 0

    error( ...
        'KAS must be positive.');

end


if KSE <= 0

    error( ...
        'KSE must be positive.');

end


if KAE < 0

    error( ...
        'KAE must be nonnegative.');

end


%% ========================================================================
% Effective branch conductance
%
% Symplastic hydraulic route consists of KAS and KSE in series.
% ========================================================================

KsymEffective = ...
    1 ...
    ./ ( ...
        1 ./ KAS ...
        + 1 ./ KSE);


%% ========================================================================
% Algebraic solution
%
% From:
%
%   JAS - JSE = dWSdt
%
%   JAE + JSE = E
%
%   JAS = KAS*(psiA - psiS)
%
%   JSE = KSE*(psiS - psiE)
%
%   JAE = KAE*(psiA - psiE)
%
%
% giving:
%
%   JSE =
%
%       E - (KAE/KAS)*dWSdt
%       ----------------------
%       1 + KAE/KAS + KAE/KSE
% ========================================================================

D = ...
    1 ...
    + KAE ./ KAS ...
    + KAE ./ KSE;


JSE = ...
    ( ...
        E ...
        - (KAE ./ KAS) .* dWSdt ...
    ) ...
    ./ D;


JAS = ...
    JSE ...
    + dWSdt;


JAE = ...
    E ...
    - JSE;


JXA = ...
    JAE ...
    + JAS;


%% ========================================================================
% Hydraulic potentials
% ========================================================================

psiA = ...
    psiS ...
    + JAS ./ KAS;


psiE = ...
    psiS ...
    - JSE ./ KSE;


psiX = ...
    psiA ...
    + JXA ./ KXA;


%% ========================================================================
% Direct constitutive-law reconstructions
% ========================================================================

JXA_check = ...
    KXA ...
    .* (psiX - psiA);


JAS_check = ...
    KAS ...
    .* (psiA - psiS);


JSE_check = ...
    KSE ...
    .* (psiS - psiE);


JAE_check = ...
    KAE ...
    .* (psiA - psiE);


%% ========================================================================
% Conservation diagnostics
% ========================================================================

apoplastClosure = ...
    JXA ...
    - JAE ...
    - JAS;


symplastClosure = ...
    JAS ...
    - JSE ...
    - dWSdt;


evaporationClosure = ...
    JAE ...
    + JSE ...
    - E;


wholeLeafClosure = ...
    JXA ...
    - E ...
    - dWSdt;


constitutiveResidualXA = ...
    JXA ...
    - JXA_check;


constitutiveResidualAS = ...
    JAS ...
    - JAS_check;


constitutiveResidualSE = ...
    JSE ...
    - JSE_check;


constitutiveResidualAE = ...
    JAE ...
    - JAE_check;


%% ========================================================================
% Advective partition diagnostics
%
% Under dynamic storage these are fractions of atmospheric water loss,
% not necessarily fractions of instantaneous xylem inflow.
%
% They may transiently fall outside [0,1] if storage charging/draining
% produces flow reversal. Do not clamp them.
% ========================================================================

if abs(E) > 1e-14

    fSymAdv = ...
        JSE ...
        ./ E;


    fApoAdv = ...
        JAE ...
        ./ E;

else

    fSymAdv = ...
        NaN;


    fApoAdv = ...
        NaN;

end


%% ========================================================================
% Steady-state branch fraction
%
% Diagnostic prediction for dWSdt = 0.
% ========================================================================

if ...
        KAE ...
        + KsymEffective > 0

    fSymSteady = ...
        KsymEffective ...
        ./ ( ...
            KAE ...
            + KsymEffective);

else

    fSymSteady = ...
        NaN;

end


%% ========================================================================
% Outputs
% ========================================================================

L = ...
    struct;


%% Potentials

L.psiXylemRequired = ...
    psiX;


L.psiApoplast = ...
    psiA;


L.psiSymplast = ...
    psiS;


L.psiEvaporationSite = ...
    psiE;


%% Fluxes

L.JxylemToApoplast = ...
    JXA;


L.JapoplastToSymplast = ...
    JAS;


L.JsymplastToEvaporation = ...
    JSE;


L.JapoplastToEvaporation = ...
    JAE;


L.symplasticStorageRate_mol_m2_s = ...
    dWSdt;


L.transpiration_mol_m2_s = ...
    E;


%% Conductance diagnostics

L.KxylemApoplast = ...
    KXA;


L.KapoplastSymplast = ...
    KAS;


L.KsymplastEvaporation = ...
    KSE;


L.KapoplastEvaporation = ...
    KAE;


L.KsymplasticBranchEffective = ...
    KsymEffective;


%% Pathway fractions

L.symplasticAdvectiveFraction = ...
    fSymAdv;


L.apoplasticAdvectiveFraction = ...
    fApoAdv;


L.symplasticAdvectiveFractionSteady = ...
    fSymSteady;


%% Conservation

L.apoplastClosure = ...
    apoplastClosure;


L.symplastClosure = ...
    symplastClosure;


L.evaporationClosure = ...
    evaporationClosure;


L.wholeLeafClosure = ...
    wholeLeafClosure;


%% Constitutive laws

L.constitutiveResidualXA = ...
    constitutiveResidualXA;


L.constitutiveResidualAS = ...
    constitutiveResidualAS;


L.constitutiveResidualSE = ...
    constitutiveResidualSE;


L.constitutiveResidualAE = ...
    constitutiveResidualAE;

end
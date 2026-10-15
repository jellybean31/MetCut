function B = ...
    prepare_root_xylem_isotope_boundary( ...
        Isoil, ...
        C, ...
        cfg)
%PREPARE_ROOT_XYLEM_ISOTOPE_BOUNDARY
%
% Convert interval-resolved common-root isotope composition from the
% resolved soil isotope model into the node representation expected by
% xylem_transport_dynamic.
%
% Interval k applies over:
%
%   time(k) -> time(k+1)
%
% and is therefore assigned to the ENDING node k+1. Xylem transport
% should use 'next' interpolation for this forcing.


if nargin < 3 || isempty(cfg)

    cfg = ...
        struct;

end


if ~isfield(cfg,'Rstandard_DH')

    cfg.Rstandard_DH = ...
        155.76e-6;

end


Rstd = ...
    cfg.Rstandard_DH;


deltaSoil = ...
    Isoil.deltaD_permil;


xSoil = ...
    Isoil.deuteriumAtomFraction;


[nCell,nT] = ...
    size(deltaSoil);


deltaInterval = ...
    Isoil.rootNodeDeltaDInterval_permil(:);


if numel(deltaInterval) ~= nT-1

    error( ...
        'Root-node interval isotope forcing has wrong length.');

end


if ...
        size(C.rootFluxGround_mol_m2_s,1) ~= nCell ...
        || size(C.rootFluxGround_mol_m2_s,2) ~= nT

    error( ...
        'Resolved root-flux dimensions do not match soil isotope state.');

end


%% ========================================================================
% Initial common-root composition
% ========================================================================

J0 = ...
    C.rootFluxGround_mol_m2_s(:,1);


uptake0 = ...
    J0 > 0;


if any(uptake0)

    xRoot0 = ...
        sum( ...
            J0(uptake0) ...
            .* xSoil(uptake0,1)) ...
        ./ sum(J0(uptake0));


    Rroot0 = ...
        xRoot0 ...
        ./ (1 - xRoot0);


    deltaRoot0 = ...
        1000 ...
        .* ( ...
            Rroot0 ...
            ./ Rstd ...
            - 1);

else

    firstFinite = ...
        find( ...
            isfinite(deltaInterval), ...
            1, ...
            'first');


    if isempty(firstFinite)

        error( ...
            'No finite root isotope source is available.');

    end


    deltaRoot0 = ...
        deltaInterval(firstFinite);

end


%% ========================================================================
% Ending-node representation
% ========================================================================

deltaNode = ...
    NaN(nT,1);


deltaNode(1) = ...
    deltaRoot0;


deltaNode(2:end) = ...
    deltaInterval;


%% ------------------------------------------------------------------------
% Isotope composition is irrelevant when root axial flow is exactly zero.
% Carry the previous finite composition through such nodes so the xylem
% ODE always receives a defined boundary value.
% -------------------------------------------------------------------------

for k = 2:nT

    if ~isfinite(deltaNode(k))

        deltaNode(k) = ...
            deltaNode(k-1);

    end

end


B.timeDays = ...
    Isoil.timeDays;


B.deltaDNode_permil = ...
    deltaNode;


B.deltaDInterval_permil = ...
    deltaInterval;


B.initialDeltaD_permil = ...
    deltaRoot0;


B.interpolationMethod = ...
    'next';


B.definition = ...
    ['Common-root isotope composition generated from resolved ' ...
     'positive soil-to-root uptake.'];

end

function C = couple_richards_soil_plant_hydraulics( ...
    timeDays, ...
    psiInitialCell_MPa, ...
    psiInitialGuess_MPa, ...
    precip_mm_d, ...
    met, ...
    rootTraitsGrid, ...
    plantState, ...
    xylemTraits, ...
    leafTraits, ...
    soilHydroGrid, ...
    G, ...
    richardsCfg, ...
    couplingCfg, ...
    initialState)
%COUPLE_RICHARDS_SOIL_PLANT_HYDRAULICS
%
% Fixed-point coupling between:
%
%   vertically resolved 1-D Richards soil hydraulics
%
% and
%
%   whole-plant hydraulics with one root pathway per Richards cell.
%
%
% SOIL/ROOT SPATIAL COUPLING
%
%   soil cell i <-> root pathway i
%
% The Richards grid may contain many hydraulic cells within each original
% physiological rooting zone. Each cell retains its own soil water
% potential and root pathway.
%
%
% AREA BASIS
%
% solve_plant_hydraulics:
%
%   rootFlux = mol H2O m^-2 LEAF s^-1
%
% Richards soil solver:
%
%   root sink/source = mol H2O m^-2 GROUND s^-1
%
% Therefore:
%
%   Jroot_ground = LAI * Jroot_leaf
%
% solve_plant_hydraulics reports both quantities. This wrapper
% verifies it before passing root flux to Richards.
%
%
% TIME BASIS
%
% Plant root fluxes are node values. Richards root sinks/sources are
% interval forcings. The wrapper therefore uses trapezoidal interval
% averages:
%
%   Jroot_interval(k) = 0.5 * [Jroot(k) + Jroot(k+1)]



%% ========================================================================
% Configuration
% ========================================================================

if nargin < 13 || isempty(couplingCfg)

    couplingCfg = ...
        struct;

end


couplingCfg = ...
    setDefault( ...
        couplingCfg, ...
        'maxIterations', ...
        30);

couplingCfg = ...
    setDefault( ...
        couplingCfg, ...
        'relaxation', ...
        0.4);

couplingCfg = ...
    setDefault( ...
        couplingCfg, ...
        'psiTolerance_MPa', ...
        1e-6);

couplingCfg = ...
    setDefault( ...
        couplingCfg, ...
        'rootFluxTolerance_mol_m2_leaf_s', ...
        1e-8);

couplingCfg = ...
    setDefault( ...
        couplingCfg, ...
        'maxAttempts', ...
        60);

couplingCfg = ...
    setDefault( ...
        couplingCfg, ...
        'minimumRelaxation', ...
        0.00625);

couplingCfg = ...
    setDefault( ...
        couplingCfg, ...
        'areaBasisTolerance_mol_m2_s', ...
        1e-12);


maxIterations = ...
    couplingCfg.maxIterations;

relaxation = ...
    couplingCfg.relaxation;

psiTolerance = ...
    couplingCfg.psiTolerance_MPa;

rootFluxTolerance = ...
    couplingCfg.rootFluxTolerance_mol_m2_leaf_s;

maxAttempts = ...
    couplingCfg.maxAttempts;

minimumRelaxation = ...
    couplingCfg.minimumRelaxation;

areaBasisTolerance = ...
    couplingCfg.areaBasisTolerance_mol_m2_s;

if nargin < 14 || isempty(initialState)

    initialState = ...
        struct;

end

%% ========================================================================
% Dimensions
% ========================================================================

t = ...
    timeDays(:).';

nTime = ...
    numel(t);


if nTime < 2 || any(diff(t) <= 0)

    error( ...
        'timeDays must contain at least two increasing values.');

end


if ~isfield(G,'nCell')

    error( ...
        'G.nCell is required.');

end


nCell = ...
    G.nCell;


psiInitialCell = ...
    psiInitialCell_MPa(:);


if numel(psiInitialCell) ~= nCell

    error( ...
        'psiInitialCell_MPa must contain one value per Richards cell.');

end


if any(~isfinite(psiInitialCell))

    error( ...
        'psiInitialCell_MPa must contain finite values.');

end


psiGuess = ...
    psiInitialGuess_MPa;


if isvector(psiGuess) && numel(psiGuess) == nCell

    psiGuess = ...
        repmat( ...
            psiGuess(:), ...
            1, ...
            nTime);

end


if ~isequal(size(psiGuess),[nCell nTime])

    error( ...
        ['psiInitialGuess_MPa must be nCell x nTime or contain ' ...
         'one value per Richards cell.']);

end


if any(~isfinite(psiGuess),'all')

    error( ...
        'psiInitialGuess_MPa contains nonfinite values.');

end


% Preserve exact initial condition.

psiGuess(:,1) = ...
    psiInitialCell;


%% ========================================================================
% Area-basis conversion
% ========================================================================

if ~isfield(plantState,'leafAreaIndex')

    error( ...
        'plantState.leafAreaIndex is required.');

end


LAI = ...
    plantState.leafAreaIndex;


if ...
        ~isscalar(LAI) ...
        || ~isfinite(LAI) ...
        || LAI <= 0

    error( ...
        'leafAreaIndex must be a positive scalar.');

end


%% ========================================================================
% Iteration diagnostics
% ========================================================================

psiResidualHistory = ...
    NaN(1,maxIterations);

rootFluxResidualHistory = ...
    NaN(1,maxIterations);

rootExtractionHistory_mm = ...
    NaN(1,maxIterations);

rootGrossUptakeHistory_mm = ...
    NaN(1,maxIterations);

rootGrossReleaseHistory_mm = ...
    NaN(1,maxIterations);

areaBasisResidualHistory = ...
    NaN(1,maxIterations);

rootFluxPrevious = ...
    NaN(nCell,nTime);

converged = ...
    false;


%% ========================================================================
% Adaptive fixed-point control
% ========================================================================

acceptedIterations = ...
    0;

attemptNumber = ...
    0;

relaxationHistory = ...
    NaN(1,maxAttempts);

failedAttempt = ...
    false(1,maxAttempts);


% Last successfully evaluated fixed-point map.

psiBase = ...
    psiGuess;

dPsiBase = ...
    NaN(size(psiGuess));

havePreviousMap = ...
    false;


%% ========================================================================
% Accepted-state bookkeeping
% ========================================================================

Haccepted = [];
Saccepted = [];

psiEvaluatedAccepted = [];

JrootLeafAccepted = [];
JrootGroundNodeAccepted = [];
JrootGroundIntervalAccepted = [];

maxAreaBasisResidualAccepted = NaN;


%% ========================================================================
% Fixed-point iteration
% ========================================================================

while ...
        acceptedIterations < maxIterations ...
        && attemptNumber < maxAttempts

    attemptNumber = ...
        attemptNumber + 1;


    relaxationHistory(attemptNumber) = ...
        relaxation;


    %% --------------------------------------------------------------------
    % Plant-facing resolved soil object
    % ---------------------------------------------------------------------

    soilPlant = ...
        struct;

    soilPlant.psiMPa = ...
        psiGuess;

    soilPlant.hydro = ...
        soilHydroGrid;


    %% --------------------------------------------------------------------
    % Whole-plant hydraulics
    % ---------------------------------------------------------------------

    H = ...
        solve_plant_hydraulics( ...
        t, ...
        soilPlant, ...
        met, ...
        rootTraitsGrid, ...
        plantState, ...
        xylemTraits, ...
        leafTraits, ...
        initialState);


    %% --------------------------------------------------------------------
    % Failed hydraulic attempt
    % ---------------------------------------------------------------------

    if any(H.solverFallback)

        failedAttempt(attemptNumber) = ...
            true;


        failureIndex = ...
            find( ...
                H.solverFallback, ...
                1, ...
                'first');


        fprintf('\n');
        fprintf('------------------------------------------------------------\n');
        fprintf('RICHARDS COUPLING ATTEMPT FAILED\n');
        fprintf('------------------------------------------------------------\n');
        fprintf('Attempt:                       %d\n', ...
            attemptNumber);
        fprintf('Accepted iterations so far:    %d\n', ...
            acceptedIterations);
        fprintf('Relaxation used:               %.6f\n', ...
            relaxation);
        fprintf('Failure cell trajectory index: %d\n', ...
            failureIndex);
        fprintf('Failure time:                  %.6f d\n', ...
            t(failureIndex));


        if ~havePreviousMap

            error([ ...
                'The hydraulic solver did not locate a valid solution for ' ...
                'the initial resolved soil coupling trajectory. Adaptive ' ...
                'relaxation cannot be applied until one successful map ' ...
                'has been evaluated.']);

        end


        relaxation = ...
            0.5 ...
            .* relaxation;


        if relaxation < minimumRelaxation

            %% ------------------------------------------------------------
            % J2b terminal coupling-feasibility diagnostic
            % ------------------------------------------------------------
            %
            % The last successful fixed-point map is
            %
            %   psi(alpha) = psiBase + alpha*dPsiBase
            %
            % and the ordinary retry logic has just shown that the next
            % allowed alpha would fall below minimumRelaxation. Diagnose
            % whether hydraulic feasibility is lost continuously along this
            % fixed-point direction by testing smaller alpha values without
            % rerunning Richards.

            fprintf('\n');
            fprintf('============================================================\n');
            fprintf('TERMINAL COUPLING-DIRECTION FEASIBILITY AUDIT\n');
            fprintf('============================================================\n');
            fprintf('Failure time:                  %.6f d\n', ...
                t(failureIndex));
            fprintf('Accepted iterations:           %d\n', ...
                acceptedIterations);
            fprintf('Last attempted relaxation:     %.8f\n', ...
                2 .* relaxation);
            fprintf('Next relaxation below minimum: %.8f\n', ...
                relaxation);
            fprintf('Minimum configured relaxation: %.8f\n', ...
                minimumRelaxation);

            alphaAudit = ...
                unique( ...
                    [ ...
                    0, ...
                    minimumRelaxation .* [0.125 0.25 0.5 0.75 1.0]]);

            nAlphaAudit = ...
                numel(alphaAudit);

            auditFeasible = ...
                false(1,nAlphaAudit);

            auditFirstFailureIndex = ...
                NaN(1,nAlphaAudit);

            auditMaxDeltaPsi_MPa = ...
                NaN(1,nAlphaAudit);

            auditFailureTime_days = ...
                NaN(1,nAlphaAudit);

            for iAudit = 1:nAlphaAudit

                alphaNow = ...
                    alphaAudit(iAudit);

                psiAudit = ...
                    psiBase;

                psiAudit(:,2:end) = ...
                    psiBase(:,2:end) ...
                    + alphaNow ...
                    .* dPsiBase(:,2:end);

                psiAudit(:,1) = ...
                    psiInitialCell;

                soilAudit = ...
                    struct;

                soilAudit.psiMPa = ...
                    psiAudit;

                soilAudit.hydro = ...
                    soilHydroGrid;

                Haudit = ...
                    solve_plant_hydraulics( ...
                    t, ...
                    soilAudit, ...
                    met, ...
                    rootTraitsGrid, ...
                    plantState, ...
                    xylemTraits, ...
                    leafTraits, ...
                    initialState);

                auditFeasible(iAudit) = ...
                    ~any(Haudit.solverFallback);

                if auditFeasible(iAudit)

                    auditFirstFailureIndex(iAudit) = ...
                        NaN;

                    auditFailureTime_days(iAudit) = ...
                        NaN;

                else

                    kAuditFail = ...
                        find( ...
                            Haudit.solverFallback, ...
                            1, ...
                            'first');

                    auditFirstFailureIndex(iAudit) = ...
                        kAuditFail;

                    auditFailureTime_days(iAudit) = ...
                        t(kAuditFail);

                end

                auditMaxDeltaPsi_MPa(iAudit) = ...
                    max( ...
                        abs( ...
                            psiAudit ...
                            - psiBase), ...
                        [], ...
                        'all', ...
                        'omitnan');

                fprintf( ...
                    ['alpha %.8f | feasible %d | first fail idx %g | ' ...
                    'fail time %.6f d | max |dPsi from base| %.12g MPa\n'], ...
                    alphaNow, ...
                    auditFeasible(iAudit), ...
                    auditFirstFailureIndex(iAudit), ...
                    auditFailureTime_days(iAudit), ...
                    auditMaxDeltaPsi_MPa(iAudit));

            end

            % Report the trial-soil movement at the original failure time.
            fprintf('\n');
            fprintf('FAILURE-NODE SOIL-TRAJECTORY MOVEMENT\n');
            fprintf('------------------------------------------------------------\n');

            dPsiFailNode = ...
                dPsiBase(:,failureIndex);

            fprintf('At time %.6f d:\n',t(failureIndex));
            fprintf('  min dPsiBase across cells: %.12g MPa\n', ...
                min(dPsiFailNode));
            fprintf('  max dPsiBase across cells: %.12g MPa\n', ...
                max(dPsiFailNode));
            fprintf('  max |dPsiBase| across cells: %.12g MPa\n', ...
                max(abs(dPsiFailNode)));

            [~,iLargestMove] = ...
                max(abs(dPsiFailNode));

            fprintf('  largest-moving soil cell:   %d\n', ...
                iLargestMove);
            fprintf('  psiBase at that cell:       %.12g MPa\n', ...
                psiBase(iLargestMove,failureIndex));
            fprintf('  dPsiBase at that cell:      %.12g MPa\n', ...
                dPsiBase(iLargestMove,failureIndex));


            %% ------------------------------------------------------------
            % Terminal residual localization / closure audit
            % ------------------------------------------------------------

            fprintf('\n');
            fprintf('TERMINAL FIXED-POINT RESIDUAL LOCALIZATION AUDIT\n');
            fprintf('------------------------------------------------------------\n');

            dPsiFail = ...
                dPsiBase(:,failureIndex);

            absDPsiFail = ...
                abs(dPsiFail);

            [absSorted,sortIdx] = ...
                sort(absDPsiFail,'descend');

            nTop = ...
                min(10,numel(absSorted));

            fprintf('Coupling psi tolerance:        %.12g MPa\n', ...
                psiTolerance);
            fprintf('Failure-node max |dPsi|:       %.12g MPa\n', ...
                max(absDPsiFail));
            fprintf('Failure-node median |dPsi|:    %.12g MPa\n', ...
                median(absDPsiFail));
            fprintf('Failure-node mean |dPsi|:      %.12g MPa\n', ...
                mean(absDPsiFail));

            qLevels = [0.50 0.75 0.90 0.95 0.99];
            qVals = quantile(absDPsiFail,qLevels);

            fprintf('Failure-node |dPsi| quantiles:\n');
            for iq = 1:numel(qLevels)
                fprintf('  q%02.0f: %.12g MPa\n', ...
                    100*qLevels(iq),qVals(iq));
            end

            fprintf('Cells above tolerance:         %d / %d\n', ...
                sum(absDPsiFail > psiTolerance),nCell);
            fprintf('Cells above 0.5*tolerance:     %d / %d\n', ...
                sum(absDPsiFail > 0.5*psiTolerance),nCell);
            fprintf('Cells above 0.25*tolerance:    %d / %d\n', ...
                sum(absDPsiFail > 0.25*psiTolerance),nCell);

            fprintf('\nTop residual cells at failure node:\n');
            fprintf('  cell    dPsiBase_MPa      |dPsi|_MPa       psiBase_MPa       psiCalc_MPa\n');

            for itop = 1:nTop
                icell = sortIdx(itop);
                fprintf('  %4d  %+ .12e  %.12e  %+ .12e  %+ .12e\n', ...
                    icell, ...
                    dPsiFail(icell), ...
                    absDPsiFail(icell), ...
                    psiBase(icell,failureIndex), ...
                    psiBase(icell,failureIndex)+dPsiFail(icell));
            end

            fprintf('\nLast accepted coupling diagnostics:\n');

            if acceptedIterations >= 1
                lastPsiResidual = ...
                    psiResidualHistory(acceptedIterations);

                lastRootFluxResidual = ...
                    rootFluxResidualHistory(acceptedIterations);

                lastAreaBasisResidual = ...
                    areaBasisResidualHistory(acceptedIterations);

                fprintf('  psi residual max:           %.12g MPa\n', ...
                    lastPsiResidual);
                fprintf('  psi residual / tolerance:   %.6f\n', ...
                    lastPsiResidual/psiTolerance);
                fprintf('  root-flux residual max:     %.12g mol m^-2 leaf s^-1\n', ...
                    lastRootFluxResidual);
                fprintf('  root-flux tolerance:        %.12g mol m^-2 leaf s^-1\n', ...
                    rootFluxTolerance);
                fprintf('  root-flux residual/tol:     %.6f\n', ...
                    lastRootFluxResidual/rootFluxTolerance);
                fprintf('  area-basis residual max:    %.12g mol m^-2 s^-1\n', ...
                    lastAreaBasisResidual);
            end

            % Print known plant closure diagnostics if present.
            fprintf('\nPlant closure diagnostics from last feasible map:\n');

            plantClosureNames = { ...
                'stemFlowResidual', ...
                'wholePlantStorageFluxClosureResidual', ...
                'massBalanceResidual', ...
                'rootFlowResidual', ...
                'leafFlowResidual'};

            anyPlantClosure = false;

            for ic = 1:numel(plantClosureNames)
                nm = plantClosureNames{ic};

                if isfield(Haccepted,nm)

                    x = Haccepted.(nm);

                    if isnumeric(x) && ~isempty(x)

                        anyPlantClosure = true;

                        fprintf('  %-38s max|.| = %.12g\n', ...
                            nm, ...
                            max(abs(x),[],'all','omitnan'));

                    end
                end
            end

            if ~anyPlantClosure
                fprintf('  No recognized plant closure fields found.\n');
            end

            % Surface scalar soil diagnostics whose names indicate mass,
            % closure, residual, storage, or water-balance information.
            fprintf('\nSoil/Richards scalar diagnostics from last feasible map:\n');

            soilNames = fieldnames(Saccepted);
            anySoilDiag = false;

            for isf = 1:numel(soilNames)

                nm = soilNames{isf};
                nmLower = lower(nm);

                isRelevant = ...
                    contains(nmLower,'residual') ...
                    || contains(nmLower,'closure') ...
                    || contains(nmLower,'mass') ...
                    || contains(nmLower,'storage') ...
                    || contains(nmLower,'water');

                if isRelevant

                    x = Saccepted.(nm);

                    if isnumeric(x) && isscalar(x) && isfinite(x)

                        anySoilDiag = true;

                        fprintf('  %-38s %.12g\n',nm,x);

                    end
                end
            end

            if ~anySoilDiag
                fprintf('  No recognized scalar soil closure fields found.\n');
            end

            % If a theta trajectory is stored, summarize how localized the
            % terminal soil-state change is in water-content space.
            thetaField = '';

            thetaCandidates = {'theta','thetaVol','thetaVolumetric','thetaCell'};

            for itc = 1:numel(thetaCandidates)
                if isfield(Saccepted,thetaCandidates{itc})
                    thetaField = thetaCandidates{itc};
                    break
                end
            end

            if ~isempty(thetaField)

                thetaStored = Saccepted.(thetaField);

                if isnumeric(thetaStored) ...
                        && isequal(size(thetaStored),[nCell nTime])

                    thetaFail = thetaStored(:,failureIndex);

                    fprintf('\nStored theta state at failure node (%s):\n',thetaField);
                    fprintf('  min theta:                  %.12g\n', ...
                        min(thetaFail,[],'omitnan'));
                    fprintf('  max theta:                  %.12g\n', ...
                        max(thetaFail,[],'omitnan'));
                    fprintf('  mean theta:                 %.12g\n', ...
                        mean(thetaFail,'omitnan'));
                end
            end

            fprintf('END TERMINAL RESIDUAL LOCALIZATION AUDIT\n');

            fprintf('\n');
            %% J2B_FULL_MAP_BRANCH_EXISTENCE_AUDIT
            %% Recompute the complete plant -> Richards map along the
            %% terminal fixed-point direction. Diagnostic only.
            fprintf('\n');
            fprintf('FULL-MAP LOCAL BRANCH-EXISTENCE AUDIT\n');
            fprintf('------------------------------------------------------------\n');

            alphaBranch=linspace(0,2*minimumRelaxation,17);
            nBranch=numel(alphaBranch);
            branchFeasible=false(1,nBranch);
            branchResidualMax=NaN(1,nBranch);
            branchResidualRMS=NaN(1,nBranch);
            branchProjection=NaN(1,nBranch);

            denomProjection=sum(dPsiBase(:).^2,'omitnan');

            for iBranch=1:nBranch
                a=alphaBranch(iBranch);
                psiBranch=psiBase;
                psiBranch(:,2:end)=psiBase(:,2:end)+a.*dPsiBase(:,2:end);
                psiBranch(:,1)=psiInitialCell;

                soilBranch=struct;
                soilBranch.psiMPa=psiBranch;
                soilBranch.hydro=soilHydroGrid;

                Hb=solve_plant_hydraulics( ...
                    t,soilBranch,met,rootTraitsGrid,plantState, ...
                    xylemTraits,leafTraits,initialState);

                branchFeasible(iBranch)=~any(Hb.solverFallback);

                if branchFeasible(iBranch)
                    Jb=Hb.rootFluxGround_mol_m2_s;
                    JbInt=0.5.*(Jb(:,1:end-1)+Jb(:,2:end));
                    Sb=soil_water_balance_richards_1d_prescribed_root( ...
                        t,precip_mm_d,JbInt,psiInitialCell, ...
                        soilHydroGrid,G,richardsCfg);
                    rb=Sb.psiMPa-psiBranch;
                    branchResidualMax(iBranch)=max(abs(rb),[],'all','omitnan');
                    branchResidualRMS(iBranch)=sqrt(mean(rb(:).^2,'omitnan'));
                    if denomProjection>0
                        branchProjection(iBranch)= ...
                            sum(rb(:).*dPsiBase(:),'omitnan')/denomProjection;
                    end
                end

                fprintf(['alpha %.9g | feasible %d | max residual %.12g MPa | ' ...
                    'RMS %.12g MPa | projection %.12g\n'], ...
                    a,branchFeasible(iBranch),branchResidualMax(iBranch), ...
                    branchResidualRMS(iBranch),branchProjection(iBranch));
            end

            finiteBranch=find(branchFeasible & isfinite(branchResidualMax));
            if ~isempty(finiteBranch)
                [minBranchResidual,jMinRel]=min(branchResidualMax(finiteBranch));
                jMin=finiteBranch(jMinRel);
                fprintf('Minimum feasible max residual: %.12g MPa at alpha %.9g\n', ...
                    minBranchResidual,alphaBranch(jMin));
                fprintf('Residual / production tolerance: %.6f\n', ...
                    minBranchResidual/psiTolerance);
                lastFeasible=finiteBranch(end);
                fprintf('Last feasible sampled alpha: %.9g\n',alphaBranch(lastFeasible));
                if lastFeasible<nBranch
                    fprintf('First infeasible sampled alpha: %.9g\n', ...
                        alphaBranch(lastFeasible+1));
                end
            end
            fprintf('END FULL-MAP LOCAL BRANCH-EXISTENCE AUDIT\n');

            fprintf('AUDIT INTERPRETATION FLAGS\n');
            fprintf('------------------------------------------------------------\n');
            fprintf('alpha=0 feasible:              %d\n', ...
                auditFeasible(1));
            fprintf('minimum alpha tested feasible: %d\n', ...
                any(auditFeasible(2:end)));
            fprintf('all positive alphas infeasible:%d\n', ...
                ~any(auditFeasible(2:end)));
            fprintf('============================================================\n');

            error( ...
                ['Hydraulic feasibility could not be recovered before ' ...
                 'relaxation fell below %.6f. terminal coupling-' ...
                 'direction audit printed above.'], ...
                minimumRelaxation);

        end


        psiGuess = ...
            psiBase;

        psiGuess(:,2:end) = ...
            psiBase(:,2:end) ...
            + relaxation ...
            .* dPsiBase(:,2:end);

        psiGuess(:,1) = ...
            psiInitialCell;


        fprintf('Retrying with relaxation:      %.6f\n', ...
            relaxation);


        continue

    end


    %% --------------------------------------------------------------------
    % Verify leaf-area -> ground-area conversion
    % ---------------------------------------------------------------------

    JrootLeaf = ...
        H.rootFlux;


    if ~isequal(size(JrootLeaf),[nCell nTime])

        error( ...
            'solve_plant_hydraulics returned an unexpected rootFlux size.');

    end


    if ~isfield(H,'rootFluxGround_mol_m2_s')

        error( ...
            ['solve_plant_hydraulics must return ' ...
             'rootFluxGround_mol_m2_s.']);

    end


    JrootGroundNode = ...
        H.rootFluxGround_mol_m2_s;


    areaBasisResidual = ...
        JrootGroundNode ...
        - LAI ...
        .* JrootLeaf;


    maxAreaBasisResidual = ...
        max( ...
            abs(areaBasisResidual), ...
            [], ...
            'all', ...
            'omitnan');


    if maxAreaBasisResidual > areaBasisTolerance

        error( ...
            ['Root-flux area-basis conversion failed: maximum residual ' ...
             '= %.12g mol m^-2 s^-1.'], ...
            maxAreaBasisResidual);

    end


    %% --------------------------------------------------------------------
    % Convert plant node fluxes to Richards interval forcing
    % ---------------------------------------------------------------------

    JrootGroundInterval = ...
        0.5 ...
        .* ( ...
            JrootGroundNode(:,1:end-1) ...
            + JrootGroundNode(:,2:end));


    %% --------------------------------------------------------------------
    % Resolved Richards soil water balance
    % ---------------------------------------------------------------------

    S = ...
        soil_water_balance_richards_1d_prescribed_root( ...
            t, ...
            precip_mm_d, ...
            JrootGroundInterval, ...
            psiInitialCell, ...
            soilHydroGrid, ...
            G, ...
            richardsCfg);


    if ~isequal(size(S.psiMPa),[nCell nTime])

        error( ...
            'Richards solver returned an unexpected psiMPa size.');

    end


    psiCalculated = ...
        S.psiMPa;


    %% --------------------------------------------------------------------
    % Fixed-point residual
    % ---------------------------------------------------------------------

    dPsi = ...
        psiCalculated ...
        - psiGuess;


    psiResidualMax = ...
        max( ...
            abs(dPsi), ...
            [], ...
            'all', ...
            'omitnan');


    %% --------------------------------------------------------------------
    % Successful map evaluation
    % ---------------------------------------------------------------------

    acceptedIterations = ...
        acceptedIterations + 1;

    acceptedIndex = ...
        acceptedIterations;


    psiResidualHistory(acceptedIndex) = ...
        psiResidualMax;

    areaBasisResidualHistory(acceptedIndex) = ...
        maxAreaBasisResidual;


    %% --------------------------------------------------------------------
    % Root-flux convergence
    % ---------------------------------------------------------------------

    if acceptedIndex == 1

        rootFluxResidualMax = ...
            Inf;

    else

        rootFluxResidualMax = ...
            max( ...
                abs( ...
                    JrootLeaf ...
                    - rootFluxPrevious), ...
                [], ...
                'all', ...
                'omitnan');

    end


    rootFluxResidualHistory(acceptedIndex) = ...
        rootFluxResidualMax;


    %% --------------------------------------------------------------------
    % Integrated root exchange on the exact forcing passed to Richards
    % ---------------------------------------------------------------------

    dtSeconds = ...
        diff(t) ...
        .* 86400;


    molarMassWater_kg_mol = ...
        18.01528e-3;


    JrootUptakeInterval = ...
        max(JrootGroundInterval,0);

    JrootReleaseInterval = ...
        max(-JrootGroundInterval,0);


    rootGrossUptakeHistory_mm(acceptedIndex) = ...
        sum( ...
            sum(JrootUptakeInterval,1) ...
            .* dtSeconds) ...
        .* molarMassWater_kg_mol;


    rootGrossReleaseHistory_mm(acceptedIndex) = ...
        sum( ...
            sum(JrootReleaseInterval,1) ...
            .* dtSeconds) ...
        .* molarMassWater_kg_mol;


    rootExtractionHistory_mm(acceptedIndex) = ...
        sum( ...
            sum(JrootGroundInterval,1) ...
            .* dtSeconds) ...
        .* molarMassWater_kg_mol;


    %% --------------------------------------------------------------------
    % Save successfully evaluated state
    % ---------------------------------------------------------------------

    Haccepted = ...
        H;

    Saccepted = ...
        S;

    psiEvaluatedAccepted = ...
        psiGuess;

    JrootLeafAccepted = ...
        JrootLeaf;

    JrootGroundNodeAccepted = ...
        JrootGroundNode;

    JrootGroundIntervalAccepted = ...
        JrootGroundInterval;

    maxAreaBasisResidualAccepted = ...
        maxAreaBasisResidual;


    %% --------------------------------------------------------------------
    % Progress
    % ---------------------------------------------------------------------

    fprintf( ...
        ['Richards coupling iteration %2d: ' ...
         'max dPsi = %.6g MPa, ' ...
         'max dJroot = %.6g, ' ...
         'net extraction = %.4f mm, ' ...
         'gross uptake = %.4f mm, ' ...
         'root release = %.4f mm, ' ...
         'relaxation = %.4f\n'], ...
        acceptedIndex, ...
        psiResidualMax, ...
        rootFluxResidualMax, ...
        rootExtractionHistory_mm(acceptedIndex), ...
        rootGrossUptakeHistory_mm(acceptedIndex), ...
        rootGrossReleaseHistory_mm(acceptedIndex), ...
        relaxation);


    %% --------------------------------------------------------------------
    % Convergence
    % ---------------------------------------------------------------------

    if ...
            psiResidualMax <= psiTolerance ...
            && rootFluxResidualMax <= rootFluxTolerance

        converged = ...
            true;

        break

    end


    %% --------------------------------------------------------------------
    % Store successful fixed-point map
    % ---------------------------------------------------------------------

    psiBase = ...
        psiGuess;

    dPsiBase = ...
        dPsi;

    havePreviousMap = ...
        true;

    rootFluxPrevious = ...
        JrootLeaf;


    %% --------------------------------------------------------------------
    % Construct next under-relaxed trajectory
    % ---------------------------------------------------------------------

    psiGuess = ...
        psiBase;

    psiGuess(:,2:end) = ...
        psiBase(:,2:end) ...
        + relaxation ...
        .* dPsiBase(:,2:end);

    psiGuess(:,1) = ...
        psiInitialCell;

end


%% ========================================================================
% Accepted state
% ========================================================================

if isempty(Haccepted) || isempty(Saccepted)

    error( ...
        'No successful resolved soil/plant hydraulic evaluation completed.');

end


if ~converged

    warning( ...
        ['Richards soil/plant coupling did not converge after %d accepted ' ...
         'iterations and %d total attempts.'], ...
        acceptedIterations, ...
        attemptNumber);

end


H = ...
    Haccepted;

S = ...
    Saccepted;


finalPsiResidual = ...
    S.psiMPa ...
    - psiEvaluatedAccepted;


%% ========================================================================
% Outputs
% ========================================================================

C.converged = ...
    converged;

C.nIterations = ...
    acceptedIterations;

C.nAttempts = ...
    attemptNumber;

C.nSoilCells = ...
    nCell;

C.grid = ...
    G;

C.hydro = ...
    H;

C.soil = ...
    S;

C.psiSoilCells_MPa = ...
    S.psiMPa;

C.rootFluxLeaf_mol_m2_s = ...
    JrootLeafAccepted;

C.rootFluxGround_mol_m2_s = ...
    JrootGroundNodeAccepted;

C.rootFluxGroundInterval_mol_m2_s = ...
    JrootGroundIntervalAccepted;

C.leafAreaIndex = ...
    LAI;

C.maxAreaBasisResidual_mol_m2_s = ...
    maxAreaBasisResidualAccepted;

C.maxFinalPsiResidual_MPa = ...
    max( ...
        abs(finalPsiResidual), ...
        [], ...
        'all', ...
        'omitnan');

C.finalPsiResidual_MPa = ...
    finalPsiResidual;

C.psiResidualHistory_MPa = ...
    psiResidualHistory(1:acceptedIterations);

C.rootFluxResidualHistory_mol_m2_leaf_s = ...
    rootFluxResidualHistory(1:acceptedIterations);

C.areaBasisResidualHistory_mol_m2_s = ...
    areaBasisResidualHistory(1:acceptedIterations);

C.rootExtractionHistory_mm = ...
    rootExtractionHistory_mm(1:acceptedIterations);

C.rootNetExtractionHistory_mm = ...
    rootExtractionHistory_mm(1:acceptedIterations);

C.rootGrossUptakeHistory_mm = ...
    rootGrossUptakeHistory_mm(1:acceptedIterations);

C.rootGrossReleaseHistory_mm = ...
    rootGrossReleaseHistory_mm(1:acceptedIterations);

C.relaxationHistory = ...
    relaxationHistory(1:attemptNumber);

C.failedAttempt = ...
    failedAttempt(1:attemptNumber);

C.finalRelaxation = ...
    relaxation;

C.psiPlantInputAccepted_MPa = ...
    psiEvaluatedAccepted;

C.areaBasisPlant = ...
    'per leaf area';

C.areaBasisSoil = ...
    'per ground area';

C.rootFluxTemporalBasisPlant = ...
    'node values';

C.rootFluxTemporalBasisSoil = ...
    'trapezoidal interval averages';

C.config = ...
    couplingCfg;

end


%% ========================================================================
function cfg = setDefault(cfg,name,value)

if ~isfield(cfg,name) || isempty(cfg.(name))

    cfg.(name) = ...
        value;

end

end

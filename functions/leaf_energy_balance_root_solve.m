function out = leaf_energy_balance_root_solve( ...
    hydraulicE_mol_m2_s,alphaSurface,shortwaveIncident_W_m2, ...
    LWincident_W_m2,epsilonLeaf,heatArgs,varargin)
% LEAF_ENERGY_BALANCE_ROOT_SOLVE
%
% Bracketed steady-state Tleaf root solver.
%
% Solves:
%   leaf_energy_balance_residual(...).residual_W_m2 = 0
%
% using the production residual helper.
%
% REQUIRED INPUTS
% ---------------
% hydraulicE_mol_m2_s
% alphaSurface
% shortwaveIncident_W_m2
% LWincident_W_m2
% epsilonLeaf
% heatArgs
%
% OPTIONAL NAME-VALUE
% -------------------
% 'Bracket_C'          [lower upper], default [-60 70]
% 'ResidualTolerance'  W m^-2, default 1e-8
% 'TemperatureTolerance_C' default 1e-8


p=inputParser;
p.addParameter('Bracket_C',[-60 70],@(x) ...
    isnumeric(x) && isreal(x) && numel(x)==2 && all(isfinite(x)));
p.addParameter('ResidualTolerance',1e-8,@(x) ...
    isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
p.addParameter('TemperatureTolerance_C',1e-8,@(x) ...
    isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
p.parse(varargin{:});

bracket_C=double(p.Results.Bracket_C(:).');
resTol=double(p.Results.ResidualTolerance);
tempTol_C=double(p.Results.TemperatureTolerance_C);

if bracket_C(1)>=bracket_C(2)
    error('leaf_energy_balance_root_solve:InvalidBracket', ...
        'Bracket_C must satisfy lower < upper.');
end

if bracket_C(1)<-100 || bracket_C(2)>100
    error('leaf_energy_balance_root_solve:BracketOutsideLatentHeatContract', ...
        'Bracket_C must remain within [-100,100] degC.');
end

lower_C=bracket_C(1);
upper_C=bracket_C(2);
lower_K=lower_C+273.15;
upper_K=upper_C+273.15;

low=evalResidualLocal(lower_K);
high=evalResidualLocal(upper_K);

requireFiniteLocal(low.residual_W_m2,'LOWER');
requireFiniteLocal(high.residual_W_m2,'UPPER');

out=struct();
out.success=false;
out.code="UNSET";
out.bracket_C=bracket_C;
out.lowerResidual_W_m2=low.residual_W_m2;
out.upperResidual_W_m2=high.residual_W_m2;
out.residualTolerance_W_m2=resTol;
out.temperatureTolerance_C=tempTol_C;
out.endpointRoot=false;
out.refinementPerformed=false;
out.refinementIterations=0;
out.iterations=NaN;
out.functionCount=NaN;
out.exitflag=NaN;
out.rootSolvePerformed=false;
out.hydraulicEAuthoritative=true;
out.hydraulicFeedback=false;

if abs(low.residual_W_m2)<=resTol
    final=low;
    root_C=lower_C;
    root_K=lower_K;
    out.endpointRoot=true;
    out.exitflag=1;
    out.iterations=0;
    out.functionCount=2;

elseif abs(high.residual_W_m2)<=resTol
    final=high;
    root_C=upper_C;
    root_K=upper_K;
    out.endpointRoot=true;
    out.exitflag=1;
    out.iterations=0;
    out.functionCount=2;

else
    if sign(low.residual_W_m2)==sign(high.residual_W_m2)
        error('leaf_energy_balance_root_solve:RootNotBracketed', ...
            ['Residual does not change sign across bracket: ' ...
             'Rlow=%.12g, Rhigh=%.12g W m^-2.'], ...
            low.residual_W_m2,high.residual_W_m2);
    end

    opts=optimset('TolX',tempTol_C,'Display','off');
    f=@(Tc) residualScalarLocal(Tc+273.15);

    [root_C,~,exitflag,output]=fzero(f,[lower_C upper_C],opts);

    if ~(isfinite(root_C) && exitflag>0)
        error('leaf_energy_balance_root_solve:RootSolveFailure', ...
            'Bracketed fzero solve failed.');
    end

    root_K=root_C+273.15;
    final=evalResidualLocal(root_K);

    out.exitflag=exitflag;
    out.iterations=output.iterations;
    out.functionCount=output.funcCount;

    % Explicitly enforce the residual tolerance independently of fzero's
    % temperature-based stopping criterion.
    if abs(final.residual_W_m2)>resTol
        out.refinementPerformed=true;

        loC=lower_C;
        hiC=upper_C;
        flo=low.residual_W_m2;
        fhi=high.residual_W_m2;

        % First tighten the sign bracket around the fzero solution.
        fr=final.residual_W_m2;
        if sign(fr)==sign(flo)
            loC=root_C;
            flo=fr;
        elseif sign(fr)==sign(fhi)
            hiC=root_C;
            fhi=fr;
        else
            % This can only occur for a numerically exact zero.
            loC=root_C;
            hiC=root_C;
        end

        maxRefine=200;
        bestC=root_C;
        best=final;

        for k=1:maxRefine
            if abs(best.residual_W_m2)<=resTol
                break
            end

            midC=0.5*(loC+hiC);

            % If floating point can no longer create a distinct midpoint,
            % residual closure cannot be improved further.
            if midC==loC || midC==hiC
                break
            end

            mid=evalResidualLocal(midC+273.15);
            requireFiniteLocal(mid.residual_W_m2,'REFINEMENT');

            if abs(mid.residual_W_m2)<abs(best.residual_W_m2)
                bestC=midC;
                best=mid;
            end

            if abs(mid.residual_W_m2)<=resTol
                bestC=midC;
                best=mid;
                break
            end

            if sign(mid.residual_W_m2)==sign(flo)
                loC=midC;
                flo=mid.residual_W_m2;
            else
                hiC=midC;
                fhi=mid.residual_W_m2;
            end

            out.refinementIterations=k;
        end

        root_C=bestC;
        root_K=root_C+273.15;
        final=best;
    end
end

requireFiniteLocal(final.residual_W_m2,'FINAL');

if abs(final.residual_W_m2)>resTol
    error('leaf_energy_balance_root_solve:ResidualClosureFailure', ...
        ['Final residual %.12g W m^-2 exceeds tolerance %.12g after ' ...
         'bracketed refinement.'], ...
        final.residual_W_m2,resTol);
end

out.success=true;
out.code="OK";
out.rootSolvePerformed=true;
out.Tleaf_K=root_K;
out.Tleaf_C=root_C;
out.residual_W_m2=final.residual_W_m2;
out.energyBalance=final;

    function r=residualScalarLocal(Tk)
        q=evalResidualLocal(Tk);
        requireFiniteLocal(q.residual_W_m2,'ITERATION');
        r=q.residual_W_m2;
    end

    function q=evalResidualLocal(Tk)
        q=leaf_energy_balance_residual( ...
            Tk,hydraulicE_mol_m2_s,alphaSurface, ...
            shortwaveIncident_W_m2,LWincident_W_m2, ...
            epsilonLeaf,heatArgs);
    end
end

function requireFiniteLocal(x,location)
if ~(isscalar(x) && isfinite(x))
    error('leaf_energy_balance_root_solve:NonfiniteResidual', ...
        '%s residual must be finite.',location);
end
end

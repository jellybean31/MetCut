
function R = full_chain_wax_only_grid_recovery(Rpilot,options)
%FULL_CHAIN_WAX_ONLY_GRID_RECOVERY
% Downstream-only 3-parameter synthetic inverse.
%
% Parameters:
%   I01 fixed
%   I02 inferred through exact upstream isotope propagation
%   L1,L2 inferred through effective-5 biochemical mapping
%
% Observations:
%   wax C27/C29/C31/C33 only


if nargin<2, options=struct; end
if ~isfield(options,'I01Fixed'), options.I01Fixed=2.03690139; end
if ~isfield(options,'truthI02'), options.truthI02=5.38306961; end
if ~isfield(options,'nI02'), options.nI02=11; end
if ~isfield(options,'nL1'), options.nL1=21; end
if ~isfield(options,'nL2'), options.nL2=21; end
if ~isfield(options,'seed'), options.seed=20260828; end

Mref=Rpilot.truth.model;
baseSample=resolve_truth_gsa_sample(Rpilot);

B=effective5_latent2_basis(Rpilot.truth);

% Use admissible latent bounds already established by the reduction.
D=effective5_latent2_bounds(B);
L1Grid=linspace(D(1).low,D(1).high,options.nL1);
L2Grid=linspace(D(2).low,D(2).high,options.nL2);

I02Grid=logspace(log10(2),log10(20),options.nI02);

% Ensure exact truth is represented.
[~,q]=min(abs(log(I02Grid)-log(options.truthI02)));
I02Grid(q)=options.truthI02;
I02Grid=sort(I02Grid);

Tfinal=full_chain_wax_only_truth(Rpilot,"final",options.seed);
Ttime=full_chain_wax_only_truth(Rpilot,"timeResolved",options.seed);

Jfinal=NaN(options.nI02,options.nL1,options.nL2);
Jtime=NaN(options.nI02,options.nL1,options.nL2);
pass=false(options.nI02,1);

for i=1:options.nI02
    fprintf('Full-chain wax-only upstream I02 %d / %d | %.6g d\n', ...
        i,options.nI02,I02Grid(i));

    try
        U=run_full_chain_upstream_I02_candidate( ...
            Mref,baseSample,options.I01Fixed,I02Grid(i));
        pass(i)=true;
    catch ME
        warning('I02 candidate %.6g failed: %s',I02Grid(i),ME.message);
        continue
    end

    for a=1:options.nL1
        for b=1:options.nL2
            L=[L1Grid(a);L2Grid(b)];

            if ~effective5_latent2_is_feasible(L,B)
                continue
            end

            try
                Y=full_chain_wax_only_candidate(Mref,U,L,B);
            catch ME
                if i==find(pass,1,'first') && a==1 && b==1
                    rethrow(ME)
                else
                    continue
                end
            end

            rf=(Y(:,Tfinal.sampleIndex)-Tfinal.observed_permil)./Tfinal.sigma_permil;
            rt=(Y(:,Ttime.sampleIndex)-Ttime.observed_permil)./Ttime.sigma_permil;

            Jfinal(i,a,b)=sum(rf(:).^2);
            Jtime(i,a,b)=sum(rt(:).^2);
        end
    end
end

R=struct;
R.version='2026-08-28_full_chain_wax_only_inverse_v1';
R.options=options;
R.I02Grid=I02Grid;
R.L1Grid=L1Grid;
R.L2Grid=L2Grid;
R.passI02=pass;
R.truthFinal=Tfinal;
R.truthTimeResolved=Ttime;
R.basis=B;
R.final=analyzeObjective(Jfinal,I02Grid,L1Grid,L2Grid, ...
    options.truthI02,effective5_to_latent2( ...
    effective5_truth_vector(Rpilot.truth),B));
R.timeResolved=analyzeObjective(Jtime,I02Grid,L1Grid,L2Grid, ...
    options.truthI02,effective5_to_latent2( ...
    effective5_truth_vector(Rpilot.truth),B));

fprintf('\n');
fprintf('============================================================\n');
fprintf('FULL-CHAIN WAX-ONLY RECOVERY COMPARISON\n');
fprintf('============================================================\n');
disp(summaryTable(R));
fprintf('============================================================\n');

end

function A=analyzeObjective(J,gI,g1,g2,trueI,trueL)
[Jmin,ix]=min(J(:),[],'omitnan');
[ii,a,b]=ind2sub(size(J),ix);

A=struct;
A.J=J;
A.Jmin=Jmin;
A.best=[gI(ii);g1(a);g2(b)];
A.truth=[trueI;trueL(:)];

% nearest truth node, but record whether truth lies in domain.
[~,it]=min(abs(log(gI)-log(trueI)));
[~,i1]=min(abs(g1-trueL(1)));
[~,i2]=min(abs(g2-trueL(2)));

A.truthGridDeltaJ=J(it,i1,i2)-Jmin;

% 3-parameter joint likelihood thresholds.
A.mask68=(J-Jmin)<=3.53;
A.mask95=(J-Jmin)<=7.81;
A.truthInside68=A.mask68(it,i1,i2);
A.truthInside95=A.mask95(it,i1,i2);

% One-dimensional profile objectives.
pI=squeeze(min(min(J,[],3,'omitnan'),[],2,'omitnan'));
p1=squeeze(min(min(J,[],3,'omitnan'),[],1,'omitnan'));
p2=squeeze(min(min(J,[],2,'omitnan'),[],1,'omitnan'));

dI=pI-min(pI,[],'omitnan');
d1=p1-min(p1,[],'omitnan');
d2=p2-min(p2,[],'omitnan');

A.I02Profile68=gI(dI<=1);
A.I02Profile95=gI(dI<=3.84);
A.L1Profile68=g1(d1<=1);
A.L1Profile95=g1(d1<=3.84);
A.L2Profile68=g2(d2<=1);
A.L2Profile95=g2(d2<=3.84);
end

function T=summaryTable(R)
Experiment=["Final wax only";"Time-resolved wax"];
bestI=[R.final.best(1);R.timeResolved.best(1)];
bestL1=[R.final.best(2);R.timeResolved.best(2)];
bestL2=[R.final.best(3);R.timeResolved.best(3)];
deltaJ=[R.final.truthGridDeltaJ;R.timeResolved.truthGridDeltaJ];
inside95=[R.final.truthInside95;R.timeResolved.truthInside95];

T=table(Experiment,bestI,bestL1,bestL2,deltaJ,inside95, ...
    'VariableNames',{'Experiment','BestI02','BestL1','BestL2', ...
    'TruthGridDeltaJ','TruthInside95'});
end

function R = example_02_parameterized_sensitivity()
%EXAMPLE_02_PARAMETERIZED_SENSITIVITY
% Compare the parameterized reference with a precipitation-dD
% perturbation, reusing the same validated hydraulic reference.


thisFile = mfilename('fullpath');
examplesDir = fileparts(thisFile);
root = fileparts(examplesDir);

fnDir = fullfile(root,'functions');
refDir = fullfile(root,'data','reference');

originalPath = path;
cleanup = onCleanup(@() path(originalPath)); 
addpath(fnDir,'-begin');
rehash

hydPath = fullfile(refDir, ...
    'hydraulic_checkpoint_76cell_8day_hourly_corrected_mass_conservative_storagefixed.mat');

assert(isfile(hydPath),'Hydraulic reference file is missing.');

P = struct;
P.initialSoilDeltaDByZone_permil = [-60; -75; -90];
P.precipDeltaD_permil = -50;

P.leafAnatomyConditioning = struct( ...
    'interveinalDistance_m',0.0005, ...
    'mesophyllThickness_m',0.00015, ...
    'effectiveFlowAreaFraction',0.02, ...
    'flowFraction',1, ...
    'tortuosity',1.5, ...
    'lateralPathFraction',0.5, ...
    'verticalPathFraction',0.5);

P.biochemistryPolicy = struct;
P.biochemistryPolicy.parameterMode = 'coreMidpointDiagnostic';
P.biochemistryPolicy.routeFractions = struct('C16',0,'C18',1);
P.biochemistryPolicy.alkaneCarbonNumber = [27 29 31 33];

fprintf('\n============================================================\n');
fprintf('METCUT EXAMPLE 02 - PARAMETERIZED SENSITIVITY\n');
fprintf('============================================================\n');

[Mref,Bref] = run_full_model_from_hydraulic_checkpoint_parameterized( ...
    hydPath,P);

Ppert = P;
Ppert.precipDeltaD_permil = -70;

[Mpert,Bpert] = run_full_model_from_hydraulic_checkpoint_parameterized( ...
    hydPath,Ppert);

Wref = localWaxMatrix(Mref);
Wpert = localWaxMatrix(Mpert);

finalRef = localLastFinite(Wref,2);
finalPert = localLastFinite(Wpert,2);

fprintf('Reference precipitation dD: %.3f per mil\n',P.precipDeltaD_permil);
fprintf('Perturbed precipitation dD: %.3f per mil\n',Ppert.precipDeltaD_permil);
fprintf('Reference final C29:        %.6f per mil\n',finalRef);
fprintf('Perturbed final C29:        %.6f per mil\n',finalPert);
fprintf('Final C29 difference:       %.6f per mil\n',finalPert-finalRef);

figure('Name','METCUT precipitation-dD sensitivity');
plot(Mref.timeDays,Wref(2,:),'DisplayName', ...
    sprintf('Precip %.0f‰',P.precipDeltaD_permil))
hold on
plot(Mpert.timeDays,Wpert(2,:),'DisplayName', ...
    sprintf('Precip %.0f‰',Ppert.precipDeltaD_permil))
hold off
xlabel('Time (days)')
ylabel('C29 instantaneous plant-wax \deltaD (‰)')
legend('Location','best')
title('Parameterized precipitation-isotope sensitivity')

R = struct;
R.referenceParameters = P;
R.perturbedParameters = Ppert;
R.referenceModel = Mref;
R.perturbedModel = Mpert;
R.referenceBudget = Bref;
R.perturbedBudget = Bpert;
R.finalC29Difference_permil = finalPert-finalRef;
end

function W = localWaxMatrix(M)
W = double(M.biochemistry.wax.deltaDMatrix);
if size(W,1) ~= 4 && size(W,2) == 4
    W = W.';
end
assert(size(W,1)==4,'Expected four plant-wax homologues.');
end

function x = localLastFinite(W,row)
v = W(row,:);
j = find(isfinite(v),1,'last');
assert(~isempty(j),'No finite plant-wax value found.');
x = v(j);
end

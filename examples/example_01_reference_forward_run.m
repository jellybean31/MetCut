function R = example_01_reference_forward_run()
%EXAMPLE_01_REFERENCE_FORWARD_RUN
% Run the METCUT reference forward model and inspect major states.


thisFile = mfilename('fullpath');
examplesDir = fileparts(thisFile);
root = fileparts(examplesDir);

fnDir = fullfile(root,'functions');
refDir = fullfile(root,'data','reference');

assert(isfolder(fnDir),'Could not locate METCUT_v1/functions.');
assert(isfolder(refDir),'Could not locate METCUT_v1/data/reference.');

originalPath = path;
cleanup = onCleanup(@() path(originalPath)); 
addpath(fnDir,'-begin');
rehash

hydPath = fullfile(refDir, ...
    'hydraulic_checkpoint_76cell_8day_hourly_corrected_mass_conservative_storagefixed.mat');

assert(isfile(hydPath),'Hydraulic reference file is missing.');

fprintf('\n============================================================\n');
fprintf('METCUT EXAMPLE 01 - REFERENCE FORWARD RUN\n');
fprintf('============================================================\n');

[M,B] = run_full_model_from_hydraulic_checkpoint(hydPath);

t = M.timeDays(:);
rootD = M.rootBoundary.deltaDNode_permil(:);
xylemD = M.xylemIsotopes.leafBoundary.deltaDWater(:);
leafD = M.leafIsotopes.deltaDSymplastBulk(:);
cyD = M.intracellular.deltaD.cytosol(:);
plD = M.intracellular.deltaD.plastid(:);

W = double(M.biochemistry.wax.deltaDMatrix);

if isfield(M.biochemistry.wax,'alkaneCarbonNumber')
    C = double(M.biochemistry.wax.alkaneCarbonNumber(:).');
else
    C = [27 29 31 33];
end

if size(W,1) ~= numel(C) && size(W,2) == numel(C)
    W = W.';
end

assert(size(W,1)==numel(C),'Unexpected plant-wax matrix shape.');

finalWax = NaN(size(C));
for i = 1:numel(C)
    row = W(i,:);
    j = find(isfinite(row),1,'last');
    finalWax(i) = row(j);
end

fprintf('Final root-source dD:   %.6f per mil\n',rootD(end));
fprintf('Final xylem dD:         %.6f per mil\n',xylemD(end));
fprintf('Final bulk-leaf dD:     %.6f per mil\n',leafD(end));
fprintf('Final cytosol dD:       %.6f per mil\n',cyD(end));
fprintf('Final plastid dD:       %.6f per mil\n',plD(end));

disp(table(C(:),finalWax(:), ...
    'VariableNames',{'AlkaneC','FinalDeltaD_permil'}))

figure('Name','METCUT reference process propagation');
plot(t,rootD,'DisplayName','Root source')
hold on
plot(t,xylemD,'DisplayName','Xylem delivered')
plot(t,leafD,'DisplayName','Leaf symplast')
plot(t,cyD,'DisplayName','Cytosol')
plot(t,plD,'DisplayName','Plastid')
plot(t,W(2,:).','DisplayName','C29 instantaneous wax')
hold off
xlabel('Time (days)')
ylabel('\deltaD (‰)')
legend('Location','best')
title('METCUT reference process propagation')

R = struct;
R.model = M;
R.budget = B;
R.carbonNumbers = C;
R.finalWaxDeltaD_permil = finalWax;
end

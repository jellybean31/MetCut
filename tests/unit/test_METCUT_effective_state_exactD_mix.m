function R = test_METCUT_effective_state_exactD_mix()
%TEST_METCUT_EFFECTIVE_STATE_EXACTD_MIX

fprintf('\n============================================================\n');
fprintf('METCUT UNIT: EFFECTIVE-STATE EXACT-D MIXING\n');
fprintf('============================================================\n');

Rstd = 155.76e-6;
dS = -200;
dW = -20;

% Endpoint identity.
d0 = effective_state_exactD_mix(dS,dW,0);
d1 = effective_state_exactD_mix(dS,dW,1);

details.endpointSourceResidual_permil = abs(d0-dS);
details.endpointWaterResidual_permil = abs(d1-dW);

assert(details.endpointSourceResidual_permil < 1e-10, ...
    'f=0 did not return the source isotope state.');
assert(details.endpointWaterResidual_permil < 1e-10, ...
    'f=1 did not return the water isotope state.');

% Independent exact-D midpoint calculation in D atom-fraction space.
f = 0.5;
d50 = effective_state_exactD_mix(dS,dW,f);

xS = localDeltaToX(dS,Rstd);
xW = localDeltaToX(dW,Rstd);
xExpected = (1-f).*xS + f.*xW;
dExpected = localXToDelta(xExpected,Rstd);

details.midpointResidual_permil = abs(d50-dExpected);
assert(details.midpointResidual_permil < 1e-10, ...
    'Exact-D midpoint does not match independent atom-fraction mixing.');

% Monotonic interpolation over the complete fraction interval.
fractions = [0 0.25 0.5 0.75 1];
d = zeros(size(fractions));
for i = 1:numel(fractions)
    d(i) = effective_state_exactD_mix(dS,dW,fractions(i));
end

details.monotonic = all(diff(d) > 0);
assert(details.monotonic, ...
    'Effective-state exact-D mixture is not monotonic between endpoints.');

R = struct;
R.name = 'test_METCUT_effective_state_exactD_mix';
R.pass = true;
R.details = details;
R.midpoint_permil = d50;

fprintf('Source endpoint residual:        %.3e per mil\n',details.endpointSourceResidual_permil);
fprintf('Water endpoint residual:         %.3e per mil\n',details.endpointWaterResidual_permil);
fprintf('Independent midpoint residual:   %.3e per mil\n',details.midpointResidual_permil);
fprintf('Monotonic 0->1 interpolation:    %d\n',details.monotonic);
fprintf('PASS\n');
fprintf('============================================================\n');
end

function x = localDeltaToX(d,Rstd)
R = Rstd.*(1+d./1000);
x = R./(1+R);
end

function d = localXToDelta(x,Rstd)
R = x./(1-x);
d = 1000.*(R./Rstd-1);
end

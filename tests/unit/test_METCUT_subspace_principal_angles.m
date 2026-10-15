function R = test_METCUT_subspace_principal_angles()
%TEST_METCUT_SUBSPACE_PRINCIPAL_ANGLES

fprintf('\n============================================================\n');
fprintf('METCUT UNIT: SUBSPACE PRINCIPAL ANGLES\n');
fprintf('============================================================\n');

details = struct;

% 1. Identical two-dimensional subspaces.
A = [1 0;0 1;0 0;0 0;0 0];
a = subspace_principal_angles_deg(A,A);
details.identicalMaxAngle_deg = max(abs(a));
assert(details.identicalMaxAngle_deg < 1e-10, ...
    'Identical subspaces did not return zero principal angles.');

% 2. Same subspace under an internal basis rotation.
theta = 37*pi/180;
Q = [cos(theta) -sin(theta); sin(theta) cos(theta)];
B = A*Q;
b = subspace_principal_angles_deg(A,B);
details.rotatedBasisMaxAngle_deg = max(abs(b));
assert(details.rotatedBasisMaxAngle_deg < 1e-10, ...
    'Principal angles changed under an internal basis rotation.');

% 3. Orthogonal one-dimensional subspaces.
c = subspace_principal_angles_deg([1;0],[0;1]);
details.orthogonalResidual_deg = abs(c-90);
assert(details.orthogonalResidual_deg < 1e-10, ...
    'Orthogonal one-dimensional subspaces did not return 90 degrees.');

% 4. Known 30-degree one-dimensional angle.
v1 = [1;0];
v2 = [cosd(30);sind(30)];
d = subspace_principal_angles_deg(v1,v2);
details.known30Residual_deg = abs(d-30);
assert(details.known30Residual_deg < 1e-10, ...
    'Known 30-degree principal angle was not recovered.');

R = struct;
R.name = 'test_METCUT_subspace_principal_angles';
R.pass = true;
R.details = details;

fprintf('Identical-subspace max angle:    %.3e deg\n',details.identicalMaxAngle_deg);
fprintf('Rotated-basis max angle:         %.3e deg\n',details.rotatedBasisMaxAngle_deg);
fprintf('Orthogonal-angle residual:       %.3e deg\n',details.orthogonalResidual_deg);
fprintf('Known-30deg residual:            %.3e deg\n',details.known30Residual_deg);
fprintf('PASS\n');
fprintf('============================================================\n');
end

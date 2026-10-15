function J = effective5_finite_difference_jacobian(Truth,x,options)
%EFFECTIVE5_FINITE_DIFFERENCE_JACOBIAN Central finite-difference 4x5 Jacobian.
if nargin<3 || isempty(options), options=struct; end
if ~isfield(options,'stepFractionOfRange'), options.stepFractionOfRange=1e-4; end
if ~isfield(options,'minimumStep_permil'), options.minimumStep_permil=1e-3; end
D=effective5_parameter_bounds(Truth.model.intracellular);
x=x(:).';
y0=effective5_predict_final_homologs(Truth,x);
J=NaN(numel(y0),5);
for j=1:5
    rangej=D(j).high-D(j).low;
    h=max(options.stepFractionOfRange*rangej,options.minimumStep_permil);
    xp=x; xm=x;
    xp(j)=min(D(j).high,x(j)+h);
    xm(j)=max(D(j).low,x(j)-h);
    if xp(j)==xm(j), error('Finite-difference step collapsed for parameter %d.',j); end
    yp=effective5_predict_final_homologs(Truth,xp);
    ym=effective5_predict_final_homologs(Truth,xm);
    J(:,j)=(yp-ym)./(xp(j)-xm(j));
end
end

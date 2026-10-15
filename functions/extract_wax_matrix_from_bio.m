
function W = extract_wax_matrix_from_bio(Bio)
%EXTRACT_WAX_MATRIX_FROM_BIO Return C27/C29/C31/C33 x time matrix.

if isfield(Bio,'wax') && isfield(Bio.wax,'deltaDMatrix')
    Y=Bio.wax.deltaDMatrix;
elseif isfield(Bio,'biochemistry') && isfield(Bio.biochemistry,'wax') && ...
        isfield(Bio.biochemistry.wax,'deltaDMatrix')
    Y=Bio.biochemistry.wax.deltaDMatrix;
else
    error('Biochemical output does not contain wax.deltaDMatrix.');
end

if isfield(Bio,'alkaneCarbonNumber')
    C=Bio.alkaneCarbonNumber(:);
elseif isfield(Bio,'biochemistry') && isfield(Bio.biochemistry,'alkaneCarbonNumber')
    C=Bio.biochemistry.alkaneCarbonNumber(:);
else
    C=[27;29;31;33];
end

if size(Y,1)~=numel(C) && size(Y,2)==numel(C)
    Y=Y.';
end

target=[27;29;31;33];
idx=zeros(4,1);

for k=1:4
    q=find(C==target(k),1);
    if isempty(q)
        error('Missing wax homolog C%d.',target(k));
    end
    idx(k)=q;
end

W=Y(idx,:);
end

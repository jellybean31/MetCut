function M = attach_wax_accumulation_to_model(M,T,waxSource)
%ATTACH_WAX_ACCUMULATION_TO_MODEL
% Add standing wax-pool accumulation without altering instantaneous wax.
%
% Interface:
%
%   M2 = attach_wax_accumulation_to_model(M,T)
%
% uses architecture-v1:
%   M.biochemistry.wax
%
% Explicit architecture-v2 interface:
%
%   M2 = attach_wax_accumulation_to_model(M,T,'waxStageExpanded')
%
% uses:
%   M.biochemistry.waxStageExpanded.deltaDMatrix
%
% while resolving time from:
%   M.biochemistry.wax.timeDays
%
% INPUTS
% ------
% M
%   full-model structure.
%
% T
%   wax accumulation traits. Default:
%     wax_accumulation_default_traits()
%
% waxSource
%   'wax'               - architecture-v1 instantaneous wax (default)
%   'waxStageExpanded'  - architecture-v2 instantaneous wax
%
% OUTPUT
% ------
% Added/updated:
%   M.biochemistry.waxAccumulation


if nargin<2 || isempty(T)
    T=wax_accumulation_default_traits();
end

if nargin<3 || isempty(waxSource)
    waxSource='wax';
end

if ~(ischar(waxSource) || isstring(waxSource))
    error('waxSource must be a character vector or string scalar.');
end
waxSource=char(waxSource);

if ~isstruct(M) || ~isfield(M,'biochemistry')
    error('Missing M.biochemistry.');
end

B=M.biochemistry;

if ~isfield(B,'wax') || ~isstruct(B.wax) || ...
        ~isfield(B.wax,'timeDays')
    error('M.biochemistry.wax.timeDays is required.');
end

timeDays=B.wax.timeDays;

switch lower(waxSource)

    case 'wax'

        if ~isfield(B.wax,'deltaDMatrix')
            error('M.biochemistry.wax.deltaDMatrix is required.');
        end

        Y=B.wax.deltaDMatrix;

        if isfield(B,'alkaneCarbonNumber')
            C=B.alkaneCarbonNumber;
        elseif isfield(B.wax,'alkaneCarbonNumber')
            C=B.wax.alkaneCarbonNumber;
        else
            error('Could not resolve alkane carbon numbers for wax.');
        end

        sourceCanonical='wax';
        architectureName='wax-biochemistry-v1';

    case 'waxstageexpanded'

        if ~isfield(B,'waxStageExpanded') || ...
                ~isstruct(B.waxStageExpanded) || ...
                ~isfield(B.waxStageExpanded,'deltaDMatrix')
            error(['M.biochemistry.waxStageExpanded.deltaDMatrix is ', ...
                'required for waxSource=''waxStageExpanded''.']);
        end

        Y=B.waxStageExpanded.deltaDMatrix;

        if isfield(B.waxStageExpanded,'alkaneCarbonNumber')
            C=B.waxStageExpanded.alkaneCarbonNumber;
        elseif isfield(B,'alkaneCarbonNumber')
            C=B.alkaneCarbonNumber;
        else
            error('Could not resolve alkane carbon numbers for waxStageExpanded.');
        end

        sourceCanonical='waxStageExpanded';

        if isfield(B,'architecture') && isfield(B.architecture,'name')
            architectureName=B.architecture.name;
        else
            architectureName='wax-biochemistry-v2';
        end

    otherwise

        error(['waxSource must be ''wax'' or ''waxStageExpanded''. ', ...
            'Received: %s'],waxSource);
end

% Interface checks before accumulation.
C=C(:);

if size(Y,1)~=numel(C)
    error('Wax isotope row count does not match resolved carbon numbers.');
end

nTime=numel(timeDays);

if size(Y,2)~=nTime && size(Y,1)==nTime
    % Let the accumulation operator support transposition only when homolog
    % count still agrees after transpose; otherwise fail explicitly.
    if size(Y,2)~=numel(C)
        error('Wax isotope dimensions are incompatible with time and homologues.');
    end
elseif size(Y,2)~=nTime
    error('Wax isotope column count does not match wax.timeDays.');
end

A=wax_accumulation_dynamic_exactD( ...
    timeDays,Y,C,T);

A.source=struct;
A.source.field=sourceCanonical;
A.source.architecture=architectureName;
A.source.timeField='biochemistry.wax.timeDays';
A.source.instantaneousPreserved=true;

M.biochemistry.waxAccumulation=A;

end

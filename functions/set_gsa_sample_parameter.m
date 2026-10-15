
function sample = set_gsa_sample_parameter(sample,id,value)
%SET_GSA_SAMPLE_PARAMETER Set retained GSA parameter in scalar struct.

id=char(string(id));

if isstruct(sample) && isfield(sample,id)
    sample.(id)=value;
    return
end

error('Could not set parameter %s in GSA sample.',id);
end

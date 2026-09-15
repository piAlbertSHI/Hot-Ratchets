function state = unpack_restart_1d_v1(p,packet)
%SMA.UNPACK_RESTART_1D_V1 Validate and unpack a self-describing restart.

assert(isstruct(p) && isstruct(packet), ...
    'sma:unpack_restart_1d_v1:Input','Parameters and restart must be structs.');
required={'schema_id','model_case_id','parameter_card_name','parameter_signature','state_width', ...
    'feature_flags','state_vector','semantic_contract'};
for k=1:numel(required)
    if ~isfield(packet,required{k})
        error('sma:unpack_restart_1d_v1:MissingMetadata', ...
            'Restart packet is missing metadata field %s.',required{k});
    end
end
if ~strcmp(packet.schema_id,'sma.restart.1d.v1')
    error('sma:unpack_restart_1d_v1:Schema','Unsupported restart schema: %s.',packet.schema_id);
end
if packet.model_case_id~=p.case_id || ~strcmp(packet.parameter_card_name,p.name)
    error('sma:unpack_restart_1d_v1:CardMismatch', ...
        'Restart card identity does not match the receiving parameter card.');
end
if ~isequaln(packet.parameter_signature,p)
    error('sma:unpack_restart_1d_v1:ParameterMismatch', ...
        'Restart parameter signature does not match the receiving card.');
end
expected_flags=feature_flags(p);
if ~isequaln(packet.feature_flags,expected_flags)
    error('sma:unpack_restart_1d_v1:FeatureMismatch', ...
        'Restart mechanism flags do not match the receiving parameter card.');
end
if packet.state_width~=sma.state_width_1d(p) || ...
        numel(packet.state_vector)~=packet.state_width
    error('sma:unpack_restart_1d_v1:Width', ...
        'Restart payload width is inconsistent with the receiving card.');
end
if ~isfield(packet.semantic_contract,'lossless') || ...
        ~isequal(packet.semantic_contract.lossless,true)
    error('sma:unpack_restart_1d_v1:SemanticContract', ...
        'Restart packet does not declare a lossless semantic contract.');
end
state=sma.unpack_state_1d(p,packet.state_vector);
end

function flags=feature_flags(p)
flags=struct( ...
    'Y_activity_population_mix',logical(p.Y_activity_population_mix>0), ...
    'Y_activity_population_bins',p.Y_activity_population_bins, ...
    'enable_explicit_r_phase_state',logical(p.enable_explicit_r_phase_state), ...
    'zhang_direct_history',logical(p.zhang_direct_Yp_coupling>0 || ...
        p.zhang_direct_Omega_coupling>0));
end

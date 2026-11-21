function state_matrix = initialize_nodes_1d(p,node_count,temperature_K)
%SMA.INITIALIZE_NODES_1D Allocate one independent state per axial node.

if nargin < 3
    temperature_K = p.T0*ones(node_count,1);
end
temperature_K = temperature_K(:);
if isscalar(temperature_K)
    temperature_K = repmat(temperature_K,node_count,1);
end
assert(numel(temperature_K)==node_count,'Temperature width mismatch.');
state_matrix = zeros(node_count,sma.state_width_1d(p));
for i = 1:node_count
    state = sma.initial_state_1d(p);
    state.T = temperature_K(i);
    state_matrix(i,:) = sma.pack_state_1d(p,state).';
end
end

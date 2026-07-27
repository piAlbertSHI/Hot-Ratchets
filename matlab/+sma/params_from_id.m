function p = params_from_id(id)
%SMA.PARAMS_FROM_ID Return a parameter set from its stable numeric ID.

names = {'niti_zhang2022','cualmn_annealed','cualmn_drawn130', ...
    'cualmn_drawn80','nimnga_textured_review','qian2019_se21', ...
    'nimnga_review_fig4'};
index = round(id);
if index < 1 || index > numel(names)
    error('sma:params_from_id:UnknownID', 'Unknown case ID: %g', id);
end
p = sma.params(names{index});
end

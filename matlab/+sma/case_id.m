function id = case_id(name)
%SMA.CASE_ID Stable IDs for Simulink parameter selection.

switch lower(char(name))
    case {'niti_zhang2022','niti','niti_s4','niti_s6', ...
            'niti_s8','niti_s10','niti_s12'}
        id = 1;
    case {'cualmn_annealed','annealed'}
        id = 2;
    case {'cualmn_drawn130','drawn130'}
        id = 3;
    case {'cualmn_drawn80','drawn80'}
        id = 4;
    case {'nimnga_textured_review','nimnga','review'}
        id = 5;
    case {'qian2019_se21','se21'}
        id = 6;
    case {'nimnga_review_fig4','review_fig4'}
        id = 7;
    otherwise
        error('sma:case_id:UnknownCase', 'Unknown case: %s', char(name));
end
end

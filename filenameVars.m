function vars = filenameVars(constants, r_num, bc_type)

if nargin < 3
    bc_type = 'mushyLayer';
end

if strcmp(bc_type, 'mushyLayer')
    
    a_fixed = constants('a_fixed');
    
    vars = strcat('Points',num2str(r_num),...
        'Rm',num2str(constants('Rm')), ...
        'H',num2str(constants('H')), ...
        'R', num2str(constants('R')));
    
    if a_fixed > 0
        vars = strcat(vars, 'a', num2str(a_fixed));
    else
        vars = strcat('_relaxed_', vars);
    end
    
    vars = strcat(vars,'b', num2str(constants('b')));

elseif strcmp(bc_type, 'heatedWire')
    vars = strcat('Rm',num2str(constants('Rm')),'r_num',num2str(r_num));
end
    
    
end
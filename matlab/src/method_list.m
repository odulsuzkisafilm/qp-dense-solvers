function methods = method_list(names)

    methods = struct('name', {}, 'fn', {});
    for i = 1:numel(names)
        switch names{i}
            case 'QR'
                fn = @(H, C, g, d) two_output(@solve_qr_null, H, C, g, d);
            case 'Schur'
                fn = @(H, C, g, d) two_output(@solve_schur, H, C, g, d);
            case 'VR'
                fn = @solve_vr_null;
            case 'LDL'
                fn = @solve_kkt_ldl;
            otherwise
                error('Unknown method %s', names{i});
        end
        methods(end+1) = struct('name', names{i}, 'fn', fn);
    end
end

function [x, lambda, info] = two_output(solver, H, C, g, d)
    [x, lambda] = solver(H, C, g, d);
    info = struct();
end

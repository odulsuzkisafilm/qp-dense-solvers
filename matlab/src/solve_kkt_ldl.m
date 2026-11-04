function [x, lambda, info] = solve_kkt_ldl(H, C, g, d)

    qp_check_inputs(H, C, g, d);
    n = size(H, 1);
    p = size(C, 1);
    K = [H, C'; C, zeros(p)];
    kmax = max(abs(K(:)));
    [L, D, pv] = ldl(K, 'vector');
    clear K
    D = sparse(D);
    [inertia, pivot_min] = ldl_inertia(D);
    info = struct('inertia', inertia, 'pivot_min', pivot_min, ...
        'l_max', max(abs(L(:))), 'd_growth', max(abs(D(:)))/kmax);
    if any(inertia ~= [n, p, 0])
        error('qp:Inertia', 'KKT inertia is (%d,%d,%d), expected (%d,%d,%d).', ...
            inertia(1), inertia(2), inertia(3), n, p, 0);
    end
    b = [g; d];
    v = zeros(n + p, 1);
    v(pv) = L' \ (D \ (L \ b(pv)));
    x = v(1:n);
    lambda = v(n+1:end);
    if any(~isfinite(x)) || any(~isfinite(lambda))
        error('qp:NonfiniteSolution', 'Nonfinite LDL solution.');
    end
end

function out = kkt_ldl_probe(H, C, g, d, alpha, delta)

    n = size(H, 1);
    p = size(C, 1);
    K = [H + alpha*eye(n), C'; C, -delta*eye(p)];
    kmax = max(abs(K(:)));
    [L, D, pv] = ldl(K, 'vector');
    clear K
    D = sparse(D);
    [inertia, pivot_min] = ldl_inertia(D);
    b = [g; d];
    v = zeros(n + p, 1);
    v(pv) = L' \ (D \ (L \ b(pv)));
    out = struct('alpha', alpha, 'delta', delta, 'inertia', inertia, ...
        'x', v(1:n), 'lambda', v(n+1:end), 'pivot_min', pivot_min, ...
        'l_max', max(abs(L(:))), 'd_growth', max(abs(D(:)))/kmax);
end

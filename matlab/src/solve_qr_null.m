function [x, lambda] = solve_qr_null(H, C, g, d)

    qp_check_inputs(H, C, g, d);
    n = size(H, 1);
    p = size(C, 1);

    [Q, R_full] = qr(C');
    Q1 = Q(:, 1:p);
    Z = Q(:, p+1:n);
    R = R_full(1:p, :);
    clear Q R_full
    if any(diag(R) == 0)
        error('qp:ConstraintRank', 'QR produced an exactly zero diagonal in R.');
    end

    x0 = Q1*(R' \ d);

    HZ = H*Z;
    G = Z'*HZ;
    clear HZ
    G = (G+G')/2;
    rhs = Z'*(g-H*x0);
    [L, flag] = chol(G, 'lower');
    if flag ~= 0
        error('qp:ReducedCholesky', 'Reduced Hessian Cholesky failed at pivot %d.', flag);
    end
    clear G
    z = L' \ (L \ rhs);
    clear L rhs

    x = x0+Z*z;
    lambda = R \ (Q1'*(g-H*x));
    if any(~isfinite(x)) || any(~isfinite(lambda))
        error('qp:NonfiniteSolution', 'Nonfinite QR/null-space solution.');
    end
end

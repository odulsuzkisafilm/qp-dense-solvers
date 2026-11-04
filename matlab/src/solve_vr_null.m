function [x, lambda, info] = solve_vr_null(H, C, g, d)

    qp_check_inputs(H, C, g, d);
    n = size(H, 1);
    p = size(C, 1);

    [L, U, perm] = lu(C', 'vector');
    if any(diag(U) == 0)
        error('qp:ConstraintRank', 'LU produced an exactly zero pivot in U.');
    end
    i2 = perm(1:p);
    i1 = perm(p+1:n);
    L1 = L(1:p, :);
    L2 = L(p+1:n, :);
    clear L

    Y = -(L1' \ L2');
    clear L2
    info = struct('y_fro', norm(Y, 'fro'), 'y_max', max(abs(Y(:))), ...
        'u_min', min(abs(diag(U))), 'u_max', max(abs(diag(U))));

    x0 = zeros(n, 1);
    x0(i2) = L1' \ (U' \ d);

    HZ = H(:, i1) + H(:, i2)*Y;
    G = HZ(i1, :) + Y'*HZ(i2, :);
    clear HZ
    G = (G+G')/2;
    r = g - H*x0;
    rhs = r(i1) + Y'*r(i2);
    [R, flag] = chol(G, 'lower');
    if flag ~= 0
        error('qp:ReducedCholesky', 'Reduced Hessian Cholesky failed at pivot %d.', flag);
    end
    clear G
    z = R' \ (R \ rhs);
    clear R rhs

    x = x0;
    x(i1) = x(i1) + z;
    x(i2) = x(i2) + Y*z;
    r = g - H*x;
    lambda = U \ (L1 \ r(i2));
    if any(~isfinite(x)) || any(~isfinite(lambda))
        error('qp:NonfiniteSolution', 'Nonfinite variable reduction solution.');
    end
end

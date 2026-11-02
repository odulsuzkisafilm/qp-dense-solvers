function [x, lambda] = solve_schur(H, C, g, d)

    qp_check_inputs(H, C, g, d);
    [L_H, flag] = chol(H, 'lower');
    if flag ~= 0
        error('qp:HessianCholesky', 'Hessian Cholesky failed at pivot %d.', flag);
    end

    u = L_H' \ (L_H \ g);
    W = L_H' \ (L_H \ C');
    clear L_H

    S = C*W;
    S = (S+S')/2;
    [L_S, flag] = chol(S, 'lower');
    if flag ~= 0
        error('qp:SchurCholesky', 'Schur Cholesky failed at pivot %d.', flag);
    end
    clear S
    lambda = L_S' \ (L_S \ (C*u-d));

    x = u-W*lambda;
    if any(~isfinite(x)) || any(~isfinite(lambda))
        error('qp:NonfiniteSolution', 'Nonfinite Schur-complement solution.');
    end
end

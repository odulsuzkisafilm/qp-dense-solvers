function P = make_indefinite_case(n, p, seed, kind, kH, kC)

    if nargin < 5
        kH = 100;
    end
    if nargin < 6
        kC = 10;
    end
    rng(seed, 'twister');
    [Q, ~] = qr(randn(n));
    V = Q(:, 1:p);
    Z = Q(:, p+1:n);
    clear Q
    q = n - p;
    [W, ~] = qr(randn(p));
    C = (W.*logspace(0, -log10(kC), p))*V';
    clear W

    a = logspace(0, -log10(kH), p);
    a(1:2:end) = -a(1:2:end);
    b = logspace(0, -log10(kH), q);
    switch kind
        case 'reduced_pd'
            m = 0;
        case 'reduced_indef'
            m = max(1, round(q/10));
            b(end-m+1:end) = -b(end-m+1:end);
        otherwise
            error('Unknown kind');
    end
    E = 0.1*randn(p, q);
    H = V*(a'.*V') + Z*(b'.*Z') + V*E*Z' + Z*E'*V';
    H = (H + H')/2;
    clear V Z E

    x = randn(n, 1);
    lambda = randn(p, 1);
    P = struct('H', H, 'C', C, 'g', H*x + C'*lambda, 'd', C*x, ...
        'x_ref', x, 'lambda_ref', lambda, 'f_constant', 0, ...
        'case_id', sprintf('INDEF_%s_n%d_p%d_seed%d', kind, n, p, seed), ...
        'family', ['INDEF_' kind], 'category', 'synthetic_indefinite_hessian', 'seed', seed, ...
        'ref_kind', 'manufactured_target_not_certified', 'kH_target', kH, 'kC_target', kC, ...
        'kH_measured', NaN, 'kC_measured', NaN, ...
        'negative_reduced_eigenvalues', m, 'expected_inertia', [n - m, p + m, 0]);
end

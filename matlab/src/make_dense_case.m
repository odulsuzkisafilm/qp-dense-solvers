function P = make_dense_case(n, p, seed, kH, kC)

    if nargin<4
        kH = 100;
    end
    if nargin<5
        kC = 10;
    end
    rng(seed, 'twister');
    [U, ~] = qr(randn(n));
    eigenvalues = logspace(0, -log10(kH), n);
    H = (U.*eigenvalues)*U';
    H = (H+H')/2;
    clear U
    [V, ~] = qr(randn(n, p), 0);
    [W, ~] = qr(randn(p));
    C = (W.*logspace(0, -log10(kC), p))*V';
    clear V W
    x = randn(n, 1);
    lambda = randn(p, 1);
    P = struct('H', H, 'C', C, 'g', H*x+C'*lambda, 'd', C*x, ...
        'x_ref', x, 'lambda_ref', lambda, 'f_constant', 0, ...
        'case_id', sprintf('DENSE_n%d_p%d_seed%d', n, p, seed), ...
        'family', 'DENSE', 'category', 'synthetic_timing', 'seed', seed, ...
        'ref_kind', 'manufactured_target_not_certified', 'kH_target', kH, 'kC_target', kC, ...
        'kH_measured', NaN, 'kC_measured', NaN);
end

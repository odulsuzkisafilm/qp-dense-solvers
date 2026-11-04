function out = run_regularize(level)

    if nargin < 1
        level = 'check';
    end
    assert(any(strcmp(level, {'check', 'full'})), 'Unknown level');
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'));
    cfg = final_options(false);
    cfg.level = level;
    out = suite_folder(root, ['regularize_' level], cfg);
    rows = [];

    deltas = [0 1e-16 1e-14 1e-12 1e-10 1e-8 1e-6];
    files = dir(fullfile(root, 'data', 'COND_*.mat'));
    if strcmp(level, 'check')
        files = files(contains({files.name}, 'seed11'));
    end
    for k = 1:numel(files)
        P = load(fullfile(files(k).folder, files(k).name));
        n = size(P.H, 1);
        p = size(P.C, 1);
        [xq, lq] = solve_qr_null(P.H, P.C, P.g, P.d);
        for delta = deltas
            r = kkt_ldl_probe(P.H, P.C, P.g, P.d, 0, delta);
            rows = [rows; summarize(P, 'constraint_block', 0, delta, r, n, p, xq, lq)];
        end
        suite_checkpoint(rows, out);
    end

    seeds = [11 22 33];
    if strcmp(level, 'check')
        seeds = 11;
    end
    for seed = seeds
        P = make_indefinite_case(1000, 500, seed, 'reduced_indef');
        n = size(P.H, 1);
        p = size(P.C, 1);
        [Q, ~] = qr(P.C');
        Z = Q(:, p+1:n);
        G = Z'*P.H*Z;
        G = (G+G')/2;
        e = eig(G);
        lam_min = min(e);
        P.negative_reduced_eigenvalues_measured = sum(e < 0);
        for c = [0 0.5 0.9 1.1 2 10]
            alpha = c*abs(lam_min);
            r = kkt_ldl_probe(P.H, P.C, P.g, P.d, alpha, 0);
            xq = nan(n, 1);
            lq = nan(p, 1);
            if alpha > abs(lam_min)
                [xq, lq] = solve_qr_null(P.H + alpha*eye(n), P.C, P.g, P.d);
            end
            row = summarize(P, 'hessian_shift', alpha, 0, r, n, p, xq, lq);
            row.shift_over_lam_min = c;
            row.lam_min_reduced = lam_min;
            row.negative_reduced_measured = P.negative_reduced_eigenvalues_measured;
            rows = [rows; row];
        end
        suite_checkpoint(rows, out);
    end
end

function row = summarize(P, kind, alpha, delta, r, n, p, xq, lq)
    rs = P.H*r.x + P.C'*r.lambda - P.g;
    rf = P.C*r.x - P.d;
    Hfro = norm(P.H, 'fro');
    Cfro = norm(P.C, 'fro');
    row = struct('case_id', string(P.case_id), 'kind', string(kind), 'n', n, 'p', p, ...
        'kH_target', P.kH_target, 'kC_target', P.kC_target, 'alpha', alpha, 'delta', delta, ...
        'inertia', string(sprintf('%d,%d,%d', r.inertia)), 'inertia_ok', double(isequal(r.inertia, [n, p, 0])), ...
        'pivot_min', r.pivot_min, 'l_max', r.l_max, 'd_growth', r.d_growth, ...
        'x_reference_error', norm(r.x - P.x_ref)/norm(P.x_ref), ...
        'lambda_reference_error', norm(r.lambda - P.lambda_ref)/norm(P.lambda_ref), ...
        'x_vs_qr', norm(r.x - xq)/norm(xq), 'lambda_vs_qr', norm(r.lambda - lq)/norm(lq), ...
        'stationarity_scaled', norm(rs)/(Hfro*norm(r.x) + Cfro*norm(r.lambda) + norm(P.g)), ...
        'feasibility_scaled', norm(rf)/(Cfro*norm(r.x) + norm(P.d)), ...
        'shift_over_lam_min', NaN, 'lam_min_reduced', NaN, 'negative_reduced_measured', NaN);
end

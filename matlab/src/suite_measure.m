function rows = suite_measure(P, cfg)

    n = size(P.H, 1);
    p = size(P.C, 1);
    base = struct('case_id', string(P.case_id), 'family', string(P.family), ...
        'category', string(P.category), 'n', n, 'p', p, 'ratio', p/n, 'seed', P.seed, ...
        'kH_target', P.kH_target, 'kC_target', P.kC_target, ...
        'kH_measured', P.kH_measured, 'kC_measured', P.kC_measured, ...
        'method', "", 'repeat', 0, 'batch_count', 0, 'seconds', NaN, 'seconds_over_n3', NaN, ...
        'status', "", 'warning_id', "", 'warning_message', "", 'error_id', "", 'error_message', "", ...
        'x_reference_error', NaN, 'lambda_reference_error', NaN, 'reference_kind', string(P.ref_kind), ...
        'kkt_fro_scaled', NaN, 'feasibility_scaled', NaN, 'stationarity_scaled', NaN, ...
        'prediction_m', NaN, 'prediction_abs_error_m', NaN, 'weight_sum_error', NaN, ...
        'objective', NaN, 'peak_estimate_GiB', cfg.peak_dense_arrays*8*n^2/2^30);
    rows = repmat(base, 0, 1);
    if base.peak_estimate_GiB>cfg.memory_budget_GiB
        for name = ["QR", "Schur"]
            row = base;
            row.method = name;
            row.status = "memory_guard_skip";
            rows(end+1, 1) = row;
        end
        return
    end
    try
        H = full(P.H);
        C = full(P.C);
        g = full(P.g);
        d = full(P.d);
    catch ME
        for name = ["QR", "Schur"]
            row = base;
            row.method = name;
            row.status = "allocation_failed";
            row.error_id = string(ME.identifier);
            row.error_message = string(ME.message);
            rows(end+1, 1) = row;
        end
        return
    end
    Hfro = norm(H, 'fro');
    Cfro = norm(C, 'fro');
    Kfro = hypot(Hfro, sqrt(2)*Cfro);
    bnorm = norm([g;d]);
    fn = {@solve_qr_null, @solve_schur};
    names = ["QR", "Schur"];
    if cfg.warmup
        for j = 1:2
            try
                [~, ~] = fn{j}(H, C, g, d);
            catch
            end
        end
    end
    for r = 1:cfg.repeats
        order = 1:2;
        if mod(r+P.seed, 2) == 0
            order = 2:-1:1;
        end
        for j = order
            row = base;
            row.method = names(j);
            row.repeat = r;
            lastwarn('');
            count = 0;
            elapsed = 0;
            try
                t = tic;
                while count<cfg.max_batch_solves && (count == 0 || elapsed<cfg.min_batch_seconds)
                    [x, lambda] = fn{j}(H, C, g, d);
                    count = count+1;
                    elapsed = toc(t);
                end
                row.seconds = elapsed/count;
                row.seconds_over_n3 = row.seconds/n^3;
                row.batch_count = count;
                [wm, wi] = lastwarn;
                row.warning_id = string(wi);
                row.warning_message = string(wm);
                row.status = "ok";
                if ~isempty(wm)
                    row.status = "returned_warning";
                end
                rs = H*x+C'*lambda-g;
                rf = C*x-d;
                row.kkt_fro_scaled = norm([rs;rf])/max(Kfro*norm([x;lambda])+bnorm, realmin);
                row.feasibility_scaled = norm(rf)/max(Cfro*norm(x)+norm(d), realmin);
                row.stationarity_scaled = norm(rs)/max(Hfro*norm(x)+Cfro*norm(lambda)+norm(g), realmin);
                row.x_reference_error = norm(x-P.x_ref)/max(norm(P.x_ref), realmin);
                row.lambda_reference_error = norm(lambda-P.lambda_ref)/max(norm(P.lambda_ref), realmin);
                row.objective = .5*x'*(H*x)-g'*x+P.f_constant;
                if isfield(cfg, 'solution_dir')
                    save(fullfile(cfg.solution_dir, sprintf('%s_%s_r%d.mat', P.case_id, names(j), r)), 'x', 'lambda');
                end
                if isfield(P, 'observations')
                    row.prediction_m = P.observations'*x;
                    row.prediction_abs_error_m = abs(row.prediction_m-P.heldout_elevation);
                    row.weight_sum_error = abs(sum(x)-1);
                end
            catch ME
                row.status = "failed";
                row.error_id = string(ME.identifier);
                row.error_message = string(ME.message);
                row.batch_count = count;
                [wm, wi] = lastwarn;
                row.warning_id = string(wi);
                row.warning_message = string(wm);
            end
            rows(end+1, 1) = row;
            fprintf('%s n=%d p=%d %s repeat=%d: %s, %.4g s, x discrepancy %.3g\n', ...
                P.case_id, n, p, names(j), r, row.status, row.seconds, row.x_reference_error);
        end
    end
end

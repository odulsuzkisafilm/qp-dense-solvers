function out = run_kriging(level)

    if nargin<1
        level = 'check';
    end
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'));
    assert(any(strcmp(level, {'check', 'pilot', 'full'})), 'Unknown level');
    cfg = final_options(strcmp(level, 'full'));
    cfg.level = level;
    ns = 200;
    seeds = 11;
    if strcmp(level, 'pilot')
        ns = [2000 12000];
    end
    if strcmp(level, 'full')
        ns = [1000 2000 4000 8000 12000];
        seeds = [11 22 33 44 55];
    end
    out = suite_folder(root, ['kriging_' level], cfg);
    rows = [];
    cfg.solution_dir = fullfile(out, 'solutions');
    mkdir(cfg.solution_dir);
    for n = ns
        assert(cfg.peak_dense_arrays*8*n^2/2^30 <= cfg.memory_budget_GiB, 'Memory guard');
        for seed = seeds
            P = make_kriging_case(n, seed, root);
            meta = rmfield(P, {'H', 'C', 'g', 'd', 'x_ref', 'lambda_ref', 'observations'});
            save(fullfile(out, [P.case_id '_metadata.mat']), 'meta');
            if strcmp(level, 'check')
                v = [P.H P.C';P.C 0]\[P.g;P.d];
                P.x_ref = v(1:n);
                P.lambda_ref = v(end);
                P.ref_kind = 'independent_double_KKT_not_high_precision';
                P.kH_measured = cond(P.H);
                assert(norm(P.H-P.H', 'fro') == 0);
                chol(P.H);
            end
            new = suite_measure(P, cfg);
            rows = [rows;new];
            suite_checkpoint(rows, out);
            if strcmp(level, 'check')
                assert(all(string({new.status}) == "ok"), 'Solver check failed; inspect CSV');
                assert(all([new.x_reference_error]<1e-8), 'KKT reference disagreement');
                assert(all([new.weight_sum_error]<1e-8), 'Weights do not sum to one');
            end
            clear P
            if any(string({new.status}) ~= "ok")
                return
            end
        end
    end
end

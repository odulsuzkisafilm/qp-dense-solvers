function out = run_large_conditioning(n)

    if nargin<1
        n = 2000;
    end
    validateattributes(n, {'double'}, {'scalar', 'integer', '>=', 100});
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'));
    cfg = suite_options();
    cfg.n = n;
    cfg.repeats = 1;
    cfg.min_batch_seconds = 0;
    cfg.max_batch_solves = 1;
    cfg.warmup = false;
    assert(cfg.peak_dense_arrays*8*n^2/2^30 <= cfg.memory_budget_GiB, 'Memory guard: reduce n');
    out = suite_folder(root, sprintf('large_conditioning_n%d', n), cfg);
    rows = [];
    for family = ["H", "C"]
        for target = [1e8 1e10 1e12]
            for seed = [11 22 33]
                kh = 100;
                kc = 10;
                if family == "H"
                    kh = target;
                else
                    kc = target;
                end
                P = make_dense_case(n, round(n/4), seed, kh, kc);
                P.case_id = sprintf('LARGE_COND_%s_n%d_log%d_seed%d', family, n, round(log10(target)), seed);
                P.family = char("COND_"+family);
                P.category = 'synthetic_large_conditioning';
                new = suite_measure(P, cfg);
                rows = [rows;new];
                suite_checkpoint(rows, out);
                clear P new
            end
        end
    end
    plot_suite_results(out);
end

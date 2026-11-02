function out = run_dense_scaling(mode, level)

    if nargin<1
        mode = 'dimension';
    end
    if nargin<2
        level = 'pilot';
    end
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'));
    cfg = suite_options();
    cfg.mode = mode;
    cfg.level = level;
    assert(any(strcmp(mode, {'dimension', 'fraction'})), 'Unknown mode');
    assert(any(strcmp(level, {'pilot', 'full', 'largest'})), 'Unknown level');
    if strcmp(level, 'pilot')
        ns = [1000 2000];
        ratios = .5;
        seeds = 11;
        cfg.repeats = 3;
        if strcmp(mode, 'fraction')
            ns = 2000;
            ratios = [.1 .5 .9];
        end
    else
        ns = [1000 2000 4000 8000 10000];
        ratios = .5;
        seeds = [11 22 33];
        if strcmp(mode, 'fraction')
            ns = 8000;
            ratios = [.1 .25 .5 .75 .9];
        end
        if strcmp(level, 'largest')
            ns = 10000;
        end
    end
    out = suite_folder(root, ['dense_' mode '_' level], cfg);
    rows = [];
    for n = ns
        for ratio = ratios
            p = round(n*ratio);
            for seed = seeds
                estimate = cfg.peak_dense_arrays*8*n^2/2^30;
                assert(estimate <= cfg.memory_budget_GiB, 'Memory guard: reduce n before generation.');
                P = make_dense_case(n, p, seed);
                new = suite_measure(P, cfg);
                rows = [rows;new];
                suite_checkpoint(rows, out);
                clear P new
            end
        end
    end
    plot_suite_results(out);
end

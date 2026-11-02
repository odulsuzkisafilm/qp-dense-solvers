function out = run_literature(level)

    if nargin<1
        level = 'check';
    end
    assert(any(strcmp(level, {'check', 'pilot', 'largepilot', 'full'})), 'Unknown level');
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'));
    cfg = final_options(strcmp(level, 'full'));
    cfg.level = level;
    out = suite_folder(root, ['literature_' level], cfg);
    rows = [];
    cfg.solution_dir = fullfile(out, 'solutions');
    mkdir(cfg.solution_dir);
    if strcmp(level, 'check')
        for f = {'DTOC3_N10.mat', 'AUG2DC_NX3_NY3.mat'}
            P = load(fullfile(root, 'data', f{1}));
            P = orthogonal_benchmark(P, 11);
            new = suite_measure(P, cfg);
            rows = [rows;new];
            suite_checkpoint(rows, out);
            assert(all(string({new.status}) == "ok"), 'Check failed; inspect CSV');
            assert(all([new.x_reference_error]<1e-8), 'Reference disagreement');
        end
        return
    end
    files = {'LSE_gemat1_p5.mat', 'LSE_gemat1_p20.mat', 'DTOC3_N1000.mat'};
    if strcmp(level, 'largepilot')
        files = {'DTOC3_N4001.mat'};
    end
    if strcmp(level, 'full')
        files = [files {'DTOC3_N4001.mat'}];
    end
    for k = 1:numel(files)
        raw = load(fullfile(root, 'data', files{k}));
        if isfield(raw, 'A_ls')
            raw = rmfield(raw, {'A_ls', 'b_ls'});
        end
        n = size(raw.H, 1);
        assert(cfg.peak_dense_arrays*8*n^2/2^30 <= cfg.memory_budget_GiB, 'Memory guard');
        P = raw;
        new = suite_measure(P, cfg);
        rows = [rows;new];
        suite_checkpoint(rows, out);
        clear P
        if startsWith(files{k}, 'DTOC3')
            seeds = 11;
            if strcmp(level, 'full')
                seeds = [11 22 33 44 55];
            end
            for seed = seeds
                P = orthogonal_benchmark(raw, seed);
                metadata = struct('source', files{k}, 'seed', seed, ...
                    'quadratic_check', P.transform_quadratic_check, ...
                    'constraint_check', P.transform_constraint_check);
                save(fullfile(out, [P.case_id '_transformation.mat']), 'metadata');
                new = suite_measure(P, cfg);
                rows = [rows;new];
                suite_checkpoint(rows, out);
                clear P
                if any(string({new.status}) ~= "ok")
                    return
                end
            end
        end
        clear raw
    end
end

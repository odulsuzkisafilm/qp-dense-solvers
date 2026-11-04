function out = run_direct3(level)

    if nargin < 1
        level = 'check';
    end
    assert(any(strcmp(level, {'check', 'pilot', 'conditioning', 'dimension'})), 'Unknown level');
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'));
    cfg = final_options(strcmp(level, 'dimension'));
    cfg.level = level;
    cfg.methods = method_list({'QR', 'Schur', 'LDL'});
    out = suite_folder(root, ['direct3_' level], cfg);
    rows = [];
    cfg.solution_dir = fullfile(out, 'solutions');
    mkdir(cfg.solution_dir);

    switch level
        case 'check'
            for f = {'DTOC3_N10.mat', 'AUG2DC_NX3_NY3.mat', 'COND_H_log00_seed11.mat'}
                P = load(fullfile(root, 'data', f{1}));
                new = suite_measure(P, cfg);
                rows = [rows; new];
                suite_checkpoint(rows, out);
                assert(all(string({new.status}) == "ok"), 'Check failed');
                assert(all([new.x_reference_error] < 1e-8), 'Reference disagreement');
            end
            P = make_dense_case(300, 120, 11);
            new = suite_measure(P, cfg);
            rows = [rows; new];
            suite_checkpoint(rows, out);
            assert(all(string({new.status}) == "ok"), 'Dense check failed');
        case 'pilot'
            for n = [2000 6000]
                P = make_dense_case(n, n/2, 11);
                new = suite_measure(P, cfg);
                rows = [rows; new];
                suite_checkpoint(rows, out);
                clear P
            end
        case 'conditioning'
            cfg.methods = method_list({'LDL'});
            files = dir(fullfile(root, 'data', 'COND_*.mat'));
            for k = 1:numel(files)
                P = load(fullfile(files(k).folder, files(k).name));
                new = suite_measure(P, cfg);
                rows = [rows; new];
                suite_checkpoint(rows, out);
            end
        case 'dimension'
            cfg.methods = method_list({'LDL'});
            for n = [1000 2000 4000 8000]
                for seed = [11 22 33 44 55]
                    P = make_dense_case(n, n/2, seed);
                    new = suite_measure(P, cfg);
                    rows = [rows; new];
                    suite_checkpoint(rows, out);
                    clear P
                    if any(string({new.status}) ~= "ok")
                        return
                    end
                end
            end
    end
end

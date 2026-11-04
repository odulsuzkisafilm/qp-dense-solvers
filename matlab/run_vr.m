function out = run_vr(level)

    if nargin < 1
        level = 'check';
    end
    assert(any(strcmp(level, {'check', 'dimension', 'fraction', 'conditioning', 'illcond', 'applications'})), 'Unknown level');
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'));
    cfg = final_options(any(strcmp(level, {'dimension', 'fraction', 'applications'})));
    cfg.level = level;
    cfg.methods = method_list({'VR'});
    out = suite_folder(root, ['vr_' level], cfg);
    rows = [];
    cfg.solution_dir = fullfile(out, 'solutions');
    mkdir(cfg.solution_dir);

    switch level
        case 'check'
            cfg.methods = method_list({'QR', 'Schur', 'VR', 'LDL'});
            for f = {'DTOC3_N10.mat', 'AUG2DC_NX3_NY3.mat', 'COND_H_log00_seed11.mat', 'COND_C_log10_seed11.mat', 'LSE_gemat1_p5.mat'}
                P = load(fullfile(root, 'data', f{1}));
                if isfield(P, 'A_ls')
                    P = rmfield(P, {'A_ls', 'b_ls'});
                end
                new = suite_measure(P, cfg);
                rows = [rows; new];
                suite_checkpoint(rows, out);
            end
            P = make_dense_case(300, 120, 11);
            new = suite_measure(P, cfg);
            rows = [rows; new];
            suite_checkpoint(rows, out);
            vr = rows(string({rows.method}) == "VR");
            assert(all(string({vr.status}) == "ok"), 'VR check failed');
            assert(all([vr.x_reference_error] < 1e-6), 'VR reference disagreement');
        case 'dimension'
            for n = [1000 2000 4000 8000 12000]
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
        case 'fraction'
            n = 12000;
            for ratio = [.1 .25 .5 .75 .9]
                for seed = [11 22 33 44 55]
                    P = make_dense_case(n, round(n*ratio), seed);
                    new = suite_measure(P, cfg);
                    rows = [rows; new];
                    suite_checkpoint(rows, out);
                    clear P
                    if any(string({new.status}) ~= "ok")
                        return
                    end
                end
            end
        case 'conditioning'
            files = dir(fullfile(root, 'data', 'COND_*.mat'));
            for k = 1:numel(files)
                P = load(fullfile(files(k).folder, files(k).name));
                new = suite_measure(P, cfg);
                rows = [rows; new];
                suite_checkpoint(rows, out);
            end
        case 'illcond'
            for n = [1000 2000 4000 8000 12000]
                for seed = [11 22]
                    P = make_dense_case(n, n/2, seed, 1e8, 1e6);
                    P.case_id = sprintf('FINAL_HC_n%d_p%d_k1e+08_seed%d', n, n/2, seed);
                    P.family = 'FINAL_HC';
                    P.category = 'illcond_dimension_uncertified_targets';
                    new = suite_measure(P, cfg);
                    rows = [rows; new];
                    suite_checkpoint(rows, out);
                    clear P
                end
            end
        case 'applications'
            cfg.methods = method_list({'QR', 'Schur', 'VR'});
            for f = {'LSE_gemat1_p5.mat', 'LSE_gemat1_p20.mat', 'DTOC3_N4001.mat'}
                raw = load(fullfile(root, 'data', f{1}));
                if isfield(raw, 'A_ls')
                    raw = rmfield(raw, {'A_ls', 'b_ls'});
                end
                new = suite_measure(raw, cfg);
                rows = [rows; new];
                suite_checkpoint(rows, out);
                if startsWith(f{1}, 'DTOC3')
                    for seed = [11 22 33 44 55]
                        P = orthogonal_benchmark(raw, seed);
                        new = suite_measure(P, cfg);
                        rows = [rows; new];
                        suite_checkpoint(rows, out);
                        clear P
                    end
                end
                clear raw
            end
            for n = [4000 12000]
                P = make_kriging_case(n, 11, root);
                new = suite_measure(P, cfg);
                rows = [rows; new];
                suite_checkpoint(rows, out);
                clear P
            end
    end
end

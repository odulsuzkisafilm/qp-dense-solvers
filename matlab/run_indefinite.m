function out = run_indefinite(level)

    if nargin < 1
        level = 'check';
    end
    assert(any(strcmp(level, {'check', 'full'})), 'Unknown level');
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'));
    cfg = final_options(false);
    cfg.level = level;
    cfg.methods = method_list({'QR', 'Schur', 'LDL'});
    out = suite_folder(root, ['indefinite_' level], cfg);
    rows = [];
    cfg.solution_dir = fullfile(out, 'solutions');
    mkdir(cfg.solution_dir);

    if strcmp(level, 'check')
        ns = 300;
        seeds = 11;
    else
        ns = [1000 4000];
        seeds = [11 22 33];
    end
    for n = ns
        for kind = {'reduced_pd', 'reduced_indef'}
            for seed = seeds
                P = make_indefinite_case(n, n/2, seed, kind{1});
                new = suite_measure(P, cfg);
                rows = [rows; new];
                suite_checkpoint(rows, out);
                if strcmp(level, 'check')
                    st = string({new.status});
                    nm = string({new.method});
                    if strcmp(kind{1}, 'reduced_pd')
                        assert(all(st(nm == "QR" | nm == "LDL") == "ok"), 'QR and LDL should solve reduced_pd');
                        assert(all(st(nm == "Schur") == "failed"), 'Schur should fail on indefinite H');
                        assert(all([new(nm == "QR" | nm == "LDL").x_reference_error] < 1e-8), 'Reference disagreement');
                    else
                        assert(all(st(nm == "LDL") == "failed"), 'LDL should report wrong inertia');
                        assert(all(st(nm == "QR") == "failed"), 'QR should fail on indefinite reduced Hessian');
                    end
                end
                clear P
            end
        end
    end
end

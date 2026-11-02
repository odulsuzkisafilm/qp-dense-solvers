function out = run_final_dense(mode)

    if nargin<1
        mode = 'pilot';
    end
    assert(any(strcmp(mode, {'pilot', 'dimension', 'fraction', 'conditioning', 'illcond'})), 'Unknown mode');
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'));
    cfg = final_options(any(strcmp(mode, {'dimension', 'fraction'})));
    cfg.mode = mode;
    ns = 12000;
    ratios = .5;
    seeds = [11 22 33 44 55];
    families = "baseline";
    targets = 100;
    if strcmp(mode, 'pilot')
        seeds = 11;
    end
    if strcmp(mode, 'dimension')
        ns = [1000 2000 4000 8000 12000];
    end
    if strcmp(mode, 'fraction')
        ratios = [.1 .25 .5 .75 .9];
    end
    if strcmp(mode, 'conditioning')
        families = ["H", "C"];
        targets = [1e4 1e8 1e10 1e12];
    end

    if strcmp(mode, 'illcond')
        ns = [1000 2000 4000 8000 12000];
        families = "HC";
        targets = 1e8;
        seeds = [11 22];
    end
    out = suite_folder(root, ['final_dense_' mode], cfg);
    rows = [];
    cfg.solution_dir = fullfile(out, 'solutions');
    mkdir(cfg.solution_dir);
    for n = ns
        assert(cfg.peak_dense_arrays*8*n^2/2^30 <= cfg.memory_budget_GiB, 'Memory guard');
        for ratio = ratios
            for family = families
                for target = targets
                    for seed = seeds
                        kh = 100;
                        kc = 10;
                        if family == "H"
                            kh = target;
                        elseif family == "C"
                            kc = target;
                        elseif family == "HC"
                            kh = 1e8;
                            kc = 1e6;
                        end
                        P = make_dense_case(n, round(n*ratio), seed, kh, kc);
                        P.case_id = sprintf('FINAL_%s_n%d_p%d_k%.0e_seed%d', family, n, round(n*ratio), target, seed);
                        P.family = char("FINAL_"+family);
                        if strcmp(mode, 'conditioning')
                            P.category = 'large_conditioning_uncertified_targets';
                        end
                        if strcmp(mode, 'illcond')
                            P.category = 'illcond_dimension_uncertified_targets';
                        end
                        new = suite_measure(P, cfg);
                        rows = [rows;new];
                        suite_checkpoint(rows, out);
                        clear P
                        if ~any(strcmp(mode, {'conditioning', 'illcond'})) && any(string({new.status}) ~= "ok")
                            return
                        end
                    end
                end
            end
        end
    end
end

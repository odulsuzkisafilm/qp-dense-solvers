function out = run_cutest_tests(profile)

    if nargin<1
        profile = 'pilot';
    end
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'));
    cfg = suite_options();
    cfg.profile = profile;
    if strcmp(profile, 'pilot')
        cfg.repeats = 3;
    elseif ~strcmp(profile, 'large')
        error('Use pilot or large.');
    end
    out = suite_folder(root, ['cutest_' profile], cfg);
    rows = [];
    files = [dir(fullfile(root, 'data', 'DTOC3_*.mat'));dir(fullfile(root, 'data', 'AUG2DC_*.mat'))];
    assert(~isempty(files), 'CUTEst MAT files are missing.');
    sizes = zeros(numel(files), 1);
    for k = 1:numel(files)
        meta = whos('-file', fullfile(files(k).folder, files(k).name), 'H');
        sizes(k) = meta.size(1);
    end
    [~, idx] = sort(sizes);
    files = files(idx);
    for k = 1:numel(files)
        P = load(fullfile(files(k).folder, files(k).name));
        if strcmp(profile, 'pilot') && size(P.H, 1)>500
            continue;
        end
        new = suite_measure(P, cfg);
        rows = [rows;new];
        suite_checkpoint(rows, out);
        clear P new
    end
    plot_suite_results(out);
end

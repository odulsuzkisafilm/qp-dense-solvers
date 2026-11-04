function out = run_conditioning_tests()

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'));
    cfg = suite_options();
    cfg.repeats = 1;
    cfg.min_batch_seconds = 0;
    cfg.max_batch_solves = 1;
    cfg.warmup = false;
    out = suite_folder(root, 'conditioning', cfg);
    rows = [];
    files = dir(fullfile(root, 'data', 'COND_*.mat'));
    assert(~isempty(files), 'Conditioning MAT files are missing.');
    for k = 1:numel(files)
        P = load(fullfile(files(k).folder, files(k).name));
        new = suite_measure(P, cfg);
        rows = [rows;new];
        suite_checkpoint(rows, out);
        clear P new
    end
end

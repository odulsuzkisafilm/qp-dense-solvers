function out = suite_folder(root, label, cfg)

    out = fullfile(root, 'results', [label '_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'))]);
    mkdir(out);
    save(fullfile(out, 'config.mat'), 'cfg');
    fid = fopen(fullfile(out, 'environment.txt'), 'w');
    fprintf(fid, 'MATLAB: %s\nPlatform: %s\nBLAS: %s\n', version, computer, version('-blas'));
    fclose(fid);
    fprintf('Results: %s\n', out);
end

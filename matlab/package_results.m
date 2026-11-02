function filename = package_results()

    root = fileparts(mfilename('fullpath'));
    filename = fullfile(root, ['results_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.zip']);
    zip(filename, {'results'}, root);
end

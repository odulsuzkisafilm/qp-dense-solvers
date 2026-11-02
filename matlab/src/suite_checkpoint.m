function suite_checkpoint(rows, out)

    T = struct2table(rows);
    writetable(T, fullfile(out, 'results.csv'));
    save(fullfile(out, 'results.mat'), 'T');
end

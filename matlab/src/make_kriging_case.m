function P = make_kriging_case(n, seed, root)

    A = load(fullfile(root, 'data', 'road_elevation.mat'));
    rng(seed, 'twister');
    ids = randperm(size(A.xy_km, 1), n+1);
    X = A.xy_km(ids(1:n), :);
    target = A.xy_km(ids(end), :);
    ell_km = 10;
    nugget = 1e-3;
    H = zeros(n);
    block = 256;
    for i = 1:block:n
        ii = i:min(i+block-1, n);
        for j = i:block:n
            jj = j:min(j+block-1, n);
            D = (X(ii, 1)-X(jj, 1)').^2+(X(ii, 2)-X(jj, 2)').^2;
            B = exp(-D/(2*ell_km^2));
            if i == j
                B = (B+B')/2;
            end
            H(ii, jj) = B;
            H(jj, ii) = B';
        end
    end
    H(1:n+1:end) = 1+nugget;
    g = exp(-sum((X-target).^2, 2)/(2*ell_km^2));
    P = struct('H', H, 'C', ones(1, n), 'g', g, 'd', 1, ...
        'x_ref', nan(n, 1), 'lambda_ref', NaN, 'f_constant', 0, ...
        'case_id', sprintf('KRIGING_n%d_seed%d', n, seed), 'family', 'KRIGING', ...
        'category', 'real_data_chosen_covariance_model', 'seed', seed, ...
        'ref_kind', 'no_independent_reference_large_case', 'kH_target', NaN, 'kC_target', 1, ...
        'kH_measured', NaN, 'kC_measured', 1, ...
        'observations', A.elevation_m(ids(1:n)), 'heldout_elevation', A.elevation_m(ids(end)), ...
        'source_rows', A.source_row(ids), 'ell_km', ell_km, 'nugget', nugget);
end

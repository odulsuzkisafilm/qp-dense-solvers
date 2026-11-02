function out = run_smoke_test()

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'src'));
    P = struct('H', diag([2 8]), 'C', [1 1], 'g', [0;0], 'd', 1);
    [x, l] = solve_qr_null(P.H, P.C, P.g, P.d);
    assert(norm(x-[.8;.2])<1e-12 && abs(l+1.6)<1e-12);
    [x, l] = solve_schur(P.H, P.C, P.g, P.d);
    assert(norm(x-[.8;.2])<1e-12 && abs(l+1.6)<1e-12);
    for name = {'DTOC3_N10.mat', 'AUG2DC_NX3_NY3.mat'}
        P = load(fullfile(root, 'data', name{1}));
        [x, l] = solve_qr_null(full(P.H), full(P.C), P.g, P.d);
        assert(norm(x-P.x_ref)/norm(P.x_ref)<1e-10);
        [x, l] = solve_schur(full(P.H), full(P.C), P.g, P.d);
        assert(norm(x-P.x_ref)/norm(P.x_ref)<1e-10);
    end
    out = run_cutest_tests('pilot');
end

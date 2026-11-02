function P = orthogonal_benchmark(P, seed)

    n = size(P.H, 1);
    rng(seed, 'twister');
    [U, ~] = qr(randn(n));
    h = full(diag(P.H));
    if nnz(P.H-spdiags(h, 0, n, n)) == 0
        Hnew = (U'.*h')*U;
    else
        Hnew = U'*(P.H*U);
    end
    Hnew = (Hnew+Hnew')/2;
    Cnew = P.C*U;
    gnew = U'*P.g;
    v = randn(n, 1);
    uv = U*v;
    qerr = abs(v'*Hnew*v-uv'*(P.H*uv))/max(1, abs(uv'*(P.H*uv)));
    cerr = norm(Cnew*v-P.C*uv)/max(1, norm(P.C*uv));
    assert(qerr<1e-10 && cerr<1e-10, 'Orthogonal transformation check failed');
    P.x_ref = U'*P.x_ref;
    P.H = Hnew;
    P.C = Cnew;
    P.g = gnew;
    clear U Hnew Cnew gnew
    P.case_id = sprintf('%s_orthogonal_seed%d', P.case_id, seed);
    P.family = [char(P.family) '_orthogonal'];
    P.seed = seed;
    P.category = 'equivalent_dense_representation';
    P.ref_kind = 'transformed_double_reference_not_certified';
    P.transform_quadratic_check = qerr;
    P.transform_constraint_check = cerr;
end

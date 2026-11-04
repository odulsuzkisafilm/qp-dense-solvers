function [inertia, pivot_min] = ldl_inertia(D)

    m = size(D, 1);
    dg = full(diag(D));
    sub = full(diag(D, -1));
    inertia = [0, 0, 0];
    pivot_min = Inf;
    i = 1;
    while i <= m
        if i < m && sub(i) ~= 0
            e = eig([dg(i), sub(i); sub(i), dg(i+1)]);
            i = i + 2;
        else
            e = dg(i);
            i = i + 1;
        end
        inertia = inertia + [sum(e > 0), sum(e < 0), sum(e == 0)];
        pivot_min = min(pivot_min, min(abs(e)));
    end
end

function qp_check_inputs(H, C, g, d)

    validateattributes(H, {'double'}, {'real', 'finite', '2d', 'nonempty', 'nonsparse'});
    validateattributes(C, {'double'}, {'real', 'finite', '2d', 'nonempty', 'nonsparse'});
    validateattributes(g, {'double'}, {'real', 'finite', 'column', 'nonempty', 'nonsparse'});
    validateattributes(d, {'double'}, {'real', 'finite', 'column', 'nonempty', 'nonsparse'});
    n = size(H, 1);
    p = size(C, 1);
    if size(H, 2) ~= n || size(C, 2) ~= n || numel(g) ~= n || numel(d) ~= p
        error('qp:Dimensions', 'Expected H(n,n), C(p,n), g(n,1), d(p,1).');
    end
    if p >= n
        error('qp:Dimensions', 'This project assumes 0 < p < n.');
    end
    if ~isequal(H, H')
        error('qp:Symmetry', 'Supply symmetric H; document any preprocessing.');
    end
end

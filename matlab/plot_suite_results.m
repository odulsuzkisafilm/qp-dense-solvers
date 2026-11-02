function plot_suite_results(out)

    T = readtable(fullfile(out, 'results.csv'), 'TextType', 'string');
    T.case_id = string(T.case_id);
    T.method = string(T.method);
    T.family = string(T.family);
    T.status = string(T.status);
    keys = unique(T.case_id, 'stable');
    summary = table();
    for k = 1:numel(keys)
        for method = ["QR", "Schur"]
            A = T(T.case_id == keys(k)&T.method == method, :);
            if isempty(A)
                continue;
            end
            ok = A.status == "ok"|A.status == "returned_warning";
            B = A(ok, :);
            med = NaN;
            lo = NaN;
            hi = NaN;
            xe = NaN;
            re = NaN;
            if ~isempty(B)
                med = median(B.seconds);
                lo = min(B.seconds);
                hi = max(B.seconds);
                xe = median(B.x_reference_error);
                re = median(B.kkt_fro_scaled);
            end
            R = table(keys(k), A.family(1), A.n(1), A.p(1), A.ratio(1), A.seed(1), method, ...
                med, lo, hi, xe, re, sum(ok), sum(A.status == "failed"|A.status == "allocation_failed"), sum(A.status == "memory_guard_skip"), ...
                'VariableNames', {'case_id', 'family', 'n', 'p', 'ratio', 'seed', 'method', 'median_seconds', ...
                'min_seconds', 'max_seconds', 'x_reference_error', 'kkt_fro_scaled', ...
                'successful_batches', 'failed_batches', 'skipped'});
            summary = [summary;R];
        end
    end
    writetable(summary, fullfile(out, 'summary.csv'));
    speed = table();
    for k = 1:numel(keys)
        A = summary(summary.case_id == keys(k), :);
        q = A.median_seconds(A.method == "QR");
        s = A.median_seconds(A.method == "Schur");
        ratio = NaN;
        if ~isempty(q) && ~isempty(s)
            ratio = q/s;
        end
        speed = [speed;table(keys(k), ratio, 'VariableNames', {'case_id', 'QR_seconds_over_Schur_seconds'})];
    end
    writetable(speed, fullfile(out, 'speed_ratios.csv'));
    families = unique(T.family, 'stable');
    slopes = table();
    for family = families'
        A = T(T.family == family, :);
        iscond = startsWith(family, "COND_");
        fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1100 450]);
        tiledlayout(1, 2);
        for panel = 1:2
            nexttile;
            hold on;
            grid on;
            if iscond
                xl = 'Measured cond_2(H)';
                if family == "COND_C"
                    xl = 'Measured cond_2(C)';
                end
                yl = 'Relative x error vs stored-input reference';
                if all(A.reference_kind == "manufactured_target_not_certified")
                    yl = 'Relative deviation from manufactured target (not certified error)';
                end
                if panel == 2
                    yl = 'KKT residual scaled by Frobenius norm';
                end
                set(gca, 'XScale', 'log', 'YScale', 'log');
            else
                xl = 'Number of variables n';
                if numel(unique(A.n)) == 1
                    xl = 'Constraint fraction p/n';
                end
                yl = 'Fresh solve time (seconds)';
                if panel == 2
                    yl = 'Seconds / n^3';
                end
                set(gca, 'YScale', 'log');
                if numel(unique(A.n))>1
                    set(gca, 'XScale', 'log');
                end
            end
            for method = ["QR", "Schur"]
                if iscond
                    B = A((A.status == "ok"|A.status == "returned_warning")&A.method == method, :);
                    if isempty(B)
                        continue;
                    end
                    xx = B.kH_measured;
                    if family == "COND_C"
                        xx = B.kC_measured;
                    end
                    if any(~isfinite(xx))
                        xx = B.kH_target;
                        xl = 'Prescribed cond_2(H)';
                        if family == "COND_C"
                            xx = B.kC_target;
                            xl = 'Prescribed cond_2(C)';
                        end
                    end
                    yy = B.x_reference_error;
                    if panel == 2
                        yy = B.kkt_fro_scaled;
                    end
                    plot(xx, yy, 'o', 'DisplayName', method);
                else
                    B = summary(summary.family == family & summary.method == method & isfinite(summary.median_seconds), :);
                    if isempty(B)
                        continue;
                    end
                    xx = B.n;
                    if numel(unique(A.n)) == 1
                        xx = B.ratio;
                    end
                    yy = B.median_seconds;
                    if panel == 2
                        yy = yy./B.n.^3;
                    end
                    ux = unique(xx);
                    med = zeros(size(ux));
                    low = med;
                    high = med;
                    for j = 1:numel(ux)
                        z = yy(xx == ux(j));
                        med(j) = median(z);
                        low(j) = min(z);
                        high(j) = max(z);
                    end
                    errorbar(ux, med, med-low, high-med, 'o-', 'DisplayName', method);
                    if panel == 1 && numel(unique(A.n))>1
                        for j = 2:numel(ux)
                            exponent = log(med(j)/med(j-1))/log(ux(j)/ux(j-1));
                            slopes = [slopes;table(family, method, ux(j-1), ux(j), exponent, ...
                                'VariableNames', {'family', 'method', 'n_low', 'n_high', 'local_time_exponent'})];
                        end
                    end
                end
            end
            xlabel(xl);
            ylabel(yl);
            legend('Location', 'best');
        end
        sgtitle(sprintf('%s: %d failed attempts, %d memory skips', family, ...
            sum(A.status == "failed"|A.status == "allocation_failed"), sum(A.status == "memory_guard_skip")), 'Interpreter', 'none');
        exportgraphics(fig, fullfile(out, char(family+".png")), 'Resolution', 180);
        savefig(fig, fullfile(out, char(family+".fig")));
        close(fig);
    end
    if ~isempty(slopes)
        writetable(slopes, fullfile(out, 'scaling_slopes.csv'));
    end
    fprintf('CSV summaries and plots saved in %s\n', out);
end

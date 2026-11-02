function cfg = final_options(timing)

    cfg = suite_options();
    cfg.memory_budget_GiB = 10;
    cfg.peak_dense_arrays = 8;
    cfg.repeats = 7;
    cfg.min_batch_seconds = .2;
    cfg.max_batch_solves = 100;
    if ~timing
        cfg.repeats = 1;
        cfg.min_batch_seconds = 0;
        cfg.max_batch_solves = 1;
        cfg.warmup = false;
    end
end

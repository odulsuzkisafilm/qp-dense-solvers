function cfg = suite_options()

    cfg.memory_budget_GiB = 10;
    cfg.peak_dense_arrays = 12;
    cfg.repeats = 5;
    cfg.min_batch_seconds = 0.20;
    cfg.max_batch_solves = 100;
    cfg.warmup = true;
end

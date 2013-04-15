using Stats

@assert norm(distances([1.0 0.0; 0.0 1.0]) - [0.0 sqrt(2); sqrt(2) 0.0]) < 1e-4
@assert norm(distances([1.0 0.0; 0.0 1.0], [0.0 1.0; 1.0 0.0]) - [sqrt(2) 0.0; 0.0 sqrt(2)]) < 1e-4

@assert invlogit(-800) > 0.0
@assert invlogit(800) <= 1.0

@assert isinf(logit(0.0))
@assert norm(invlogit(logit(0.01)) - 0.01) < 1e-4
@assert norm(invlogit(logit(0.50)) - 0.50) < 1e-4
@assert norm(invlogit(logit(0.99)) - 0.99) < 1e-4
@assert isinf(logit(1.0))

@assert isfinite(logsumexp([1.0, 1.0, 1.0]))
@assert isfinite(logsumexp([1.0, 1.0, 1.0e7]))
@assert isfinite(logsumexp([-1.0, -1.0, -7.5]))
@assert isfinite(logsumexp([-10.0, -10.0, -75000.0]))

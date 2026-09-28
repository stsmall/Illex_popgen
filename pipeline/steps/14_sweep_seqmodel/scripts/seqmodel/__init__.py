STATES = ("N", "LL", "C", "LR")
STATE_IDX = {s: i for i, s in enumerate(STATES)}
EPS = 1e-12
TRAINING_PRIORS = (0.2, 0.4, 0.4)  # (N, L, C)

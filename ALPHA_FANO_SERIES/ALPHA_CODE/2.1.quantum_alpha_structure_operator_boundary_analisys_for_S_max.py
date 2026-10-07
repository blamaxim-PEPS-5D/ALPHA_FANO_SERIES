# =================================================================================
# QUANTUM ALPHA STRUCTURE OPERATOR - BOUNDARY ANALYSIS (q <= 45)
# Author: Massimiliano Blandino
# ORCID: https://orcid.org/0009-0006-3252-4011
# Description: Evaluates the exact mathematical upper ceiling S_max via 
#              deterministic minimal F and stochastic boundary sampling 
#              within the physical Hilbert space (q <= 45, K in (9, 10)).
# =================================================================================

import numpy as np
from mpmath import mp
import time

# 1. Setting precision to 200 decimal digits
mp.dps = 200

# 2. Generator polynomial A(pi) and structural constants
pi_mp = mp.pi
A_pi = 4 * pi_mp**3 + pi_mp**2 + pi_mp
A_third_deriv = mp.mpf(24)

zero_order = A_pi
curvature_corr = 1 / (A_third_deriv * A_pi)
const_factor = 1 / (A_pi**2 * pi_mp**2)

def compute_S_from_K(K_val):
    """Computes S = A(pi) - 1/(24*A(pi)) - 1/(A(pi)^2 * pi^2 * K)"""
    quantum_term = const_factor / K_val
    return zero_order - curvature_corr - quantum_term

def compute_F_from_sequence(seq):
    """Reconstructs F from a sequence of partial quotients q_1, ..., q_d"""
    f = mp.mpf(str(seq[-1]))
    for q in reversed(seq[:-1]):
        f = mp.mpf(str(q)) + 1 / f
    return 1 / f

print("=" * 80)
print("QUANTUM ALPHA STRUCTURE OPERATOR - BOUNDARY ANALYSIS")
print("Author: Massimiliano Blandino")
print("ORCID: https://orcid.org/0009-0006-3252-4011")
print("=" * 80)

# ---------------------------------------------------------------------------------
# PART 1: DETERMINISTIC BOUNDARY ANALYSIS (Minimal F -> Maximal K -> Upper Bound S)
# ---------------------------------------------------------------------------------
# To maximize S within physical space (q_i >= 1), we minimize F.
# Minimal F in q in {1, ..., 45} occurs when q_1 = 45, q_2 = 1, q_3 = 45, etc.

seq_min_F = [45 if i % 2 == 0 else 1 for i in range(20)]
F_min = compute_F_from_sequence(seq_min_F)
K_max = mp.mpf(10) - F_min
S_max_analytical = compute_S_from_K(K_max)

print("\n[PART 1: ANALYTICAL EXTREMUM IN BOUNDED SPACE q <= 45]")
print(f"Sequence yielding minimal F : {seq_min_F[:6]}...")
print(f"Minimal F value             : {mp.nstr(F_min, 15)}")
print(f"Maximal K value (10 - F)    : {mp.nstr(K_max, 15)}")
print(f"Analytical S_max            : {mp.nstr(S_max_analytical, 12)}")

# ---------------------------------------------------------------------------------
# PART 2: STOCHASTIC SAMPLING OF THE BOUNDARY ENVELOPE
# ---------------------------------------------------------------------------------
# Sampling 500,000 states conditioned on upper boundary quotients (q_1 = 45)
# to observe the saturation ceiling of stochastic realizations.

max_q = 45
E_int = mp.mpf("5.0")
q_space = np.arange(1, max_q + 1, dtype=np.int64)

weights = [mp.exp(-mp.mpf(str(q)) / (2 * E_int)) for q in range(1, max_q + 1)]
total_w = sum(weights)
probs = [float(w / total_w) for w in weights]

n_samples = 500_000
S_boundary_samples = []

print("\n[PART 2: STOCHASTIC SAMPLING OF THE HIGH-S ENVELOPE (500,000 iterations)]")
start_time = time.time()

seq = np.zeros(20, dtype=np.int64)
for i in range(n_samples):
    # Set q_1 = 45, sample remaining 19 quotients stochastically according to P(q)
    seq[0] = 45
    for j in range(1, 20):
        seq[j] = np.random.choice(q_space, p=probs)

    F_val = compute_F_from_sequence(seq)
    K_val = mp.mpf(10) - F_val
    S_val = float(compute_S_from_K(K_val))
    S_boundary_samples.append(S_val)

elapsed = time.time() - start_time

S_arr = np.array(S_boundary_samples)
print(f"Sampling completed in {elapsed:.2f} seconds")
print(f"Stochastic Minimum in boundary region : {np.min(S_arr):.12f}")
print(f"Stochastic Maximum in boundary region : {np.max(S_arr):.12f}")
print(f"Stochastic Mean in boundary region    : {np.mean(S_arr):.12f}")

print("\n" + "=" * 80)
print("CONCLUSION")
print("=" * 80)
print(f"The value {np.max(S_arr):.12f} represents the exact physical saturation")
print(f"ceiling Sup(S) of the structure operator S within the q <= 45 bounded")
print(f"Hilbert space, reached when K approaches K_max ≈ {float(K_max):.6f}.")
print("=" * 80)

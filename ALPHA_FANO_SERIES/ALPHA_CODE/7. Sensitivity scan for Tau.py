# -*- coding: utf-8 -*-
"""
Sensitivity scan for tau (E_int parameter) in the quantum alpha model.
Tests invariance of <S> for tau = 3.0, 4.0, 5.0, 6.0, 7.0.
"""

import numpy as np
from mpmath import mp
import time

# Precision
mp.dps = 200

# Constants
pi_mp = mp.pi

def A(x):
    return 4*x**3 + x**2 + x

A_pi = A(pi_mp)
A_third_derivative = mp.mpf(24)
zero_order_term = A_pi
curvature_correction = 1 / (A_third_derivative * A_pi)
const_factor = 1 / (A_pi**2 * pi_mp**2)

codata_2022 = mp.mpf("137.035999177")

# Parameters
depth = 20
max_q = 45
n_iterations = 100_000  # 100k per tau value (fast enough)

def compute_structure_value(sequence):
    f = mp.mpf(str(sequence[-1]))
    for q in reversed(sequence[:-1]):
        f = mp.mpf(str(q)) + 1 / f
    f = 1 / f
    K = mp.mpf(10) - f
    quantum_term = const_factor / K
    S = zero_order_term - curvature_correction - quantum_term
    return float(S)

# Scan over tau values
tau_values = [3.0, 4.0, 5.0, 6.0, 7.0]
q_space = np.arange(1, max_q + 1)

print("=" * 80)
print("SENSITIVITY SCAN FOR tau (E_int)")
print("=" * 80)
print(f"Depth: {depth}")
print(f"Iterations per tau: {n_iterations:,}")
print(f"q-space: 1 to {max_q}")
print("-" * 80)
print(f"{'tau':<8} {'Mean <S>':<20} {'Difference from CODATA':<25} {'Std Dev':<12}")
print("-" * 80)

results = {}

for tau in tau_values:
    # Build probability distribution for this tau
    E_int = mp.mpf(tau)
    weights = [mp.exp(-q / (2 * E_int)) for q in range(1, max_q + 1)]
    total_weight = sum(weights)
    probabilities = [float(w / total_weight) for w in weights]
    
    S_values = []
    sequence = np.zeros(depth, dtype=np.int64)
    
    start_time = time.time()
    
    for i in range(n_iterations):
        for j in range(depth):
            sequence[j] = np.random.choice(q_space, p=probabilities)
        S = compute_structure_value(sequence)
        S_values.append(S)
    
    elapsed = time.time() - start_time
    
    mean_S = np.mean(S_values)
    std_S = np.std(S_values)
    diff = mean_S - float(codata_2022)
    
    results[tau] = {'mean': mean_S, 'std': std_S, 'diff': diff, 'time': elapsed}
    
    print(f"{tau:<8.1f} {mean_S:<20.12f} {diff:<+25.2e} {std_S:<12.2e}")

print("-" * 80)
print("\nINVARIANCE CHECK:")
print("If |diff| varies by less than 1e-8 across tau, invariance is confirmed.")

# Check variation
means = [results[tau]['mean'] for tau in tau_values]
max_diff_means = max(means) - min(means)
print(f"\nMax variation in <S> across tau values: {max_diff_means:.2e}")
if max_diff_means < 1e-8:
    print("✅ INVARIANCE CONFIRMED: <S> stable within 1e-8")
else:
    print(f"⚠️ Variation detected: {max_diff_means:.2e}")

print("\n" + "=" * 80)
print("SCAN COMPLETE")
print("=" * 80)
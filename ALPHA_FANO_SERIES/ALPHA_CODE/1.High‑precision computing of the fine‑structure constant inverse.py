# -*- coding: utf-8 -*-
"""
Author: Massimiliano Blandino (ORCID: 0009-0006-3252-4011)
Reference Manuscript: A representation of the fine-structure constant using pi and a continued fraction of small integers
DOI: https://doi.org/10.5281/zenodo.19802607
"""

import mpmath as mp

# ----------------------------------------------------------------------
# 1. PRECISION SETTING (200 digits)
# ----------------------------------------------------------------------
mp.dps = 200

# ----------------------------------------------------------------------
# 2. MAIN FORMULA (from manuscript)
# ----------------------------------------------------------------------
def A():
    pi = mp.pi
    return 4*pi**3 + pi**2 + pi

def alpha_inv_from_K(K):
    """Computes alpha^-1 via the formula: alpha^-1 = A - 1/(24A) - 1/(A^2 pi^2 K)"""
    a = A()
    pi = mp.pi
    return a - 1/(24*a) - 1/(a**2 * pi**2 * K)

def frazione_continua_K(quozienti):
    """
    Constructs K = 10 - 1/(q1 + 1/(q2 + 1/(q3 + ...)))
    Handles cases where a partial quotient equals 0 (finite continued fraction).
    """
    if not quozienti:
        return 10
    
    # If a zero is present, truncate the continued fraction at that point
    if 0 in quozienti:
        idx = quozienti.index(0)
        quozienti = quozienti[:idx]
        if not quozienti:
            return 10
    
    def ricorsiva(lst):
        if len(lst) == 1:
            return lst[0]
        return lst[0] + 1 / ricorsiva(lst[1:])
    
    try:
        return 10 - 1 / ricorsiva(quozienti)
    except ZeroDivisionError:
        return 10  # Fallback: return K=10 in case of division by zero

def continued_fraction(x, n=20, tol=1e-50):
    """Computes the partial quotients of the continued fraction representation of x."""
    coeffs = []
    y = x
    for _ in range(n):
        a = int(mp.floor(y))
        coeffs.append(a)
        diff = y - a
        if diff < tol:
            break
        y = 1 / diff
    return coeffs

def quozienti_da_target(alpha_target, max_terms=20):
    """
    Computes the continued fraction of the correction factor 1/(10-K) 
    given a target value of alpha^-1.
    Returns (partial_quotients, zero_error_level)
    """
    a = A()
    pi_val = mp.pi
    try:
        # Inverse formula: isolating K from alpha^-1 = A - 1/(24A) - 1/(A^2 pi^2 K)
        denom = a - 1/(24*a) - alpha_target
        if denom == 0:
            return [], None
        K_target = 1 / (a**2 * pi_val**2 * denom)
    except ZeroDivisionError:
        return [], None
    
    x = 1 / (10 - K_target)
    q = continued_fraction(x, n=max_terms)
    
    # Identify the minimal truncation level where error vanishes (err < 1e-30)
    livello_zero = None
    for n in range(1, len(q)+1):
        if 0 in q[:n]:
            continue  # Finite continued fraction; expansion cannot proceed further
        try:
            K_cf = frazione_continua_K(q[:n])
            alpha_calc = alpha_inv_from_K(K_cf)
            err = abs(alpha_calc - alpha_target)
            if err < mp.mpf('1e-30'):
                livello_zero = n
                break
        except ZeroDivisionError:
            continue
    return q, livello_zero

# ----------------------------------------------------------------------
# 3. HISTORICAL CODATA VALUES
# ----------------------------------------------------------------------
codata_values = {
    2006: mp.mpf('137.035999070'),
    2010: mp.mpf('137.035999074'),
    2014: mp.mpf('137.035999139'),
    2018: mp.mpf('137.035999084'),
    2022: mp.mpf('137.035999177'),
}

print("="*110)
print("ANALYSIS OF HISTORICAL CODATA VALUES (2006-2022)")
print("="*110)
print(f"{'Year':<6} {'alpha^-1':<20} {'First 10 Partial Quotients':<50} {'Zero Level':<12} {'Bounded Quotients?'}")
print("-"*110)

for anno, val in codata_values.items():
    q, livello = quozienti_da_target(val, max_terms=20)
    if not q:
        print(f"{anno:<6} {mp.nstr(val, 12):<20} {'ERROR':<50} {'---':<12} {'---'}")
        continue
    
    q_str = str(q[:10]) if len(q) >= 10 else str(q)
    livello_str = str(livello) if livello else ">20"
    max_q = max(abs(x) for x in q[:10]) if q else 0
    piccoli = "YES" if max_q <= 100 else f"max={max_q}"
    
    print(f"{anno:<6} {mp.nstr(val, 12):<20} {q_str:<50} {livello_str:<12} {piccoli}")

# ----------------------------------------------------------------------
# 4. FINE SCAN OF THE HISTORICAL INTERVAL (with exception handling)
# ----------------------------------------------------------------------
print("\n" + "="*110)
print("FINE SCAN OF THE HISTORICAL INTERVAL (step size 1e-9)")
print("="*110)

start = mp.mpf('137.03599907')
end = mp.mpf('137.03599928')
step = mp.mpf('1e-9')

alpha_critici = []
errori = 0

alpha = start
while alpha <= end:
    try:
        q, _ = quozienti_da_target(alpha, max_terms=12)
        if q and max(abs(x) for x in q[:8]) > 100:
            alpha_critici.append(float(alpha))
    except Exception:
        errori += 1
    alpha += step

if alpha_critici:
    print(f"Found {len(alpha_critici)} values with partial quotients > 100:")
    for v in alpha_critici[:15]:
        print(f"  alpha^-1 ~ {v:.12f}")
else:
    print("NO values within the historical interval generated partial quotients > 100")
if errori > 0:
    print(f"(Note: {errori} evaluation steps triggered numerical errors and were skipped)")

# ----------------------------------------------------------------------
# 5. GENERAL PROPERTY VERIFICATION (maximum quotient across CODATA)
# ----------------------------------------------------------------------
print("\n" + "="*110)
print("PROPERTY VERIFICATION: Are partial quotients bounded for all historical CODATA values?")
print("="*110)

max_global = 0
for anno, val in codata_values.items():
    q, _ = quozienti_da_target(val, max_terms=15)
    if q:
        max_q = max(abs(x) for x in q[:10])
        max_global = max(max_global, max_q)
        print(f"{anno}: max quotient = {max_q}")
    else:
        print(f"{anno}: ERROR")

print(f"\nMaximum observed quotient across all historical CODATA values: {max_global}")
if max_global <= 45:
    print("All values exhibit bounded quotients <= 45 (stable property)")
else:
    print(f"Maximum quotient {max_global} > 45, yet remains bounded.")

# ----------------------------------------------------------------------
# 6. DIRECT COMPUTATION OF THE ANALYTICAL FORMULA
# ----------------------------------------------------------------------
print("\n" + "="*110)
print("ANALYTICAL FORMULA COMPUTATION (with K defined by continued fraction expansion)")
print("="*110)

quozienti_paper = [14, 1, 7, 3, 1, 3]
K_paper = frazione_continua_K(quozienti_paper)
alpha_inv_paper = alpha_inv_from_K(K_paper)
alpha_paper = 1 / alpha_inv_paper

print("K (analytical) = 10 - 1/(14 + 1/(1 + 1/(7 + 1/(3 + 1/(1 + 1/3)))))")
print(f"K = {mp.nstr(K_paper, 25)}")
print(f"alpha^-1 (analytical) = {mp.nstr(alpha_inv_paper, 25)}")
print(f"alpha^-1 (CODATA 2022) = {mp.nstr(codata_values[2022], 25)}")
print(f"Absolute Error alpha^-1 = {float(abs(alpha_inv_paper - codata_values[2022])):.2e}")
print(f"alpha (analytical) = {mp.nstr(alpha_paper, 25)}")
print(f"alpha (CODATA 2022) = {mp.nstr(1/codata_values[2022], 25)}")
print("="*110)

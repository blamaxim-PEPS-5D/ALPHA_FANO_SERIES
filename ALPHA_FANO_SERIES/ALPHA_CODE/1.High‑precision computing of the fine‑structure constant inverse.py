# -*- coding: utf-8 -*-
import mpmath as mp

# ----------------------------------------------------------------------
# 1. PRECISIONE (200 cifre)
# ----------------------------------------------------------------------
mp.dps = 200

# ----------------------------------------------------------------------
# 2. FORMULA PRINCIPALE (dal paper)
# ----------------------------------------------------------------------
def A():
    pi = mp.pi
    return 4*pi**3 + pi**2 + pi

def alpha_inv_from_K(K):
    """Calcola alpha^-1 dalla formula: alpha^-1 = A - 1/(24A) - 1/(A^2 pi^2 K)"""
    a = A()
    pi = mp.pi
    return a - 1/(24*a) - 1/(a**2 * pi**2 * K)

def frazione_continua_K(quozienti):
    """
    Costruisce K = 10 - 1/(q1 + 1/(q2 + 1/(q3 + ...)))
    Gestisce il caso in cui un quoziente e' 0 (frazione continua finita).
    """
    if not quozienti:
        return 10
    
    # Se c'e' uno zero, tronchiamo la frazione continua a quel punto
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
        return 10  # fallback: se c'e' ancora errore, restituisci K=10

def continued_fraction(x, n=20, tol=1e-50):
    """Calcola i quozienti della frazione continua di x."""
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
    Calcola la frazione continua della correzione 1/(10-K) a partire da un valore target di alpha^-1.
    Restituisce (quozienti, livello_zero)
    """
    a = A()
    pi_val = mp.pi
    try:
        # Formula inversa: isolo K da alpha^-1 = A - 1/(24A) - 1/(A^2 pi^2 K)
        denom = a - 1/(24*a) - alpha_target
        if denom == 0:
            return [], None
        K_target = 1 / (a**2 * pi_val**2 * denom)
    except ZeroDivisionError:
        return [], None
    
    x = 1 / (10 - K_target)
    q = continued_fraction(x, n=max_terms)
    
    # Trova il livello minimo in cui l'errore si annulla (err < 1e-30)
    livello_zero = None
    for n in range(1, len(q)+1):
        if 0 in q[:n]:
            continue  # frazione continua finita, non possiamo andare oltre
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
# 3. VALORI CODATA STORICI
# ----------------------------------------------------------------------
codata_values = {
    2006: mp.mpf('137.035999070'),
    2010: mp.mpf('137.035999074'),
    2014: mp.mpf('137.035999139'),
    2018: mp.mpf('137.035999084'),
    2022: mp.mpf('137.035999177'),
}

print("="*110)
print("ANALISI DEI VALORI CODATA STORICI (2006-2022)")
print("="*110)
print(f"{'Anno':<6} {'alpha^-1':<20} {'Primi 10 quozienti':<50} {'Livello zero':<12} {'Quozienti piccoli?'}")
print("-"*110)

for anno, val in codata_values.items():
    q, livello = quozienti_da_target(val, max_terms=20)
    if not q:
        print(f"{anno:<6} {mp.nstr(val, 12):<20} {'ERRORE':<50} {'---':<12} {'---'}")
        continue
    
    q_str = str(q[:10]) if len(q) >= 10 else str(q)
    livello_str = str(livello) if livello else ">20"
    max_q = max(abs(x) for x in q[:10]) if q else 0
    piccoli = "YES" if max_q <= 100 else f"max={max_q}"
    
    print(f"{anno:<6} {mp.nstr(val, 12):<20} {q_str:<50} {livello_str:<12} {piccoli}")

# ----------------------------------------------------------------------
# 4. SCANSIONE FINE DELL'INTERVALLO STORICO (con gestione errori)
# ----------------------------------------------------------------------
print("\n" + "="*110)
print("SCANSIONE FINE DELL'INTERVALLO STORICO (step 1e-9)")
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
    print(f"Trovati {len(alpha_critici)} valori con quozienti > 100:")
    for v in alpha_critici[:15]:
        print(f"  alpha^-1 ~ {v:.12f}")
else:
    print("NESSUN valore nell'intervallo storico ha prodotto quozienti > 100")
if errori > 0:
    print(f"(Nota: {errori} valori hanno causato errori e sono stati saltati)")

# ----------------------------------------------------------------------
# 5. VERIFICA PROPRIETA GENERALE (massimo quoziente tra i CODATA)
# ----------------------------------------------------------------------
print("\n" + "="*110)
print("VERIFICA PROPRIETA: tutti i valori CODATA storici hanno quozienti limitati?")
print("="*110)

max_global = 0
for anno, val in codata_values.items():
    q, _ = quozienti_da_target(val, max_terms=15)
    if q:
        max_q = max(abs(x) for x in q[:10])
        max_global = max(max_global, max_q)
        print(f"{anno}: max quoziente = {max_q}")
    else:
        print(f"{anno}: ERRORE")

print(f"\nMassimo quoziente osservato tra tutti i valori CODATA storici: {max_global}")
if max_global <= 45:
    print("Tutti i valori hanno quozienti <= 45 (proprieta stabile)")
else:
    print(f"Valore massimo {max_global} > 45, ma comunque piccolo.")

# ----------------------------------------------------------------------
# 6. CALCOLO DIRETTO DELLA FORMULA DEL PAPER
# ----------------------------------------------------------------------
print("\n" + "="*110)
print("CALCOLO DELLA FORMULA ESATTA DEL PAPER (con K definito dalla frazione continua)")
print("="*110)

quozienti_paper = [14, 1, 7, 3, 1, 3]
K_paper = frazione_continua_K(quozienti_paper)
alpha_inv_paper = alpha_inv_from_K(K_paper)
alpha_paper = 1 / alpha_inv_paper

print("K (dal paper) = 10 - 1/(14 + 1/(1 + 1/(7 + 1/(3 + 1/(1 + 1/3)))))")
print(f"K = {mp.nstr(K_paper, 25)}")
print(f"alpha^-1 (formula) = {mp.nstr(alpha_inv_paper, 25)}")
print(f"alpha^-1 (CODATA 2022) = {mp.nstr(codata_values[2022], 25)}")
print(f"Errore alpha^-1    = {float(abs(alpha_inv_paper - codata_values[2022])):.2e}")
print(f"alpha (formula) = {mp.nstr(alpha_paper, 25)}")
print(f"alpha (CODATA 2022) = {mp.nstr(1/codata_values[2022], 25)}")
print("="*110)
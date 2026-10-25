_CURRENT_WEIGHTS = None
_GB_CACHE = {}
STOP_AFTER_FILTER = True

# Auto-generated from LND_colab.ipynb
from sage.all import *
import time, sys, traceback
print(f"SAGE {version()}", flush=True)

print("### CELL 0 START ###", flush=True)
from itertools import product as iproduct
print("### CELL 0 END ###", flush=True)

print("### CELL 1 START ###", flush=True)
from sage.arith.misc import gcd
print("### CELL 1 END ###", flush=True)

print("### CELL 2 START ###", flush=True)
def monomials_of_weighted_degree(R, weights, target_deg, max_deg_per_var=30):
    n = R.ngens()
    def gen(idx, remaining, exp):
        if idx == n - 1:
            w = weights[idx]
            if w > 0:
                if remaining % w == 0:
                    e = remaining // w
                    if e >= 0:
                        yield R.monomial(*(exp + [e]))
            elif w == 0:
                if remaining == 0:
                    yield R.monomial(*(exp + [0]))
            return
        w = weights[idx]
        if w > 0:
            max_e = min(remaining // w, max_deg_per_var)
        else:
            max_e = 0 if remaining != 0 else max_deg_per_var
        for e in range(max_e + 1):
            yield from gen(idx + 1, remaining - e * w, exp + [e])
    return list(gen(0, target_deg, []))
print("### CELL 2 END ###", flush=True)

print("### CELL 3 START ###", flush=True)
# ------------------------------------------------------------------
# 1. Multi-graded monomial generator (FIXED)
# ------------------------------------------------------------------
def monomials_of_multigraded_degree(R, weight_vecs, target_deg_vec, max_deg_per_var=30):
    n = R.ngens()
    q = len(target_deg_vec)

    def gen(idx, remaining_vec, exp):
        if idx == n - 1:
            w = weight_vecs[idx]
            e = None
            valid = True
            for k in range(q):
                if w[k] == 0:
                    if remaining_vec[k] != 0:
                        valid = False
                        break
                else:
                    if remaining_vec[k] % w[k] != 0:
                        valid = False
                        break
                    candidate = remaining_vec[k] // w[k]
                    if e is None:
                        e = candidate
                    elif e != candidate:
                        valid = False
                        break
            if valid and e is not None and e >= 0:
                yield R.monomial(*(exp + [e]))
            return

        w = weight_vecs[idx]
        max_e = max_deg_per_var
        for k in range(q):
            if w[k] > 0:
                max_e = min(max_e, remaining_vec[k] // w[k])
        if max_e < 0:
            return
        for e in range(max_e + 1):
            new_remaining = tuple(remaining_vec[k] - e * w[k] for k in range(q))
            if all(r >= 0 for r in new_remaining):
                yield from gen(idx + 1, new_remaining, exp + [e])

    return list(gen(0, tuple(target_deg_vec), []))
print("### CELL 3 END ###", flush=True)

print("### CELL 4 START ###", flush=True)
class LND:
    def __init__(self, images, ring=None):
        self.images = images
        if ring is None:
            self.R = images[0].parent()
        else:
            self.R = ring
        self.n = len(images)
        self.vars = self.R.gens()
    def __call__(self, f):
        if f == 0:
            return self.R(0)
        result = self.R(0)
        for i, var in enumerate(self.vars):
            deriv = f.derivative(var)
            if deriv != 0:
                result += deriv * self.images[i]
        return result
    def iterate(self, f, k):
        result = f
        for _ in range(k):
            result = self(result)
            if result == 0:
                break
        return result
    def nilpotency_index(self, f, max_iter=50):
        cur = f
        for k in range(max_iter + 1):
            if cur == 0:
                return k
            cur = self(cur)
        return max_iter + 1
    def is_in_kernel(self, f):
        return self(f) == 0
    def __repr__(self):
        return f"LND({', '.join(str(img) for img in self.images)})"
print("### CELL 4 END ###", flush=True)

print("### CELL 5 START ###", flush=True)
def maubach_kernel(D, weights, D_degree, max_degree):
    R = D.R
    kernel_by_degree = {}
    base_ring = QQ
    for d in range(max_degree + 1):
        monomials_d = monomials_of_weighted_degree(R, weights, d)
        if not monomials_d:
            kernel_by_degree[d] = []
            continue
        target_deg = d + D_degree
        if target_deg < 0:
            kernel_by_degree[d] = list(monomials_d)
            continue
        monomials_t = monomials_of_weighted_degree(R, weights, target_deg)
        if not monomials_t:
            kernel_by_degree[d] = list(monomials_d)
            continue
        rows = []
        for m in monomials_d:
            Dm = D(m)
            row = [Dm.monomial_coefficient(tm) for tm in monomials_t]
            row = [base_ring(c) for c in row]
            rows.append(row)
        M = matrix(base_ring, rows)
        K = M.left_kernel()
        kernel_polys = []
        for vec in K.basis():
            poly = sum(c * m for c, m in zip(vec, monomials_d))
            if poly != 0:
                kernel_polys.append(poly)
        kernel_by_degree[d] = kernel_polys
    return kernel_by_degree
print("### CELL 5 END ###", flush=True)

print("### CELL 6 START ###", flush=True)
# ------------------------------------------------------------------
# 2. Multi-graded Maubach algorithm
# ------------------------------------------------------------------
def maubach_kernel_multigraded(D, weight_vecs, D_degree_vec, max_degree_vec):
    R = D.R
    q = len(D_degree_vec)
    base_ring = QQ

    ranges = [range(int(m) + 1) for m in max_degree_vec]
    kernel_by_degree = {}

    for deg_vec in iproduct(*ranges):
        monomials_d = monomials_of_multigraded_degree(R, weight_vecs, deg_vec)
        if not monomials_d:
            continue

        target_deg_vec = tuple(deg_vec[k] + D_degree_vec[k] for k in range(q))

        if any(t < 0 for t in target_deg_vec):
            kernel_by_degree[deg_vec] = list(monomials_d)
            continue

        monomials_t = monomials_of_multigraded_degree(R, weight_vecs, target_deg_vec)
        if not monomials_t:
            kernel_by_degree[deg_vec] = list(monomials_d)
            continue

        rows = []
        for m in monomials_d:
            Dm = D(m)
            row = [Dm.monomial_coefficient(tm) for tm in monomials_t]
            row = [base_ring(c) for c in row]
            rows.append(row)

        M = matrix(base_ring, rows)
        K = M.left_kernel()
        kernel_polys = []
        for vec in K.basis():
            poly = sum(c * m for c, m in zip(vec, monomials_d))
            if poly != 0:
                kernel_polys.append(poly)
        if kernel_polys:
            kernel_by_degree[deg_vec] = kernel_polys

    return kernel_by_degree
print("### CELL 6 END ###", flush=True)

print("### CELL 7 START ###", flush=True)
def filter_new_generators(kernel_by_degree, R):
    """
    Keep only new generators (not products of previous ones).
    """
    generators = []
    new_by_degree = {}
    for d in sorted(kernel_by_degree.keys()):
        new_at_d = []
        for p in kernel_by_degree[d]:
            if p == 1:
                continue
            is_new = True
            for (gd, g) in generators:
                if g != 0 and p in R.ideal(g):
                    is_new = False
                    break
            if is_new:
                new_at_d.append(p)
                generators.append((d, p))
        if new_at_d:
            new_by_degree[d] = new_at_d
    return new_by_degree
print("### CELL 7 END ###", flush=True)

print("### CELL 8 START ###", flush=True)
# ------------------------------------------------------------------
# 3. Filter: keep only NEW generators
# ------------------------------------------------------------------
def filter_new_generators_multigraded(kernel_by_degree, R):
    generators = []
    new_by_degree = {}
    sorted_degrees = sorted(kernel_by_degree.keys(), key=lambda d: (sum(d), d))

    for d in sorted_degrees:
        new_at_d = []
        for p in kernel_by_degree[d]:
            if p == 1:
                continue
            is_new = True
            for (gd, g) in generators:
                if g != 0 and p in R.ideal(g):
                    is_new = False
                    break
            if is_new:
                new_at_d.append(p)
                generators.append((d, p))
        if new_at_d:
            new_by_degree[d] = new_at_d
    return new_by_degree
print("### CELL 8 END ###", flush=True)

print("### CELL 9 START ###", flush=True)
def filter_new_generators_multigraded(kernel_by_degree, R, main_vars):
    """
    main_vars: list of variables (e.g., [X,Y,Z]) that must appear in a generator.
    Filters out:
      - Base ring elements (no main variable appears)
      - Polynomials in the ideal generated by all previous generators
    """
    generators = []
    new_by_degree = {}
    sorted_degrees = sorted(kernel_by_degree.keys(), key=lambda d: (sum(d), d))

    for d in sorted_degrees:
        new_at_d = []
        for p in kernel_by_degree[d]:
            if p == 1:
                continue
            # Skip base ring elements (no main variable involved)
            if not any(p.degree(v) > 0 for v in main_vars):
                continue
            # Check if p is in the ideal generated by ALL previous generators
            is_new = True
            if generators:
                I = R.ideal([g for (_, g) in generators])
                if p in I:
                    is_new = False
            if is_new:
                new_at_d.append(p)
                generators.append((d, p))
        if new_at_d:
            new_by_degree[d] = new_at_d
    return new_by_degree
print("### CELL 9 END ###", flush=True)

print("### CELL 10 START ###", flush=True)
# ----------------------------------------------------------------------
# 4. Essen's algorithm (corrected coefficient extraction)
# ----------------------------------------------------------------------
def essen_generators(D):
    R = D.R
    X, Y, Z = D.vars
    f, g, h = D.images[0], D.images[1], D.images[2]
    F = g*X - f*Y
    K = R.fraction_field()
    T = PolynomialRing(K, 'T').gen()
    h_T = K(h).subs({X: T/g, Y: 0})
    poly_T = h_T.polynomial(T)
    coeffs = poly_T.coefficients(sparse=False)
    coeffs_R = []
    for c in coeffs:
        if c.denominator() == 1:
            coeffs_R.append(c.numerator())
        else:
            coeffs_R.append(c)
    all_vals = [g] + coeffs_R
    d = gcd(all_vals)
    g_prime = g / d
    a_prime = [c / d for c in coeffs_R]
    G = f*Z - h*X
    H_num = g_prime * G + sum(a_prime[i] * F**(i+1) for i in range(len(a_prime)))
    H = H_num / f
    return [F, G, H]
print("### CELL 10 END ###", flush=True)

print("### CELL 11 START ###", flush=True)
# ================================================================
# TEST 1: D(X)=0, D(Y)=X, D(Z)=Y
# ================================================================
print("=== TEST 1: D(X)=0, D(Y)=X, D(Z)=Y ===")
R.<X,Y,Z> = PolynomialRing(QQ, order='degrevlex')
D = LND([R(0), X, Y])
weights = (1,1,1)
D_degree = 0
max_deg = 5

kernel = maubach_kernel(D, weights, D_degree, max_deg)
new_gens = filter_new_generators(kernel, R)

print("New generators by degree (minimal generating set):")
for d in range(max_deg + 1):
    polys = new_gens.get(d, [])
    if polys:
        print(f"  Degree {d}:")
        for p in polys:
            print(f"      {p}")
    else:
        print(f"  Degree {d}: (no new generators)")
print("### CELL 11 END ###", flush=True)

print("### CELL 12 START ###", flush=True)
# ============================================================
# TEST 1: Homogeneous derivation (correct grading)
# ============================================================
print("=== TEST 1: D(X)=0, D(Y)=X, D(Z)=Y ===")
R.<X,Y,Z> = PolynomialRing(QQ, order='degrevlex')
D = LND([R(0), X, Y])
weights = (1,1,1)      # standard grading
D_degree = 0           # D preserves total degree
max_deg = 5

kernel = maubach_kernel(D, weights, D_degree, max_deg)
new_gens = filter_new_generators(kernel, R)

for d in range(max_deg + 1):
    polys = new_gens.get(d, [])
    print(f"Degree {d}: {len(polys)} elements")
    if polys:
        for idx, p in enumerate(polys, 1):
            print(f"  {idx}: {p}")
    else:
        print("  (none)")
print("### CELL 12 END ###", flush=True)

print("### CELL 13 START ###", flush=True)
# ============================================================
# TEST 1 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Å“ Correct D_degree = 0, print all polynomials
# ============================================================
print("=== TEST 1: D(X)=0, D(Y)=X, D(Z)=Y ===")
R.<X,Y,Z> = PolynomialRing(QQ, order='degrevlex')
D = LND([R(0), X, Y])
weights = (1,1,1)
D_degree = 0   # correct: D preserves total degree
max_deg = 5
kernel = maubach_kernel(D, weights, D_degree, max_deg)
new_gens = filter_new_generators(kernel, R)

for d in range(max_deg + 1):
    polys = new_gens.get(d, [])
    print(f"Degree {d}: {len(polys)} elements")
    if polys:
        # Print each polynomial, nicely formatted
        for idx, p in enumerate(polys, 1):
            print(f"  {idx}: {p}")
    else:
        print("  (none)")
print("### CELL 13 END ###", flush=True)

print("### CELL 14 START ###", flush=True)
# ============================================================
# TEST 2: Essen generator construction
# ============================================================
print("\n=== TEST 2: D(X)=2, D(Y)=3, D(Z)=(3X-2Y)^2 ===")
R.<X,Y,Z> = PolynomialRing(QQ, order='degrevlex')
f = R(2)
g = R(3)
F_expr = g*X - f*Y
h = F_expr^2
D = LND([f, g, h])
F, G, H = essen_generators(D)

print("Essen generators:")
print(f"F = {F}")
print(f"G = {G}")
print(f"H = {H}")

print("\nVerification (all should be 0):")
for name, poly in [("F", F), ("G", G), ("H", H)]:
    Dpoly = D(poly)
    print(f"D({name}) = {Dpoly} {'ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“' if Dpoly == 0 else 'ÃƒÂ¢Ã…â€œÃ¢â‚¬â€'}")
print("### CELL 14 END ###", flush=True)

print("### CELL 15 START ###", flush=True)
# ============================================================
# TEST 3: Slice derivation (thesis example)
# ============================================================
print("\n=== TEST 3: D(X)=1, D(Y)=X, D(Z)=X^2 ===")
R.<X,Y,Z> = PolynomialRing(QQ, order='degrevlex')
D = LND([R(1), X, X^2])
weights = (1, 2, 3)      # deg(X)=1, deg(Y)=2, deg(Z)=3
D_degree = -1            # D decreases degree by 1
max_deg = 4

kernel = maubach_kernel(D, weights, D_degree, max_deg)
new_gens = filter_new_generators(kernel, R)

for d in range(max_deg + 1):
    polys = new_gens.get(d, [])
    print(f"Weighted degree {d}: {len(polys)} elements")
    if polys:
        for idx, p in enumerate(polys, 1):
            print(f"  {idx}: {p}")
    else:
        print("  (none)")
print("### CELL 15 END ###", flush=True)

print("### CELL 16 START ###", flush=True)
# ============================================================
# TEST 3: Slice derivation (thesis example)
# ============================================================
print("\n=== TEST 3: D(X)=1, D(Y)=X, D(Z)=X^2 ===")
R.<X,Y,Z> = PolynomialRing(QQ, order='degrevlex')
D = LND([R(1), X, X^2])
weights = (1, 2, 3)      # deg(X)=1, deg(Y)=2, deg(Z)=3
D_degree = -1            # D decreases degree by 1
max_deg = 4

kernel = maubach_kernel(D, weights, D_degree, max_deg)

for d in range(max_deg + 1):
    polys = kernel.get(d, [])
    print(f"Weighted degree {d}: {len(polys)} elements")
    if polys:
        for idx, p in enumerate(polys, 1):
            print(f"  {idx}: {p}")
    else:
        print("  (none)")
print("### CELL 16 END ###", flush=True)

print("### CELL 17 START ###", flush=True)
# ================================================================
# TEST 3: Slice derivation
# ================================================================
print("\n=== TEST 3: D(X)=1, D(Y)=X, D(Z)=X^2 ===")
R.<X,Y,Z> = PolynomialRing(QQ, order='degrevlex')
D3 = LND([R(1), X, X^2])
weights = (1,2,3)
D_degree = -1
max_deg = 4

kernel3 = maubach_kernel(D3, weights, D_degree, max_deg)
new_gens3 = filter_new_generators(kernel3, R)

print("New generators by degree (minimal generating set):")
for d in range(max_deg + 1):
    polys = new_gens3.get(d, [])
    if polys:
        print(f"  Degree {d}:")
        for p in polys:
            print(f"      {p}")
    else:
        print(f"  Degree {d}: (no new generators)")
print("### CELL 17 END ###", flush=True)

print("### CELL 18 START ###", flush=True)
# ================================================================
# TEST 4 (corrected): D(X)=0, D(Y)=X, D(Z)=X^2
# ================================================================
print("\n=== TEST 4: D(X)=0, D(Y)=X, D(Z)=X^2 ===")
R.<X,Y,Z> = PolynomialRing(QQ, order='degrevlex')
D4 = LND([R(0), X, X^2])

# Correct grading: deg(X)=1, deg(Y)=1, deg(Z)=2, D_degree = 0
weights = (1, 1, 2)
D_degree = 0
max_deg = 5

kernel4 = maubach_kernel(D4, weights, D_degree, max_deg)
new_gens4 = filter_new_generators(kernel4, R)

print("New generators by degree (minimal generating set):")
for d in range(max_deg + 1):
    polys = new_gens4.get(d, [])
    if polys:
        print(f"  Degree {d}:")
        for p in polys:
            print(f"      {p}")
    else:
        print(f"  Degree {d}: (no new generators)")

# Verification
print("\nVerification:")
print(f"  D(X) = {D4(X)}")
for d, polys in new_gens4.items():
    for p in polys:
        print(f"  D({p}) = {D4(p)} {'ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“' if D4(p)==0 else 'ÃƒÂ¢Ã…â€œÃ¢â‚¬â€'}")
print("### CELL 18 END ###", flush=True)

print("### CELL 19 START ###", flush=True)
# ================================================================
# Example 4.1: D(X)=a, D(Y)=b, D(Z)=c*X^2 + a*Y^2   (non-homogenous)
# ================================================================
print("=== Example 4.1 ===")
R.<a,b,c,X,Y,Z> = PolynomialRing(QQ, order='degrevlex')

D = LND([R(0), R(0), R(0), a, b, c*X^2 + a*Y^2])

weight_vecs = [
    (1, 0),   # a
    (1, 1),   # b
    (1, 2),   # c
    (1, 1),   # X
    (1, 2),   # Y
    (3, 5),   # Z
]
D_degree_vec = (0, -1)
max_degree_vec = (10, 10)

kernel = maubach_kernel_multigraded(D, weight_vecs, D_degree_vec, max_degree_vec)
new_gens = filter_new_generators_multigraded(kernel, R, main_vars=[X, Y, Z])

print("New generators (minimal set):")
for d in sorted(new_gens.keys(), key=lambda x: (sum(x), x)):
    print(f"  Degree {d}:")
    for p in new_gens[d]:
        print(f"      {p}")
print("### CELL 19 END ###", flush=True)

print("### CELL 20 START ###", flush=True)
# ================================================================
# Example 4.2: D(X)=U, D(Y)=V, D(Z)=X^2 + Y^3   (non-homogenous)
# ================================================================
print("\n=== Example 4.2 ===")
R.<U,V,X,Y,Z> = PolynomialRing(QQ, order='degrevlex')

D = LND([R(0), R(0), U, V, X^2 + Y^3])

weight_vecs = [
    (3, 2),   # U
    (2, 1),   # V
    (3, 3),   # X
    (2, 2),   # Y
    (6, 7),   # Z
]
D_degree_vec = (0, -1)
max_degree_vec = (15, 13)

kernel = maubach_kernel_multigraded(D, weight_vecs, D_degree_vec, max_degree_vec)
new_gens = filter_new_generators_multigraded(kernel, R, main_vars=[X, Y, Z])

print("New generators (minimal set):")
for d in sorted(new_gens.keys(), key=lambda x: (sum(x), x)):
    print(f"  Degree {d}:")
    for p in new_gens[d]:
        print(f"      {p}")
print("### CELL 20 END ###", flush=True)

print("### CELL 21 START ###", flush=True)
# ================================================================
# Homogeneous Test: D(X)=0, D(Y)=X, D(Z)=X+Y
# ================================================================
print("=== Homogeneous Test: D(X)=0, D(Y)=X, D(Z)=X+Y ===")
R.<X,Y,Z> = PolynomialRing(QQ, order='degrevlex')
D = LND([R(0), X, X + Y])

weights = (1, 1, 1)   # standard grading
D_degree = 0          # D preserves total degree
max_deg = 5

kernel = maubach_kernel(D, weights, D_degree, max_deg)
new_gens = filter_new_generators(kernel, R)

print("New generators (minimal set):")
for d in range(max_deg + 1):
    polys = new_gens.get(d, [])
    if polys:
        print(f"  Degree {d}:")
        for p in polys:
            print(f"      {p}")
    else:
        print(f"  Degree {d}: (no new generators)")

# Verification
print("\nVerification:")
for d, polys in new_gens.items():
    for p in polys:
        print(f"  D({p}) = {D(p)} {'ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“' if D(p)==0 else 'ÃƒÂ¢Ã…â€œÃ¢â‚¬â€'}")
print("### CELL 21 END ###", flush=True)

print("### CELL 22 START ###", flush=True)
def run_df_example(m, r, max_degree_vec):
    """Run multi-graded Maubach on D_(m,r)."""
    print(f"\n=== D_({m},{r}) on k[X,S,T,U,V] ===")
    R.<X,S,T,U,V> = PolynomialRing(QQ, order='degrevlex')

    # Derivation: X^m ÃƒÂ¢Ã‹â€ Ã¢â‚¬Å¡_S + S ÃƒÂ¢Ã‹â€ Ã¢â‚¬Å¡_T + T ÃƒÂ¢Ã‹â€ Ã¢â‚¬Å¡_U + X^r ÃƒÂ¢Ã‹â€ Ã¢â‚¬Å¡_V
    D = LND([R(0), X^m, S, T, X^r])

    # Multi-grading from the thesis
    weight_vecs = [
        (1, 0),   # X
        (m, 1),   # S
        (m, 2),   # T
        (m, 3),   # U
        (r, 1),   # V
    ]
    D_degree_vec = (0, -1)

    kernel = maubach_kernel_multigraded(D, weight_vecs, D_degree_vec, max_degree_vec)
    #new_gens = filter_new_generators_multigraded(kernel, R, main_vars=[X, S, T, U, V])

    print(f"New generators (minimal set): {sum(len(v) for v in kernel.values())} elements")
    for d in sorted(kernel.keys(), key=lambda x: (sum(x), x)):
        for p in kernel[d]:
            print(f"  [{d}]  {p}")
    return kernel
print("### CELL 22 END ###", flush=True)

print("### CELL 23 START ###", flush=True)
run_df_example(3, 1, (12, 4))
print("### CELL 23 END ###", flush=True)

print("### CELL 24 START ###", flush=True)
run_df_example(4, 2, (12, 4))
print("### CELL 24 END ###", flush=True)

print("### CELL 25 START ###", flush=True)
run_df_example(4, 2, (12, 4))
print("### CELL 25 END ###", flush=True)

print("### CELL 26 START ###", flush=True)
def subalgebra_membership(f_list, g, R=None):
    """
    Check whether g belongs to the subalgebra k[f_list] ÃƒÂ¢Ã…Â Ã¢â‚¬Â  R.
    Implements Proposition C.2.3 (van den Essen).
    Returns (True, P) if g = P(f_list) for some P, else (False, None).
    """
    if R is None:
        R = g.parent()
    n = R.ngens()
    x_names = list(R.variable_names())
    m = len(f_list)
    y_names = [f'Y_{i}' for i in range(m)]

    # New ring with X's first (higher priority) then Y's ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â eliminates X's in GB
    S = PolynomialRing(R.base_ring(), names=x_names + y_names, order='lex')
    Xs = S.gens()[:n]
    Ys = S.gens()[n:]

    # Ideal J = (Y_i - f_i)
    J = S.ideal([Ys[i] - S(f_list[i]) for i in range(m)])
    G = J.groebner_basis()

    # Normal form of g modulo J
    nf = S(g).reduce(G)

    # g ÃƒÂ¢Ã‹â€ Ã‹â€  k[f_list]  iff  nf involves only Y variables
    if all(nf.degree(x) == 0 for x in Xs):
        return True, nf
    return False, None
print("### CELL 26 END ###", flush=True)

print("### CELL 27 START ###", flush=True)
def filter_via_subalgebra(kernel_by_degree, R, main_vars, verbose=True):
    """
    Fast multi-tier filter:
      1) skip constants and base-ring elements
      2) skip if divisible by any single existing generator (FAST)
      3) subalgebra membership check (SLOW, but correct)
    """
    def sort_key(d):
        if isinstance(d, tuple):
            return (sum(d), d)
        return (d, (d,))

    generators = []
    new_by_degree = {}
    sorted_degrees = sorted(kernel_by_degree.keys(), key=sort_key)

    total_deg = len(sorted_degrees)
    for i, d in enumerate(sorted_degrees):
        if verbose:
            print(f"  [Filter] degree {d}  ({i+1}/{total_deg})  "
                  f"candidates = {len(kernel_by_degree[d])}")
        new_at_d = []
        for p in kernel_by_degree[d]:
            if p == 1:
                continue
            if not any(p.degree(v) > 0 for v in main_vars):
                continue

            is_new = True

            if generators:
                # --- Fast check: divisible by any single generator ---
                for (_, g) in generators:
                    if g.divides(p):
                        is_new = False
                        break

                # --- Slow check: subalgebra membership ---
                if is_new:
                    gen_polys = [g for (_, g) in generators]
                    try:
                        in_subalg, _ = subalgebra_membership(gen_polys, p, R)
                        if in_subalg:
                            is_new = False
                    except Exception as e:
                        # If Singular fails, fall back to "not new"
                        # (conservative: avoids false positives)
                        pass

            if is_new:
                new_at_d.append(p)
                generators.append((d, p))

        if new_at_d:
            new_by_degree[d] = new_at_d
            if verbose:
                print(f"    ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“ NEW GENERATORS at {d}:")
                for p in new_at_d:
                    print(f"        {p}")

    return new_by_degree
print("### CELL 27 END ###", flush=True)

print("### CELL 28 START ###", flush=True)
def subalgebra_membership(f_list, g, R=None):
    if R is None: R = g.parent()
    n = R.ngens()
    x_names = list(R.variable_names())
    m = len(f_list)
    y_names = [f'Y_{i}' for i in range(m)]
    S = PolynomialRing(R.base_ring(), names=x_names + y_names, order='lex')
    Xs = S.gens()[:n]
    Ys = S.gens()[n:]
    J = S.ideal([Ys[i] - S(f_list[i]) for i in range(m)])
    G = J.groebner_basis()
    nf = S(g).reduce(G)
    result = all(nf.degree(x) == 0 for x in Xs)
    # Restore active ring in Singular
    try:
        R.ring()
    except Exception:
        pass
    if result:
        return True, nf
    return False, None
print("### CELL 28 END ###", flush=True)

print("### CELL 29 START ###", flush=True)
# --- (4) Full pipeline ---
def maubach_with_membership(D, weights, D_degree, max_degree,
                            R, main_vars, multigraded=False):
    """
    Complete pipeline: Maubach ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ Subalgebra filter.
    """
    if multigraded:
        kernel = maubach_kernel_multigraded(D, weights, D_degree, max_degree)
    else:
        kernel = maubach_kernel(D, weights, D_degree, max_degree)

    return filter_via_subalgebra(kernel, R, main_vars)
print("### CELL 29 END ###", flush=True)

print("### CELL 30 START ###", flush=True)
# --- Test 1: D(X)=0, D(Y)=X, D(Z)=Y (single weight) ---
print("=== Test 1 (single weight) ===")
R.<X,Y,Z> = PolynomialRing(QQ, order='degrevlex')
D = LND([R(0), X, Y])
new_gens = maubach_with_membership(D, (1,1,1), 0, 5, R, [X,Y,Z], multigraded=False)
for d, polys in new_gens.items():
    print(f"  Degree {d}:")
    for p in polys:
        print(f"      {p}")
print("### CELL 30 END ###", flush=True)

print("### CELL 31 START ###", flush=True)
print("\n=== Example 4.1 (multi-graded) ===")
R.<a,b,c,X,Y,Z> = PolynomialRing(QQ, order='degrevlex')
D = LND([R(0), R(0), R(0), a, b, c*X^2 + a*Y^2])
weight_vecs = [(1,0),(1,1),(1,2),(1,1),(1,2),(3,5)]
new_gens = maubach_with_membership(D, weight_vecs, (0,-1), (10,10), R, [X,Y,Z], multigraded=True)
for d, polys in new_gens.items():
    print(f"  Degree {d}:")
    for p in polys: print(f"      {p}")
print("### CELL 31 END ###", flush=True)

print("### CELL 32 START ###", flush=True)
print("\n=== New Test (single weight, D(X)=0, D(Y)=X, D(Z)=X+Y) ===")
R.<X,Y,Z> = PolynomialRing(QQ, order='degrevlex')
D = LND([R(0), X, X+Y])
new_gens = maubach_with_membership(D, (1,1,1), 0, 5, R, [X,Y,Z], multigraded=False)
for d, polys in new_gens.items():
    print(f"  Degree {d}:")
    for p in polys: print(f"      {p}")
print("### CELL 32 END ###", flush=True)

print("### CELL 33 START ###", flush=True)
# ================================================================
# Step 2 only: Two separate lists
#   List 1 = kernel generators (in R = k[X,S,T,U,V])
#   List 2 = relations among them (in RY = k[Y_0,...,Y_{m-1}])
# ================================================================

def compute_relations(f_list, R):
    """
    f_list: list of polynomials in R (the kernel generators)
    R:      the polynomial ring (e.g. k[X,S,T,U,V])
    
    Returns: list of relations in a NEW ring k[Y0, Y1, ...]
    """
    # 1. Names of the original variables + new Y variables
    x_names = list(R.variable_names())
    n = len(x_names)
    m = len(f_list)
    y_names = [f'Y{i}' for i in range(m)]

    # 2. Big ring: k[X, S, T, U, V, Y0, Y1, ...] with lex order
    S = PolynomialRing(R.base_ring(), names=x_names + y_names, order='lex')
    X = S.gens()[:n]          # original variables
    Y = S.gens()[n:]          # new Y variables

    # 3. Ideal: J = (Y0 - f1, Y1 - f2, ..., Y_{m-1} - fm)
    J = S.ideal([Y[i] - S(f_list[i]) for i in range(m)])

    # 4. Groebner basis of J
    GB = J.groebner_basis()

    # 5. Keep only the GB elements that involve ONLY the Y variables
    #    (these are the relations)
    relations = [g for g in GB if all(g.degree(x) == 0 for x in X)]

    # 6. Move relations into a smaller ring k[Y0, Y1, ...]
    RY = PolynomialRing(R.base_ring(), names=y_names, order='degrevlex')
    relations = [RY(r) for r in relations]

    return relations, RY
print("### CELL 33 END ###", flush=True)

print("### CELL 34 START ###", flush=True)
def run_df_single(m, r, max_deg, verbose_filter=True):
    print(f"\n{'='*60}")
    print(f"=== D_({m},{r}) ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â single-weight ===")
    print(f"{'='*60}")
    R.<X,S,T,U,V> = PolynomialRing(QQ, order='degrevlex')
    # Derivation: X^m ÃƒÂ¢Ã‹â€ Ã¢â‚¬Å¡_S + S ÃƒÂ¢Ã‹â€ Ã¢â‚¬Å¡_T + T ÃƒÂ¢Ã‹â€ Ã¢â‚¬Å¡_U + X^r ÃƒÂ¢Ã‹â€ Ã¢â‚¬Å¡_V
    D = LND([R(0), X^m, S, T, X^r])

    # Single weight that makes D homogeneous
    weights = (1, m, m, m, r)   # X, S, T, U, V
    D_degree = 0                # all shifts are 0

    # Step 1: Maubach
    print(f"Weights: {weights},  D_degree = {D_degree},  max_deg = {max_deg}")
    kernel = maubach_kernel(D, weights, D_degree, max_deg)
    
    #step 2 code for finding elements in kernel
    total = sum(len(v) for v in kernel.values())

    """# Step 2: Filter using algebra membership
    new_gens = filter_via_subalgebra(kernel, R, [X, S, T, U, V],
                                     verbose=verbose_filter)"""

    # Step 3: Relation 
    # --- Step 2: Compute relation ideal for the CURRENT generators ---
    # Start with an empty generator set; we will grow it as we find new elements.
    # LIST 1 (in R = k[X,S,T,U,V])
    # ------------------------------------------------------------
    # LIST 1: kernel generators (polynomials in R = k[X,S,T,U,V])
    # LIST 2: relations among them (polynomials in RY)
    # ------------------------------------------------------------
    kernel_generators = []          # LIST 1  (in R)
    relation_generators = []        # LIST 2  (in RY)
    RY = None                       # the ring where LIST 2 lives

    # Take the FIRST kernel element at each degree as a generator candidate.
    # (This is a simple selection; no membership test yet.)
    for d in sorted(kernel.keys()):
        for p in kernel[d]:
            if p == 1:
                continue
            kernel_generators.append(p)

    # ------------------------------------------------------------
    # Step 2: compute the relation ideal for List 1
    # ------------------------------------------------------------
    print(f"\n[Step 2] Relations among the {len(kernel_generators)} generators")
    if len(kernel_generators) >= 2:
        relation_generators, RY = compute_relations(kernel_generators, R)

    # ------------------------------------------------------------
    # Print the TWO lists separately
    # ------------------------------------------------------------
    print("\n" + "="*60)
    print("LIST 1 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Kernel generators  (in k[X,S,T,U,V])")
    print("="*60)
    for i, p in enumerate(kernel_generators, 1):
        print(f"  {i}. {p}")

    print("\n" + "="*60)
    if RY is not None:
        print(f"LIST 2 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Relations  (in {RY})")
    else:
        print("LIST 2 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Relations  (ring not defined)")
    print("="*60)
    if relation_generators:
        for i, rel in enumerate(relation_generators, 1):
            print(f"  {i}. {rel}")
    else:
        print("  (empty)")

    return kernel_generators, relation_generators, RY
        # --- Final summary ---
    print(f"\n[Final] Total minimal generators: {len(generators)}")
    for d, p in generators:
        print(f"        [deg {d}]  {p}")

    return generators, kernel        

    """print(f"\n>>> TOTAL GENERATORS: {total}")
    for d in sorted(new_gens.keys()):
        for p in new_gens[d]:
            print(f"  [deg {d}]  {p}")"""
    #return kernel, total       #new_gens
print("### CELL 34 END ###", flush=True)

print("### CELL 35 START ###", flush=True)
def run_df_single(m, r, max_deg):
    """
    Step 1: Maubach  ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ kernel elements
    Step 2: Relations ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ ker(ÃƒÂÃ¢â‚¬Â ) among them
    
    Returns two separate lists:
      List 1: kernel generators   (in k[X,S,T,U,V])
      List 2: relations           (in k[Y0,Y1,...])
    """
    print(f"\n{'='*60}")
    print(f"=== D_({m},{r}) ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Maubach + Relations ===")
    print(f"{'='*60}")

    # --- Set up ring and derivation ---
    R.<X,S,T,U,V> = PolynomialRing(QQ, order='degrevlex')
    D = LND([R(0), X^m, S, T, X^r])
    weights = (1, m, m, m, r)      # X, S, T, U, V
    D_degree = 0

    # ============================================================
    # STEP 1: Maubach
    # ============================================================
    print(f"\n[Step 1] Maubach  (weights={weights}, D_degree={D_degree})")
    kernel = maubach_kernel(D, weights, D_degree, max_deg)
    total = sum(len(v) for v in kernel.values())
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ {total} kernel elements up to degree {max_deg}")

    # ============================================================
    # LIST 1: kernel generators (in R = k[X,S,T,U,V])
    # ============================================================
    kernel_generators = []
    for d in sorted(kernel.keys()):
        for p in kernel[d]:
            if p == 1:
                continue
            kernel_generators.append(p)

    # ============================================================
    # STEP 2: Relations among List 1  ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢  LIST 2
    # ============================================================
    print(f"\n[Step 2] Computing relations among {len(kernel_generators)} generators")
    relation_generators = []
    RY = None
    if len(kernel_generators) >= 2:
        relation_generators, RY = compute_relations(kernel_generators, R)

    # ============================================================
    # Print the TWO lists separately
    # ============================================================
    print("\n" + "="*60)
    print("LIST 1 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Kernel generators  (in k[X,S,T,U,V])")
    print("="*60)
    for i, p in enumerate(kernel_generators, 1):
        print(f"  {i}. {p}")

    print("\n" + "="*60)
    if RY is not None:
        print(f"LIST 2 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Relations  (in {RY})")
    else:
        print("LIST 2 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Relations")
    print("="*60)
    if relation_generators:
        for i, rel in enumerate(relation_generators, 1):
            print(f"  {i}. {rel}")
    else:
        print("  (empty)")

    return kernel_generators, relation_generators, RY


print("### CELL 35 END ###", flush=True)

print("### CELL 36 START ###", flush=True)
def compute_relations(f_list, R):
    """
    Compute relations among f_list via GrÃƒÆ’Ã‚Â¶bner elimination.
    Defensive: uses explicit algorithm and restores ring state.
    """
    x_names = list(R.variable_names())
    n = len(x_names)
    m = len(f_list)
    y_names = [f'Y{i}' for i in range(m)]

    # Big ring: k[X,S,T,U,V, Y0, Y1, ...] with lex order
    S = PolynomialRing(R.base_ring(), names=x_names + y_names, order='lex')
    X = S.gens()[:n]
    Y = S.gens()[n:]

    # Ideal J = (Y0 - f1, Y1 - f2, ...)
    J = S.ideal([Y[i] - S(f_list[i]) for i in range(m)])

    # GrÃƒÆ’Ã‚Â¶bner basis ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â force the standard algorithm explicitly
    try:
        GB = J.groebner_basis(algorithm='singular:std')
    except Exception as e:
        print(f"  [Warning] singular:std failed: {e}")
        GB = J.groebner_basis()

    # Keep only the elements that involve ONLY the Y variables
    relations = [g for g in GB if all(g.degree(x) == 0 for x in X)]

    # Move relations into a smaller ring k[Y0, Y1, ...]
    RY = PolynomialRing(R.base_ring(), names=y_names, order='degrevlex')
    relations = [RY(r) for r in relations]

    # Restore active ring (best effort)
    try:
        R.ring()
    except Exception:
        pass

    return relations, RY
print("### CELL 36 END ###", flush=True)

print("### CELL 37 START ###", flush=True)
# ----------------------------------------------------------------
# Helper 1: build the GrÃƒÆ’Ã‚Â¶bner basis of the substitution ideal
# ----------------------------------------------------------------
def build_gb(f_list, R):
    x_names = list(R.variable_names())
    n = len(x_names)
    m = len(f_list)
    y_names = [f'Y{i}' for i in range(m)]

    S = PolynomialRing(R.base_ring(), names=x_names + y_names, order='lex')

    # *** Force Singular to activate this ring ***
    S._singular_()          # creates the Singular ring and sets it active

    Xs = S.gens()[:n]
    Ys = S.gens()[n:]

    # Convert generators into S first (separate step) to avoid state confusion
    f_in_S = [S(f) for f in f_list]

    J = S.ideal([Ys[i] - f_in_S[i] for i in range(m)])

    # Explicitly specify the algorithm to avoid the buggy default path
    GB = J.groebner_basis(algorithm='singular:std')

    return S, Xs, Ys, GB
print("### CELL 37 END ###", flush=True)

print("### CELL 38 START ###", flush=True)
def build_gb(f_list, R, Q_list=None):
    """
    Build the GrÃƒÆ’Ã‚Â¶bner basis of the substitution ideal.
    
    Parameters:
      f_list : list of generators (in R)
      R      : the polynomial ring
      Q_list : list of polynomials defining R, if R is a quotient ring.
               For a free polynomial ring, leave as None.
    
    The ideal is:
      J = (Y0 - F1, Y1 - F2, ..., Ym - Fm, Q1, Q2, ..., Qs)
    """
    if Q_list is None:
        Q_list = []      # no defining relations (free ring)
    
    x_names = list(R.variable_names())
    n = len(x_names)
    m = len(f_list)
    y_names = [f'Y{i}' for i in range(m)]
    
    # Big ring: k[X, Y]
    S = PolynomialRing(R.base_ring(), names=x_names + y_names, order='lex')
    Xs = S.gens()[:n]
    Ys = S.gens()[n:]
    
    # Substitution part: Y_i - F_i
    f_in_S = [S(f) for f in f_list]
    subst_gens = [Ys[i] - f_in_S[i] for i in range(m)]
    
    # Defining relations of R: Q_1, ..., Q_s (if any)
    q_gens = [S(q) for q in Q_list]
    
    # Combine both parts
    J = S.ideal(subst_gens + q_gens)
    
    # Compute GrÃƒÆ’Ã‚Â¶bner basis
    GB = J.groebner_basis(algorithm='libsingular:std')
    
    return S, Xs, Ys, GB
print("### CELL 38 END ###", flush=True)

print("### CELL 39 START ###", flush=True)
def build_gb(f_list, R, Q_list=None):
    """
    Build the GrÃƒÆ’Ã‚Â¶bner basis of the substitution ideal.
    J = (Y0 - F1, ..., Ym - Fm, Q1, ..., Qs)
    """
    if Q_list is None:
        Q_list = []

    x_names = list(R.variable_names())
    n = len(x_names)
    m = len(f_list)
    y_names = [f'Y{i}' for i in range(m)]

    # *** CRITICAL: rename Python variable to avoid shadowing ring variable 'S' ***
    S_ring = PolynomialRing(R.base_ring(), names=x_names + y_names, order='lex')
    Xs = S_ring.gens()[:n]
    Ys = S_ring.gens()[n:]

    f_in_S_ring = [S_ring(f) for f in f_list]
    subst_gens = [Ys[i] - f_in_S_ring[i] for i in range(m)]
    q_gens = [S_ring(q) for q in Q_list]

    J = S_ring.ideal(subst_gens + q_gens)
    GB = J.groebner_basis()

    return S_ring, Xs, Ys, GB
print("### CELL 39 END ###", flush=True)

print("### CELL 40 START ###", flush=True)
# ----------------------------------------------------------------
# Helper 2: extract relations (pure-Y elements) from GB
# ----------------------------------------------------------------
def relations_from_gb(GB, Xs, y_names, R):
    rels = [g for g in GB if all(g.degree(x) == 0 for x in Xs)]
    RY = PolynomialRing(R.base_ring(), names=y_names, order='degrevlex')
    RY._singular_()                          # force activation
    rels_RY = [RY(r) for r in rels]
    return rels_RY, RY
print("### CELL 40 END ###", flush=True)

print("### CELL 41 START ###", flush=True)
# ----------------------------------------------------------------
# Helper 3: membership test using a precomputed GB
def membership_using_gb(GB, Xs, g, S):
    nf = S(g).reduce(GB)
    is_member = all(nf.degree(x) == 0 for x in Xs)
    return (True, nf) if is_member else (False, None)
print("### CELL 41 END ###", flush=True)

print("### CELL 42 START ###", flush=True)
def substitute_Y_with_F(rel, LIST1, RY, R):
    """Manually substitute Y_i -> F_i (avoids ring coercion)."""
    result = R(0)
    Ys = RY.gens()
    for coeff, monomial in zip(rel.coefficients(), rel.monomials()):
        term = R(coeff)
        for i, Yv in enumerate(Ys):
            deg = monomial.degree(Yv)
            if deg > 0:
                term *= LIST1[i] ** deg
        result += term
    return result
print("### CELL 42 END ###", flush=True)

print("### CELL 43 START ###", flush=True)
# ================================================================
# FULL PIPELINE ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â SageMath
# Maubach kernel -> membership filter -> relations -> verification
# ================================================================

def to_big_ring(f, R, BigRing):
    """
    Convert a polynomial from R into BigRing by mapping the
    original variables to the first n variables of BigRing.
    Avoid implicit cross-ring polynomial coercions.
    """
    n = R.ngens()
    X_alt = BigRing.gens()[:n]

    result = BigRing.zero()

    for coeff, mon in zip(f.coefficients(), f.monomials()):
        term = BigRing(coeff)

        for j, orig_var in enumerate(R.gens()):
            exponent = mon.degree(orig_var)
            if exponent:
                term *= X_alt[j] ** exponent

        result += term

    return result


def build_gb(f_list, R, Q_list=None):
    """
    Build the ideal
        J = <Y_i - F_i, Q_j>
    in a lexicographic polynomial ring.

    Original-variable copies come first, followed by Y variables.
    """
    if Q_list is None:
        Q_list = []

    n = R.ngens()
    k = len(f_list)

    if len(Q_list):
        Q_list = list(Q_list)
    else:
        Q_list = []

    copies = ['a', 'b', 'c', 'd', 'e'][:n]
    y_names = [f'Y{i}' for i in range(k)]

    BigRing = PolynomialRing(
        QQ,
        names=copies + y_names,
        order='lex'
    )

    X_alt = BigRing.gens()[:n]
    Y_alt = BigRing.gens()[n:]

    f_big = [to_big_ring(f, R, BigRing) for f in f_list]
    q_big = [to_big_ring(q, R, BigRing) for q in Q_list]

    generators = [
        Y_alt[i] - f_big[i]
        for i in range(k)
    ] + q_big

    J = BigRing.ideal(generators)
    GB = J.groebner_basis(algorithm='libsingular:std')

    return BigRing, X_alt, Y_alt, GB


def relations_from_gb(GB, X_alt, y_names, R, BigRing):
    """
    Extract pure-Y polynomials from the elimination GrÃƒÆ’Ã‚Â¶bner basis
    and manually convert them into a polynomial ring RY.
    """
    RY = PolynomialRing(
        R.base_ring(),
        names=y_names,
        order='lex'
    )

    Y_RY = RY.gens()
    n = len(X_alt)

    rels_big = [
        g for g in GB
        if all(g.degree(x) == 0 for x in X_alt)
        and g != 0
    ]

    rels_RY = []

    for rel_big in rels_big:
        poly = RY.zero()

        for coeff, mon in zip(
            rel_big.coefficients(),
            rel_big.monomials()
        ):
            term = RY(coeff)

            for j, Y_big in enumerate(BigRing.gens()[n:]):
                exponent = mon.degree(Y_big)
                if exponent:
                    term *= Y_RY[j] ** exponent

            poly += term

        if poly != 0:
            rels_RY.append(poly)

    return rels_RY, RY


def membership_using_gb(GB, X_alt, g, BigRing, R):
    """
    Test whether g belongs to the subalgebra generated by LIST1,
    modulo the constraints encoded in GB.

    A zero remainder means membership. A nonzero remainder that
    contains only Y variables also means membership.
    """
    g_big = to_big_ring(g, R, BigRing)
    nf = g_big.reduce(GB)

    is_member = all(nf.degree(x) == 0 for x in X_alt)

    return is_member, nf if is_member else None


def substitute_Y_with_F(rel, LIST1, RY, R):
    """
    Substitute Y_i -> LIST1[i].

    Use explicit monomial construction rather than cross-ring
    coercion of the whole polynomial.
    """
    if RY.ngens() != len(LIST1):
        raise ValueError(
            f"Expected {len(LIST1)} Y variables; "
            f"got {RY.ngens()}."
        )

    result = R.zero()
    Ys = RY.gens()

    for coeff, mon in zip(rel.coefficients(), rel.monomials()):
        term = R(coeff)

        for i, Yv in enumerate(Ys):
            exponent = mon.degree(Yv)

            if exponent:
                term *= LIST1[i] ** exponent

        result += term

    return result
print("### CELL 43 END ###", flush=True)

print("### CELL 44 START ###", flush=True)
"""def run_pipeline(m, r, max_deg, Q_list=None):
    print(f"\n{'='*60}")
    print(f"=== D_({m},{r})  ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â  Full Pipeline ===")
    if Q_list:
        print(f"    Q_list (Step 2 only) = {Q_list}")
    print(f"{'='*60}")

    R.<X,S,T,U,V> = PolynomialRing(QQ, order='degrevlex')
    D = LND([R(0), X^m, S, T, X^r])
    weights = (1, m, m, m, r)
    D_degree = 0

    # --- Step 1: Maubach ---
    print(f"\n[Step 1] Maubach")
    kernel = maubach_kernel(D, weights, D_degree, max_deg)
    total = sum(len(v) for v in kernel.values())
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ {total} raw kernel elements")

    # --- Step 1.5: Membership filter (Q_list NOT used here) ---
    print(f"\n[Step 1.5] Membership filter (free ring)")
    LIST1 = []
    for d in sorted(kernel.keys()):
        for p in kernel[d]:
            if p == 1:
                continue
            if LIST1:
                S_ring, Xs, Ys, GB = build_gb(LIST1, R, Q_list=None)  # ÃƒÂ¢Ã¢â‚¬Â Ã‚Â NO Q
                is_member, _ = membership_using_gb(GB, Xs, p, S_ring)
                if is_member:
                    continue
            LIST1.append(p)
            print(f"        ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“ NEW [deg {d}]  {p}")
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ LIST 1: {len(LIST1)} minimal generators")

    # --- Step 2: Relations (Q_list IS used here) ---
    print(f"\n[Step 2] Relations among LIST 1 (with Q_list)")
    LIST2 = []
    RY = None
    if len(LIST1) >= 2:
        S_ring, Xs, Ys, GB = build_gb(LIST1, R, Q_list=Q_list)  # ÃƒÂ¢Ã¢â‚¬Â Ã‚Â USE Q
        y_names = [str(y) for y in Ys]
        LIST2, RY = relations_from_gb(GB, Xs, y_names, R)
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ LIST 2: {len(LIST2)} relations")

    # --- Step 3: Verify completeness (Q_list NOT used here) ---
    print(f"\n[Step 3] Verify completeness (free ring)")
    S_ring, Xs, Ys, GB = build_gb(LIST1, R, Q_list=None)  # ÃƒÂ¢Ã¢â‚¬Â Ã‚Â NO Q
    missed = []
    for d in sorted(kernel.keys()):
        for p in kernel[d]:
            if p == 1 or p in LIST1:
                continue
            is_member, _ = membership_using_gb(GB, Xs, p, S_ring)
            if not is_member:
                missed.append(p)
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ {len(missed)} missed" if missed else "        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ STOP ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“")

    # --- Print ---
    print("\n" + "="*60)
    print("LIST 1 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Generators")
    print("="*60)
    for i, p in enumerate(LIST1, 1):
        print(f"  {i}. {p}")

    print("\n" + "="*60)
    print(f"LIST 2 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Relations  (in {RY})")
    print("="*60)
    if LIST2:
        for i, rel in enumerate(LIST2, 1):
            print(f"  {i}. {rel}")
    else:
        print("  (empty)")

    return LIST1, LIST2, RY"""
"""def run_pipeline(m, r, max_deg, Q_list=None):
    print(f"\n{'='*60}")
    print(f"=== D_({m},{r})  ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â  Full Pipeline ===")
    if Q_list:
        print(f"    Q_list (Step 2 only) = {Q_list}")
    print(f"{'='*60}")

    R.<X,S,T,U,V> = PolynomialRing(QQ, order='degrevlex')
    D = LND([R(0), X^m, S, T, X^r])
    weights = (1, m, m, m, r)
    D_degree = 0

    # --- Step 1: Maubach ---
    print(f"\n[Step 1] Maubach")
    kernel = maubach_kernel(D, weights, D_degree, max_deg)
    total = sum(len(v) for v in kernel.values())
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ {total} raw kernel elements")

    # --- Step 1.5: Membership filter (free ring) ---
    print(f"\n[Step 1.5] Membership filter (free ring)")
    LIST1 = []
    for d in sorted(kernel.keys()):
        for p in kernel[d]:
            if p == 1:
                continue
            if LIST1:
                S_ring, Xs, Ys, GB = build_gb(LIST1, R, Q_list=None)
                is_member, _ = membership_using_gb(GB, Xs, p, S_ring)
                if is_member:
                    continue
            LIST1.append(p)
            print(f"        ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“ NEW [deg {d}]  {p}")
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ LIST 1: {len(LIST1)} minimal generators")

    # --- Step 2: Relations (with Q_list) ---
    print(f"\n[Step 2] Relations among LIST 1 (with Q_list)")
    LIST2 = []
    RY = None
    if len(LIST1) >= 2:
        S_ring, Xs, Ys, GB = build_gb(LIST1, R, Q_list=Q_list)
        y_names = [str(y) for y in Ys]
        LIST2, RY = relations_from_gb(GB, Xs, y_names, R)
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ LIST 2: {len(LIST2)} relations")

    # =========================================================
    # STEP 2.5 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Substitute Y_i ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ F_i in each relation
    # =========================================================
    print(f"\n[Step 2.5] Substituting Y_i ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ F_i in each relation")
    print(f"        (result shown in R and R/(X))")

    # Build substitution map: Y_i ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ F_i
    # RY has variables Y0, Y1, ..., Y_{k-1} where k = len(LIST1)
    Y_vars = RY.gens()
    subs_map = {Y_vars[i]: LIST1[i] for i in range(len(LIST1))}

    relations_in_R = []       # relations after substitution (free ring)
    relations_in_Rmod = []    # relations after substitution + X=0 (quotient)

    for i, rel in enumerate(LIST2, 1):
        try:
            # Apply substitution Y_i ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ F_i
            result = rel.subs(subs_map)
            result = R(result)      # ensure it's in R

            # Also compute modulo X = 0 (for quotient)
            result_mod = result.subs({X: 0})

            relations_in_R.append(result)
            relations_in_Rmod.append(result_mod)

            print(f"\n  Relation {i}: {rel}")
            print(f"    In R:        {result}")
            print(f"    In R/(X):    {result_mod}")
            if result_mod == 0:
                print(f"    ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“ Vanishes in quotient R/(X)")
            else:
                print(f"    ÃƒÂ¢Ã…Â¡Ã‚Â  Nonzero in quotient")
        except Exception as e:
            print(f"\n  Relation {i}: {rel}")
            print(f"    [Substitution failed: {e}]")

    # --- Step 3: Verify completeness (free ring) ---
    print(f"\n[Step 3] Verify completeness (free ring)")
    S_ring, Xs, Ys, GB = build_gb(LIST1, R, Q_list=None)
    missed = []
    for d in sorted(kernel.keys()):
        for p in kernel[d]:
            if p == 1 or p in LIST1:
                continue
            is_member, _ = membership_using_gb(GB, Xs, p, S_ring)
            if not is_member:
                missed.append(p)
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ {len(missed)} missed" if missed else "        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ STOP ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“")

    # --- Print final lists ---
    print("\n" + "="*60)
    print("LIST 1 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Generators")
    print("="*60)
    for i, p in enumerate(LIST1, 1):
        print(f"  {i}. {p}")

    print("\n" + "="*60)
    print(f"LIST 2 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Relations  (in {RY})")
    print("="*60)
    if LIST2:
        for i, rel in enumerate(LIST2, 1):
            print(f"  {i}. {rel}")
    else:
        print("  (empty)")

    return LIST1, LIST2, RY, relations_in_R, relations_in_Rmod"""
"""def run_pipeline(m, r, max_deg, Q_list=None):
    print(f"\n{'='*60}")
    print(f"=== D_({m},{r})  ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â  Full Pipeline ===")
    if Q_list:
        print(f"    Q_list (Step 2 only) = {Q_list}")
    print(f"{'='*60}")

    R.<X,S,T,U,V> = PolynomialRing(QQ, order='degrevlex')
    D = LND([R(0), X^m, S, T, X^r])
    weights = (1, m, m, m, r)
    D_degree = 0

    # --- Step 1 ---
    print(f"\n[Step 1] Maubach")
    kernel = maubach_kernel(D, weights, D_degree, max_deg)
    total = sum(len(v) for v in kernel.values())
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ {total} raw kernel elements")

    # --- Step 1.5 ---
    print(f"\n[Step 1.5] Membership filter (free ring)")
    LIST1 = []
    for d in sorted(kernel.keys()):
        for p in kernel[d]:
            if p == 1:
                continue
            if LIST1:
                S_ring, Xs, Ys, GB = build_gb(LIST1, R, Q_list=None)
                is_member, _ = membership_using_gb(GB, Xs, p, S_ring)
                if is_member:
                    continue
            LIST1.append(p)
            print(f"        ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“ NEW [deg {d}]  {p}")
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ LIST 1: {len(LIST1)} minimal generators")

    # --- Step 2 ---
    print(f"\n[Step 2] Relations among LIST 1 (with Q_list)")
    LIST2 = []
    RY = None
    if len(LIST1) >= 2:
        S_ring, Xs, Ys, GB = build_gb(LIST1, R, Q_list=Q_list)
        y_names = [str(y) for y in Ys]
        LIST2, RY = relations_from_gb(GB, Xs, y_names, R)
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ LIST 2: {len(LIST2)} relations")

    # =========================================================
    # STEP 2.5 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Substitute + Factor
    # =========================================================
    print(f"\n[Step 2.5] Substituting Y_i ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ F_i and factoring")

    Y_vars = RY.gens()
    subs_map = {Y_vars[i]: LIST1[i] for i in range(len(LIST1))}

    relations_in_R = []
    relations_factored = []
    relations_in_Rmod = []

    for i, rel in enumerate(LIST2, 1):
        try:
            result = R(rel.subs(subs_map))
            factored = result.factor()

            # Minimum exponent of X across all monomials
            if result != 0:
                min_X = min(mm.degree(X) for mm in result.monomials())
                if min_X > 0:
                    quot = result / (X ** min_X)
                    x_power_str = f"X^{min_X} * ({quot})"
                else:
                    x_power_str = str(result)
            else:
                min_X = 0
                x_power_str = "0"

            result_mod = result.subs({X: 0})

            relations_in_R.append(result)
            relations_factored.append(factored)
            relations_in_Rmod.append(result_mod)

            print(f"\n  Relation {i}: {rel}")
            print(f"    Raw:         {result}")
            print(f"    Factored:    {factored}")
            print(f"    X-power:     {x_power_str}")
            print(f"    In R/(X):    {result_mod}")
            if result_mod == 0:
                print(f"    ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“ Vanishes in quotient R/(X)")
            else:
                print(f"    ÃƒÂ¢Ã…Â¡Ã‚Â  Nonzero in quotient")
        except Exception as e:
            print(f"\n  Relation {i}: {rel}")
            print(f"    [Failed: {e}]")

    # Summary
    print(f"\n  [Summary of X-powers]")
    for i, res_R in enumerate(relations_in_R, 1):
        if res_R != 0:
            min_X = min(mm.degree(X) for mm in res_R.monomials())
        else:
            min_X = "ÃƒÂ¢Ã‹â€ Ã…Â¾"
        print(f"    Relation {i} vanishes to order X^{min_X} at X = 0")

    # --- Step 3 ---
    print(f"\n[Step 3] Verify completeness (free ring)")
    S_ring, Xs, Ys, GB = build_gb(LIST1, R, Q_list=None)
    missed = []
    for d in sorted(kernel.keys()):
        for p in kernel[d]:
            if p == 1 or p in LIST1:
                continue
            is_member, _ = membership_using_gb(GB, Xs, p, S_ring)
            if not is_member:
                missed.append(p)
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ {len(missed)} missed" if missed else "        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ STOP ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“")

    # --- Print final lists ---
    print("\n" + "="*60)
    print("LIST 1 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Generators")
    print("="*60)
    for i, p in enumerate(LIST1, 1):
        print(f"  {i}. {p}")

    print("\n" + "="*60)
    print(f"LIST 2 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Relations  (in {RY})")
    print("="*60)
    if LIST2:
        for i, rel in enumerate(LIST2, 1):
            print(f"  {i}. {rel}")
    else:
        print("  (empty)")

    return LIST1, LIST2, RY, relations_in_R, relations_factored, relations_in_Rmod"""
# ================================================================
# FULL PIPELINE ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â corrected
# Maubach -> membership filter -> relations -> substitute/factor/normalize/normal-form
# ================================================================

def to_big_ring(f, R, BigRing):
    """Manual conversion: R -> BigRing via monomial walking."""
    n = R.ngens()
    X_alt = BigRing.gens()[:n]
    result = BigRing.zero()
    for coeff, mon in zip(f.coefficients(), f.monomials()):
        term = BigRing(coeff)
        for j, orig_var in enumerate(R.gens()):
            e = mon.degree(orig_var)
            if e:
                term *= X_alt[j] ** e
        result += term
    return result


_CURRENT_WEIGHTS = None
_GB_CACHE = {}

def build_gb(f_list, R, Q_list=None):
    if Q_list is None:
        Q_list = []
    n = R.ngens()
    k = len(f_list)
    use_modstd = (k >= 10)
    key = (tuple(str(f) for f in f_list), tuple(str(q) for q in Q_list), use_modstd)
    if key in _GB_CACHE:
        print(f"        BUILD_GB_CACHED k={k} modstd={use_modstd}", flush=True)
        return _GB_CACHE[key]
    names = ['a', 'b', 'c', 'd', 'e'][:n] + [f'Y{i}' for i in range(k)]
    if use_modstd:
        weights = _CURRENT_WEIGHTS if _CURRENT_WEIGHTS is not None else (1,) * n
        def _wdeg(f):
            ds = [sum(e * w for e, w in zip(mon.degrees(), weights)) for mon in f.monomials()]
            return max(ds) if ds else 0
        y_weights = tuple(_wdeg(f) for f in f_list)
        from sage.rings.polynomial.term_order import TermOrder
        order = TermOrder('wdegrevlex', weights) + TermOrder('wdegrevlex', y_weights)
        print(f"        BUILD_GB_MODSTD k={k} Xw={weights} Yw={y_weights}", flush=True)
    else:
        order = f'degrevlex({n}),degrevlex({k})'
        print(f"        BUILD_GB_SLIMGB k={k}", flush=True)
    BigRing = PolynomialRing(QQ, names=names, order=order)
    X_alt = BigRing.gens()[:n]
    Y_alt = BigRing.gens()[n:]
    f_big = [to_big_ring(f, R, BigRing) for f in f_list]
    q_big = [to_big_ring(q, R, BigRing) for q in Q_list]
    J = BigRing.ideal([Y_alt[i] - f_big[i] for i in range(k)] + q_big)
    t0 = time.time()
    if use_modstd:
        from sage.libs.singular.function import singular_function, lib as singular_lib
        singular_lib('modstd.lib')
        modStd = singular_function('modStd')
        GB = list(modStd(J))
    else:
        GB = J.groebner_basis(algorithm='libsingular:slimgb')
    print(f"        BUILD_GB_DONE k={k} modstd={use_modstd} time={time.time()-t0:.2f} len={len(GB)}", flush=True)
    result = (BigRing, X_alt, Y_alt, GB)
    _GB_CACHE[key] = result
    return result


def relations_from_gb(GB, X_alt, y_names, R, BigRing):
    """Extract pure-Y relations and move to RY = QQ[Y0,...,Y_{k-1}]."""
    RY = PolynomialRing(R.base_ring(), names=y_names, order='invlex')
    Y_RY = RY.gens()
    n = len(X_alt)
    Y_big = BigRing.gens()[n:]

    rels_big = [g for g in GB if all(g.degree(x) == 0 for x in X_alt) and g != 0]

    rels_RY = []
    for rel_big in rels_big:
        poly = RY.zero()
        for coeff, mon in zip(rel_big.coefficients(), rel_big.monomials()):
            term = RY(coeff)
            for j, Y_b in enumerate(Y_big):
                e = mon.degree(Y_b)
                if e:
                    term *= Y_RY[j] ** e
            poly += term
        if poly != 0:
            rels_RY.append(poly)
    return rels_RY, RY


def to_ry(r, RY, BigRing, X_alt):
    """Convert a pure-Y polynomial from BigRing to RY."""
    n = len(X_alt)
    Y_big = BigRing.gens()[n:]
    Y_RY = RY.gens()
    poly = RY.zero()
    for coeff, mon in zip(r.coefficients(), r.monomials()):
        term = RY(coeff)
        for j, Y_b in enumerate(Y_big):
            e = mon.degree(Y_b)
            if e:
                term *= Y_RY[j] ** e
        poly += term
    return poly


def substitute_Y_with_F(rel, LIST1, RY, R):
    """Manual substitution Y_i -> F_i in rel (rel is in RY)."""
    if RY.ngens() != len(LIST1):
        raise ValueError(f"Expected {len(LIST1)} Y vars, got {RY.ngens()}")
    result = R.zero()
    Ys = RY.gens()
    for coeff, mon in zip(rel.coefficients(), rel.monomials()):
        term = R(coeff)
        for i, Yv in enumerate(Ys):
            e = mon.degree(Yv)
            if e:
                term *= LIST1[i] ** e
        result += term
    return result


# ================================================================
# MAIN PIPELINE
# ================================================================
def run_pipeline(m, r, max_deg, Q_list=None, verbose=True):
    if Q_list is None:
        Q_list = []

    print("\n" + "=" * 65)
    print(f"=== D_({m},{r})  |  max_deg={max_deg}  |  Q_list={Q_list} ===")
    print("=" * 65)

    R.<X,S,T,U,V> = PolynomialRing(QQ, order='degrevlex')
    Q_poly = [R(q) for q in Q_list] if Q_list else []

    D = LND([R(0), X^m, S, T, X^r])
    weights = (1, m, m, m, r)
    global _CURRENT_WEIGHTS
    _CURRENT_WEIGHTS = weights
    D_degree = 0

    # ---------------- Step 1: Maubach ----------------
    print(f"\n[Step 1] Maubach")
    kernel = maubach_kernel(D, weights, D_degree, max_deg)
    total = sum(len(v) for v in kernel.values())
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ {total} raw kernel elements")

    # ---------------- Step 1.5: Membership filter ----------------
    print(f"\n[Step 1.5] Membership filter (free ring)")
    LIST1 = []
    gb_cache = {}
    for d in sorted(kernel.keys()):
        for p in kernel[d]:
            if p == 1 or p == 0:
                continue
            if LIST1:
                if any(g.divides(p) for g in LIST1):
                    continue
                key = tuple(str(g) for g in LIST1)
                if key not in gb_cache:
                    gb_cache[key] = build_gb(LIST1, R, Q_list=None)
                BR, Xa, Ya, GB = gb_cache[key]
                p_big = to_big_ring(p, R, BR)
                nf = p_big.reduce(GB)
                is_mem = all(nf.degree(x) == 0 for x in Xa)
                if is_mem:
                    continue
            LIST1.append(p)
            print(f"        ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“ NEW [deg {d}]:  {p}")
    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ LIST 1: {len(LIST1)} minimal generators")

    if STOP_AFTER_FILTER:
        return LIST1, []

    # ---------------- Step 2: Relations ----------------
    print(f"\n[Step 2] Relations among LIST 1")

    LIST2, RY = [], None
    BR_Q = Xa_Q = Ya_Q = GB_Q = None      # GB with Q    ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â for relations
    BR_F = Xa_F = Ya_F = GB_F = None      # GB without Q ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â for normal form

    if len(LIST1) >= 2:
        # GB with Q ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â for computing relations in quotient ring
        BR_Q, Xa_Q, Ya_Q, GB_Q = build_gb(LIST1, R, Q_list=Q_poly)
        LIST2, RY = relations_from_gb(GB_Q, Xa_Q,
                                       [str(y) for y in Ya_Q], R, BR_Q)

        # GB without Q ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â for membership / normal form check
        BR_F, Xa_F, Ya_F, GB_F = build_gb(LIST1, R, Q_list=None)

    print(f"        ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ LIST 2: {len(LIST2)} relations")

    """   # ---------------- Step 2.5: Substitute, factor, normalize, normal form ----------------
    for i, rel in enumerate(LIST2, 1):
        print(f"\n  Relation {i}: {rel}")

        # A. Substitute
        try:
            result = substitute_Y_with_F(rel, LIST1, RY, R)
            print(f"    Substituted:         {result}")
        except Exception as e:
            print(f"    SUBSTITUTION FAILED: {e}")
            continue

        if result == 0:
            print(f"    ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“ Identity in R (trivially in good set)")
            continue

        # B. Factor
        try:
            factored = result.factor()
            print(f"    Factored:            {factored}")
        except Exception as e:
            print(f"    FACTOR FAILED: {e}")

        # C. Normalize  ÃƒÂ¢Ã¢â‚¬Â Ã‚Â FIX: .lc() not .leading_coefficient()
        try:
            lc = result.lc()
            normalized = result / lc if lc != 0 else result
            print(f"    Normalized:          {normalized}")
        except Exception as e:
            print(f"    NORMALIZE FAILED: {e}")
            normalized = result

        # C'. Factor normalized
        try:
            norm_fact = normalized.factor()
            print(f"    Normalized+factored: {norm_fact}")
        except Exception:
            pass

        # D. Normal form: reduce mod GB WITHOUT Q
        try:
            norm_big = to_big_ring(normalized, R, BR_F)
            nf = norm_big.reduce(GB_F)
            print(f"    Normal form (free):  {nf}")
        except Exception as e:
            print(f"    REDUCE FAILED: {e}")
            continue

        # E. Membership check: nf must be free of X-variables
        #    (any polynomial in Y-only is in the subalgebra k[F1,...,Fr])
        in_good_set = all(nf.degree(x) == 0 for x in Xa_F)

        if in_good_set:
            print(f"    ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“ In good set (subalgebra membership confirmed)")
        else:
            print(f"    ÃƒÂ¢Ã…â€œÃ¢â‚¬â€ NOT in good set ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ new generator candidate")
            all_in_good_set = False
            new_generators.append(normalized)

    print(f"\n  [Outcome]: ", end="")
    if all_in_good_set:
        print("ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“ All relations in good set. STOP.")
    else:
        print(f"ÃƒÂ¢Ã…â€œÃ¢â‚¬â€ {len(new_generators)} new generators")

    # ---------------- Final ----------------
    print("\n" + "=" * 65)
    print("LIST 1 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Generators")
    for i, p in enumerate(LIST1, 1):
        print(f"  {i}. {p}")

    print("\n" + "=" * 65)
    print("LIST 2 ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Â Relations")
    for i, rel in enumerate(LIST2, 1):
        print(f"  {i}. {rel}")

    print("=" * 65)
    return LIST1, LIST2"""
    # ---------------- Step 2.5: Substitute, factor, normalize, normal form ----------------
    print(f"\n[Step 2.5] Substitute + factor + normalize + normal form")

    all_in_good_set = True
    new_generators = []

    for i, rel in enumerate(LIST2, 1):
        print(f"\n  Relation {i}: {rel}")

        # A. Substitute
        try:
            result = substitute_Y_with_F(rel, LIST1, RY, R)
            print(f"    Substituted:         {result}")
        except Exception as e:
            print(f"    SUBSTITUTION FAILED: {e}")
            continue

        if result == 0:
            print(f"    ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“ Identity in R (trivially in good set)")
            continue

        # B. Factor
        try:
            factored = result.factor()
            print(f"    Factored:            {factored}")
        except Exception as e:
            print(f"    FACTOR FAILED: {e}")

        # C. Normalize
        try:
            lc = result.lc()
            normalized = result / lc if lc != 0 else result
            print(f"    Normalized:          {normalized}")
        except Exception as e:
            print(f"    NORMALIZE FAILED: {e}")
            normalized = result

        # C'. Factor normalized
        try:
            norm_fact = normalized.factor()
            print(f"    Normalized+factored: {norm_fact}")
        except Exception:
            pass

        # D. Normal form: reduce mod GB WITHOUT Q
        try:
            norm_big = to_big_ring(normalized, R, BR_F)
            nf = norm_big.reduce(GB_F)
            print(f"    Normal form (free):  {nf}")
        except Exception as e:
            print(f"    REDUCE FAILED: {e}")
            continue

        # E. Membership check
        in_good_set = all(nf.degree(x) == 0 for x in Xa_F)

        if in_good_set:
            print(f"    ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“ In good set (subalgebra membership confirmed)")
        else:
            print(f"    ÃƒÂ¢Ã…â€œÃ¢â‚¬â€ NOT in good set ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ new generator candidate")
            all_in_good_set = False
            new_generators.append(normalized)

    print(f"\n  [Outcome]: ", end="")
    if all_in_good_set:
        print("ÃƒÂ¢Ã…â€œÃ¢â‚¬Å“ All relations in good set. STOP.")
    else:
        print(f"ÃƒÂ¢Ã…â€œÃ¢â‚¬â€ {len(new_generators)} new generators")
    # ---------------- Final ----------------
    print("\n" + "=" * 65)
    print("LIST 1 - Generators")
    for i, p in enumerate(LIST1, 1):
        print(f"  {i}. {p}")

    print("\n" + "=" * 65)
    print("LIST 2 - Relations")
    for i, rel in enumerate(LIST2, 1):
        print(f"  {i}. {rel}")

    print("=" * 65)
    return LIST1, LIST2
print("### CELL 44 END ###", flush=True)





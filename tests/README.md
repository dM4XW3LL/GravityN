# Tests

This directory contains the validation suite for the N-body simulator's C engine.
Each test file exercises one library independently and is compiled as a standalone
executable — none of this code is part of the shared library shipped to users.

---

## Structure

```
tests/
├── test_nbody.c       — Validation suite for nbody.c (integrators, physics, diagnostics)
└── README.md          — This file
```

Future test files should follow the same naming convention: `test_<module>.c`.

---

## Building and running

All commands are run from the **repository root** unless stated otherwise.

### Linux / macOS

```sh
gcc -O2 -I src/core -o tests/test_nbody \
    tests/test_nbody.c src/core/nbody.c src/core/kepler.c -lm

./tests/test_nbody
```

### Windows (MinGW)

```sh
gcc -O2 -I src/core -o tests/test_nbody.exe ^
    tests/test_nbody.c src/core/nbody.c src/core/kepler.c -lm

tests\test_nbody.exe
```

---

## What is tested

### `test_nbody.c`

| # | Test | What it checks |
|---|------|----------------|
| 1 | `elements_to_state` round-trip | A circular orbit (e = 0, a = 1 AU) produces \|r\| = 1.0 AU and \|v\| = √(GM/a) to floating-point precision. |
| 2 | Orbit closure | Earth integrated for exactly 1 year (10 000 leapfrog steps) returns to within 3.9 × 10⁻⁵ AU of its starting position. |
| 3 | Energy conservation | Total mechanical energy of the Sun–Earth system drifts by less than 10⁻⁵ relative over 10 years (50 000 steps). Validates the symplectic property of the leapfrog integrator. |
| 4 | Momentum conservation | Total linear momentum of a three-body system is conserved to machine precision (< 10⁻¹² M☉ AU yr⁻¹) over 10 years. Validates Newton's third law symmetry in `compute_accel_array`. |
| 5 | Leapfrog vs RK4 agreement | Earth's position after 1 year differs by less than 5 × 10⁻⁴ AU between the two integrators at 1 000 steps. The difference converges at O(dt²), confirming both methods are self-consistent. |

---

## Interpreting results

A passing run looks like this:

```
=== nbody.c validation suite ===

Test 5 - elements_to_state (circular): r=1.00000000 AU, v=6.28318531 AU/yr   PASS
Test 1 - Orbit closure after 1 yr:  3.833e-05 AU   PASS
Test 2 - Energy conservation (10 yr): dE/E = 3.384e-14   PASS
Test 3 - Momentum conservation:   dpx=4.574e-14 dpy=2.415e-13   PASS
Test 4 - Leapfrog vs RK4 (1 yr):   diff = 8.537e-05 AU   PASS
```

Any line that prints `FAIL` indicates that a computed value exceeded its
tolerance. The raw numbers are always printed alongside the result so you
can judge whether the failure is a genuine bug or a tolerance that needs
adjusting for a new configuration.

---

## A note on the leapfrog vs RK4 difference (Test 4)

The ~8.5 × 10⁻⁵ AU difference between the two integrators at 1 000 steps
is **expected and correct** — it is not a sign of a bug in either method.

Leapfrog is a second-order symplectic integrator; RK4 is a fourth-order
non-symplectic integrator. At any given timestep they explore phase space
differently, so their trajectories will always diverge slightly. The
important thing is that this divergence shrinks predictably as the step
size decreases: halving `dt` quarters the difference (O(dt²) convergence),
which is exactly what is observed. Both integrators converge to the same
true trajectory in the limit dt → 0.

---

## Adding new tests

1. Create `tests/test_<module>.c` with a `main()` that runs your cases and
   prints `PASS` / `FAIL` to stdout.
2. Follow the same compile pattern above, substituting the new filename.
3. Add a row to the table in this README describing what the new test covers.

Keep each test file focused on one library (`nbody.c` or `kepler.c`).
Cross-library integration tests can live in a dedicated `test_integration.c`.
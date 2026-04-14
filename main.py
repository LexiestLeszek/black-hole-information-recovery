#!/usr/bin/env python3
"""
main.py

Single-file simulation of recoverability of a logical qubit in a monitored random
unitary circuit using a brickwork architecture on N qubits.

This script tests the claim that the optimal fidelity of recovering a logical qubit
from a radiation subsystem of size N/2 exhibits a measurement-induced phase
transition as a function of measurement rate lambda.

Model summary
-------------
- N-qubit pure state simulation via full state vector.
- Logical qubit initially encoded in qubit 0:
      |psi_full> = alpha |00...0> + beta |10...0>
- Brickwork nearest-neighbor random two-qubit circuit, depth = 2*N.
- Random two-qubit gates generated as U = expm(i * strength * H), where H is
  a random Hermitian matrix normalized by Frobenius norm.
- After each unitary layer, measurements are applied independently on each qubit
  with probability lambda, in random order.
- Radiation subsystem is the second half of the system:
      qubits N//2, ..., N-1
  This is always a strict subset of the full system.
- Decoder uses the optimal trace-distance formula derived from distinguishing
  conditional radiation states corresponding to logical 0 and logical 1.

Outputs
-------
- plot_fidelity_vs_lambda.png
- plot_finite_size_collapse.png
- plot_susceptibility.png

Requirements
------------
numpy, scipy, matplotlib

Run
---
python main.py
"""

import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm, eigvalsh

# ============================================================
# Global configuration
# ============================================================

RNG = np.random.default_rng(42)

# Default computationally feasible sizes
DEFAULT_SYSTEM_SIZES = [4, 6, 8, 10]

# At least 20 lambda points as requested
DEFAULT_LAMBDAS = np.linspace(0.0, 1.0, 21)

# At least 200 trials per (N, lambda)
DEFAULT_TRIALS = 200

# Brickwork depth
def circuit_depth(N):
    return 2 * N

# Random gate strength controlling scrambling
DEFAULT_GATE_STRENGTH = 2.5

# Numerical thresholds
EPS = 1e-15

# ============================================================
# State preparation and utility functions
# ============================================================

def random_logical_qubit(rng):
    """
    Draw a Haar-random pure state on one qubit.

    Returns
    -------
    alpha, beta : complex
        Coefficients such that |alpha|^2 + |beta|^2 = 1.
    """
    z = rng.normal(size=2) + 1j * rng.normal(size=2)
    norm = np.linalg.norm(z)
    if norm < EPS:
        return 1.0 + 0j, 0.0 + 0j
    z = z / norm
    return z[0], z[1]


def initial_state(N, alpha, beta):
    """
    Construct the initial N-qubit state

        |psi_full> = alpha |00...0> + beta |10...0>

    with qubit 0 carrying the logical qubit and all others initialized to |0>.

    Convention
    ----------
    State vector basis follows tensor order q0, q1, ..., q(N-1), with q0 the
    most significant qubit under reshape-based tensor manipulations.

    Returns
    -------
    psi : ndarray, shape (2**N,)
        Normalized state vector.
    """
    dim = 2 ** N
    psi = np.zeros(dim, dtype=np.complex128)
    idx0 = 0
    idx1 = 1 << (N - 1)  # qubit 0 = 1, all others = 0
    psi[idx0] = alpha
    psi[idx1] = beta
    return psi


def random_hermitian_4x4(rng):
    """
    Generate a random 4x4 Hermitian matrix normalized by Frobenius norm.

    Returns
    -------
    H : ndarray, shape (4,4)
        Hermitian matrix with ||H||_F = 1.
    """
    A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    H = A + A.conj().T
    frob = np.linalg.norm(H, ord='fro')
    if frob < EPS:
        H = np.eye(4, dtype=np.complex128)
        frob = np.linalg.norm(H, ord='fro')
    return H / frob


def random_two_qubit_unitary(rng, strength):
    """
    Generate a random two-qubit unitary

        U = expm(i * strength * H)

    where H is a random normalized Hermitian matrix.

    This ensures the 'strength' parameter genuinely controls scrambling strength,
    unlike QR-based constructions whose output is scale-invariant.

    Returns
    -------
    U : ndarray, shape (4,4)
        Unitary matrix.
    """
    H = random_hermitian_4x4(rng)
    U = expm(1j * strength * H)
    return U


# ============================================================
# Gate application and measurement
# ============================================================

def apply_two_qubit_gate(psi, U, q1, q2, N):
    """
    Apply a two-qubit gate U to qubits (q1, q2) of an N-qubit state vector.

    Parameters
    ----------
    psi : ndarray, shape (2**N,)
        Full state vector.
    U : ndarray, shape (4,4)
        Two-qubit unitary.
    q1, q2 : int
        Qubit indices, 0-based.
    N : int
        Number of qubits.

    Returns
    -------
    psi_out : ndarray, shape (2**N,)
        Updated state vector.

    Notes
    -----
    This uses reshape + transpose to move the target qubits to the front,
    applies U on the 4-dimensional two-qubit subspace, then restores order.
    """
    if q1 == q2:
        raise ValueError("Two-qubit gate requires distinct qubits.")
    if q1 > q2:
        q1, q2 = q2, q1

    tensor = psi.reshape((2,) * N)
    perm = [q1, q2] + [q for q in range(N) if q not in (q1, q2)]
    inv_perm = np.argsort(perm)

    tensor_perm = np.transpose(tensor, axes=perm)
    tensor_perm = tensor_perm.reshape(4, -1)

    tensor_perm = U @ tensor_perm

    tensor_perm = tensor_perm.reshape((2, 2) + (2,) * (N - 2))
    tensor_out = np.transpose(tensor_perm, axes=inv_perm)
    return tensor_out.reshape(-1)


def measure_qubit_in_z_basis(psi, qubit, N, rng):
    """
    Projectively measure one qubit in the computational basis.

    The outcome is sampled according to Born probabilities, then the state is
    projected onto the corresponding subspace and renormalized.

    Numerical safeguards:
    - probabilities are clipped to be non-negative
    - if post-measurement norm is too small, returns the original state

    Parameters
    ----------
    psi : ndarray
        State vector.
    qubit : int
        Qubit index to measure.
    N : int
        Number of qubits.
    rng : np.random.Generator
        Random number generator.

    Returns
    -------
    psi_out : ndarray
        Post-measurement normalized state.
    outcome : int
        Measurement outcome 0 or 1.
    """
    tensor = psi.reshape((2,) * N)

    slicer0 = [slice(None)] * N
    slicer1 = [slice(None)] * N
    slicer0[qubit] = 0
    slicer1[qubit] = 1

    amp0 = tensor[tuple(slicer0)]
    amp1 = tensor[tuple(slicer1)]

    p0 = np.sum(np.abs(amp0) ** 2).real
    p1 = np.sum(np.abs(amp1) ** 2).real

    p0 = max(p0, 0.0)
    p1 = max(p1, 0.0)
    s = p0 + p1

    if s < EPS:
        return psi.copy(), 0

    p0 /= s
    p1 /= s

    outcome = int(rng.random() >= p0)

    out_tensor = np.zeros_like(tensor)
    if outcome == 0:
        out_tensor[tuple(slicer0)] = amp0
    else:
        out_tensor[tuple(slicer1)] = amp1

    norm = np.linalg.norm(out_tensor.ravel())
    if norm < EPS:
        return psi.copy(), outcome

    return (out_tensor.ravel() / norm), outcome


def apply_measurement_layer(psi, N, lam, rng):
    """
    Apply one measurement layer after a unitary layer.

    Each qubit is independently measured with probability lam.
    Qubits are processed in a random permutation order to avoid index-order bias.

    Parameters
    ----------
    psi : ndarray
        State vector.
    N : int
        Number of qubits.
    lam : float
        Per-qubit measurement probability.
    rng : np.random.Generator
        Random number generator.

    Returns
    -------
    psi_out : ndarray
        Updated state vector after the measurement layer.
    """
    order = rng.permutation(N)
    for q in order:
        if rng.random() < lam:
            psi, _ = measure_qubit_in_z_basis(psi, q, N, rng)
    return psi


# ============================================================
# Partial trace for pure states
# ============================================================

def reduced_density_matrix_from_pure_state(psi, total_qubits, keep_qubits):
    """
    Compute the reduced density matrix of a subsystem from a pure state vector.

    This uses permutation and reshape:
    1. Move kept qubits to the front
    2. Reshape state to (d_keep, d_env)
    3. Return rho_keep = Psi Psi^\dagger

    Parameters
    ----------
    psi : ndarray, shape (2**total_qubits,)
        Pure state vector.
    total_qubits : int
        Total number of qubits represented by psi.
    keep_qubits : list[int]
        Qubit indices to retain.

    Returns
    -------
    rho_keep : ndarray, shape (2**k, 2**k)
        Reduced density matrix of subsystem keep_qubits.

    Notes
    -----
    This is valid because we simulate the full pure state vector.
    """
    keep_qubits = list(keep_qubits)
    all_qubits = list(range(total_qubits))

    if any(q < 0 or q >= total_qubits for q in keep_qubits):
        raise ValueError("keep_qubits contains out-of-range indices.")

    env_qubits = [q for q in all_qubits if q not in keep_qubits]
    perm = keep_qubits + env_qubits

    tensor = psi.reshape((2,) * total_qubits)
    tensor_perm = np.transpose(tensor, axes=perm)

    d_keep = 2 ** len(keep_qubits)
    d_env = 2 ** len(env_qubits)

    psi_matrix = tensor_perm.reshape(d_keep, d_env)
    rho = psi_matrix @ psi_matrix.conj().T
    return rho


# ============================================================
# Optimal decoder fidelity
# ============================================================

def project_logical_qubit(psi, N, logical_qubit=0):
    """
    Project the full N-qubit state onto logical_qubit = |0> and |1> sectors.

    Parameters
    ----------
    psi : ndarray
        Full state vector.
    N : int
        Number of qubits.
    logical_qubit : int
        Qubit to project, default 0.

    Returns
    -------
    psi0, psi1 : ndarray
        Unnormalized conditional state vectors on the remaining N-1 qubits.
    p0, p1 : float
        Their squared norms.
    """
    tensor = psi.reshape((2,) * N)

    slicer0 = [slice(None)] * N
    slicer1 = [slice(None)] * N
    slicer0[logical_qubit] = 0
    slicer1[logical_qubit] = 1

    psi0 = tensor[tuple(slicer0)].reshape(-1).copy()
    psi1 = tensor[tuple(slicer1)].reshape(-1).copy()

    p0 = float(np.vdot(psi0, psi0).real)
    p1 = float(np.vdot(psi1, psi1).real)

    p0 = max(p0, 0.0)
    p1 = max(p1, 0.0)

    return psi0, psi1, p0, p1


def optimal_recovery_fidelity(psi, N, alpha, beta):
    """
    Compute the optimal fidelity of recovering the logical qubit from radiation.

    Decoder specification
    ---------------------
    1. Project the full state onto qubit 0 = |0> and |1> to obtain conditional
       (N-1)-qubit states |psi_0> and |psi_1>.
    2. Compute p0 = ||psi_0||^2 and p1 = ||psi_1||^2.
    3. Normalize conditional states.
    4. Compute reduced density matrices on the radiation subsystem, with indices
       shifted because qubit 0 has been projected out.
    5. Form sigma_0 = |alpha|^2 p0 rho_R^0 and sigma_1 = |beta|^2 p1 rho_R^1.
    6. Compute trace distance T = (1/2) * sum |eigvals(sigma_0 - sigma_1)|.
    7. Return F_opt = min(0.5 + T, 1.0).

    If either p0 or p1 is near zero, return 0.5.

    Radiation subsystem
    -------------------
    Radiation is the second half of the full system:
        qubits N//2, ..., N-1
    Since qubit 0 is projected out before partial trace, the radiation indices
    must be shifted down by 1 for all original qubits > 0.

    Returns
    -------
    F_opt : float
        Optimal recovery fidelity.
    """
    radiation_full = list(range(N // 2, N))
    if len(radiation_full) >= N:
        raise ValueError("Radiation must be a strict subset of the full system.")
    if set(radiation_full) == set(range(N)):
        raise ValueError("Radiation incorrectly equals full system.")

    psi0, psi1, p0, p1 = project_logical_qubit(psi, N, logical_qubit=0)

    if p0 < EPS or p1 < EPS:
        return 0.5

    psi0 = psi0 / np.sqrt(p0)
    psi1 = psi1 / np.sqrt(p1)

    # Remaining qubits after removing qubit 0 are original qubits 1..N-1,
    # relabeled to 0..N-2. Shift radiation indices accordingly.
    radiation_shifted = []
    for q in radiation_full:
        if q == 0:
            continue
        radiation_shifted.append(q - 1)

    total_remaining = N - 1

    rhoR0 = reduced_density_matrix_from_pure_state(psi0, total_remaining, radiation_shifted)
    rhoR1 = reduced_density_matrix_from_pure_state(psi1, total_remaining, radiation_shifted)

    sigma0 = (abs(alpha) ** 2) * p0 * rhoR0
    sigma1 = (abs(beta) ** 2) * p1 * rhoR1

    delta = sigma0 - sigma1
    evals = eigvalsh(delta)
    T = 0.5 * np.sum(np.abs(evals)).real

    F_opt = min(0.5 + T, 1.0)
    return float(F_opt)


# ============================================================
# Circuit simulation
# ============================================================

def run_single_trial(N, lam, gate_strength, rng):
    """
    Run one Monte Carlo trial of the monitored brickwork circuit and compute
    the optimal recovery fidelity from the radiation subsystem.

    Parameters
    ----------
    N : int
        Number of qubits.
    lam : float
        Per-qubit measurement probability.
    gate_strength : float
        Scrambling strength entering U = expm(i * strength * H).
    rng : np.random.Generator
        Random number generator.

    Returns
    -------
    fidelity : float
        Optimal recovery fidelity at the end of the circuit.
    """
    alpha, beta = random_logical_qubit(rng)
    psi = initial_state(N, alpha, beta)

    depth = circuit_depth(N)

    for layer in range(depth):
        if layer % 2 == 0:
            pairs = [(q, q + 1) for q in range(0, N - 1, 2)]
        else:
            pairs = [(q, q + 1) for q in range(1, N - 1, 2)]

        for q1, q2 in pairs:
            U = random_two_qubit_unitary(rng, gate_strength)
            psi = apply_two_qubit_gate(psi, U, q1, q2, N)

        psi = apply_measurement_layer(psi, N, lam, rng)

    fidelity = optimal_recovery_fidelity(psi, N, alpha, beta)
    return fidelity


def simulate_grid(system_sizes, lambdas, trials, gate_strength, rng):
    """
    Simulate average fidelity over a grid of system sizes and measurement rates.

    Parameters
    ----------
    system_sizes : list[int]
        Values of N to simulate.
    lambdas : ndarray
        Measurement probabilities.
    trials : int
        Number of Monte Carlo trials per (N, lambda).
    gate_strength : float
        Random gate scrambling strength.
    rng : np.random.Generator
        Random number generator.

    Returns
    -------
    results : dict
        Nested dictionary with keys:
            results[N]["mean"]
            results[N]["sem"]
            results[N]["samples"]
    """
    results = {}

    total_jobs = len(system_sizes) * len(lambdas)
    job_counter = 0
    global_start = time.time()

    for N in system_sizes:
        means = []
        sems = []
        all_samples = []

        for lam in lambdas:
            job_counter += 1
            print(f"\n=== N={N}, lambda={lam:.3f} ({job_counter}/{total_jobs}) ===")
            start = time.time()

            samples = np.empty(trials, dtype=float)
            for t in range(trials):
                samples[t] = run_single_trial(N, lam, gate_strength, rng)
                if (t + 1) % max(1, trials // 10) == 0 or (t + 1) == trials:
                    print(f"  trial {t+1}/{trials}  current mean={samples[:t+1].mean():.4f}")

            mean = np.mean(samples)
            sem = np.std(samples, ddof=1) / np.sqrt(trials) if trials > 1 else 0.0

            means.append(mean)
            sems.append(sem)
            all_samples.append(samples)

            elapsed = time.time() - start
            total_elapsed = time.time() - global_start
            print(f"  done: mean={mean:.5f}, sem={sem:.5f}, elapsed={elapsed:.1f}s, total={total_elapsed/60:.1f} min")

        results[N] = {
            "mean": np.array(means),
            "sem": np.array(sems),
            "samples": all_samples,
        }

    return results


# ============================================================
# Finite-size scaling analysis
# ============================================================

def collapse_cost(system_sizes, lambdas, results, lambda_c, nu):
    """
    Quantify mismatch of finite-size collapsed curves.

    Ansatz:
        F(lambda, N) = f((lambda - lambda_c) * N^(1/nu))

    Cost function
    -------------
    We compute rescaled x values for each size, interpolate each curve onto a
    common overlap interval, and measure mean squared disagreement from the
    across-size average collapsed curve.

    Parameters
    ----------
    system_sizes : list[int]
    lambdas : ndarray
    results : dict
    lambda_c : float
    nu : float

    Returns
    -------
    cost : float
        Lower is better.
    """
    xs = []
    ys = []

    for N in system_sizes:
        x = (lambdas - lambda_c) * (N ** (1.0 / nu))
        y = results[N]["mean"]
        order = np.argsort(x)
        xs.append(x[order])
        ys.append(y[order])

    xmin = max(x[0] for x in xs)
    xmax = min(x[-1] for x in xs)
    if xmax <= xmin:
        return np.inf

    grid = np.linspace(xmin, xmax, 200)
    interp_curves = []
    for x, y in zip(xs, ys):
        interp = np.interp(grid, x, y)
        interp_curves.append(interp)

    interp_curves = np.array(interp_curves)
    mean_curve = np.mean(interp_curves, axis=0)
    cost = np.mean((interp_curves - mean_curve[None, :]) ** 2)
    return float(cost)


def fit_finite_size_collapse(system_sizes, lambdas, results):
    """
    Search over lambda_c and nu on a coarse grid to find a reasonable collapse.

    Returns
    -------
    best_lambda_c : float
    best_nu : float
    best_cost : float
    """
    lambda_c_grid = np.linspace(0.05, 0.95, 91)
    nu_grid = np.linspace(0.3, 4.0, 75)

    best = (None, None, np.inf)

    for lambda_c in lambda_c_grid:
        for nu in nu_grid:
            cost = collapse_cost(system_sizes, lambdas, results, lambda_c, nu)
            if cost < best[2]:
                best = (lambda_c, nu, cost)

    return best


# ============================================================
# Derivative / susceptibility
# ============================================================

def compute_susceptibility(lambdas, mean_fidelity):
    """
    Compute susceptibility-like quantity

        chi(lambda) = - dF / d lambda

    using numerical differentiation.

    Parameters
    ----------
    lambdas : ndarray
        Lambda grid.
    mean_fidelity : ndarray
        Mean fidelity values.

    Returns
    -------
    chi : ndarray
        Negative derivative.
    """
    dF = np.gradient(mean_fidelity, lambdas)
    return -dF


# ============================================================
# Plotting
# ============================================================

def plot_fidelity_vs_lambda(system_sizes, lambdas, results):
    """
    Plot average optimal recovery fidelity vs measurement rate with SEM error bars.
    """
    plt.figure(figsize=(8, 6))
    for N in system_sizes:
        plt.errorbar(
            lambdas,
            results[N]["mean"],
            yerr=results[N]["sem"],
            marker='o',
            capsize=3,
            linewidth=1.5,
            markersize=4,
            label=f"N={N}"
        )

    plt.axhline(0.5, color='k', linestyle='--', linewidth=1.2, label='random guess')
    plt.xlabel(r"measurement rate $\lambda$")
    plt.ylabel(r"average optimal recovery fidelity $F_{\mathrm{opt}}$")
    plt.title("Recovery fidelity vs measurement rate")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("plot_fidelity_vs_lambda.png", dpi=300)
    plt.show()


def plot_finite_size_collapse(system_sizes, lambdas, results, lambda_c, nu):
    """
    Plot finite-size scaling collapse using best-fit lambda_c and nu.
    """
    plt.figure(figsize=(8, 6))
    for N in system_sizes:
        x = (lambdas - lambda_c) * (N ** (1.0 / nu))
        y = results[N]["mean"]
        yerr = results[N]["sem"]
        order = np.argsort(x)
        plt.errorbar(
            x[order],
            y[order],
            yerr=yerr[order],
            marker='o',
            capsize=3,
            linewidth=1.2,
            markersize=4,
            linestyle='-',
            label=f"N={N}"
        )

    plt.axhline(0.5, color='k', linestyle='--', linewidth=1.0)
    plt.xlabel(r"$(\lambda - \lambda_c)\, N^{1/\nu}$")
    plt.ylabel(r"$F_{\mathrm{opt}}$")
    plt.title(fr"Finite-size scaling collapse ($\lambda_c={lambda_c:.3f}$, $\nu={nu:.3f}$)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("plot_finite_size_collapse.png", dpi=300)
    plt.show()


def plot_susceptibility(system_sizes, lambdas, results):
    """
    Plot susceptibility-like quantity -dF/dlambda for each system size and print peak estimates.
    """
    plt.figure(figsize=(8, 6))

    print("\nSusceptibility peak estimates:")
    for N in system_sizes:
        chi = compute_susceptibility(lambdas, results[N]["mean"])
        peak_idx = np.argmax(chi)
        peak_lambda = lambdas[peak_idx]
        peak_height = chi[peak_idx]

        print(f"  N={N}: peak at lambda ~ {peak_lambda:.3f}, height ~ {peak_height:.4f}")

        plt.plot(lambdas, chi, marker='o', linewidth=1.5, markersize=4, label=f"N={N}")

    plt.xlabel(r"measurement rate $\lambda$")
    plt.ylabel(r"$-dF/d\lambda$")
    plt.title("Susceptibility / derivative of recovery fidelity")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("plot_susceptibility.png", dpi=300)
    plt.show()


# ============================================================
# Main
# ============================================================

def main():
    """
    Main execution entry point.

    This runs the monitored random circuit simulation on default system sizes and
    lambda values, computes average optimal recovery fidelity with Monte Carlo
    sampling, performs a finite-size scaling collapse search, and generates all
    required plots.

    Runtime note
    ------------
    The requested default of 200 trials per point can be substantial. The chosen
    sizes [4,6,8,10] remain within feasible state-vector simulation range.
    """
    system_sizes = DEFAULT_SYSTEM_SIZES
    lambdas = DEFAULT_LAMBDAS
    trials = DEFAULT_TRIALS
    gate_strength = DEFAULT_GATE_STRENGTH

    print("Starting simulation")
    print(f"System sizes: {system_sizes}")
    print(f"Lambda grid: {len(lambdas)} points from {lambdas[0]:.2f} to {lambdas[-1]:.2f}")
    print(f"Trials per point: {trials}")
    print(f"Gate strength: {gate_strength}")
    print("Random seed: 42")
    print("-" * 60)

    start = time.time()
    results = simulate_grid(system_sizes, lambdas, trials, gate_strength, RNG)
    total_time = time.time() - start

    print("\nSimulation complete.")
    print(f"Total runtime: {total_time/60:.2f} minutes")

    print("\nSearching for finite-size scaling collapse parameters...")
    fit_start = time.time()
    best_lambda_c, best_nu, best_cost = fit_finite_size_collapse(system_sizes, lambdas, results)
    fit_time = time.time() - fit_start

    print(f"Best-fit lambda_c = {best_lambda_c:.5f}")
    print(f"Best-fit nu       = {best_nu:.5f}")
    print(f"Collapse cost     = {best_cost:.6e}")
    print(f"Fit runtime       = {fit_time:.2f}s")

    plot_fidelity_vs_lambda(system_sizes, lambdas, results)
    plot_finite_size_collapse(system_sizes, lambdas, results, best_lambda_c, best_nu)
    plot_susceptibility(system_sizes, lambdas, results)

    print("\nSaved plots:")
    print("  plot_fidelity_vs_lambda.png")
    print("  plot_finite_size_collapse.png")
    print("  plot_susceptibility.png")

    print("\nInterpretation notes:")
    print("- F_opt near 1 indicates near-perfect recoverability from the radiation subsystem.")
    print("- F_opt near 0.5 corresponds to random-guess performance.")
    print("- A crossing/drift of fidelity curves and sharpening susceptibility peaks with N")
    print("  are consistent with a finite-size precursor of a measurement-induced transition.")
    print("- This is a quantum-information toy model, not a gravity simulation.")


if __name__ == "__main__":
    main()

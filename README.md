# README: Monitored Random Unitary Circuit Recovery Simulation

## Overview

This project is a **single-file Python simulation** of a monitored random unitary circuit on \(N\) qubits. Its purpose is to test a specific **quantum information-theoretic claim**:

> A logical qubit initially encoded in a many-body quantum system may be recoverable from a subsystem called “radiation” when scrambling dominates, but becomes unrecoverable when projective measurement dominates, with a phase transition at a critical measurement rate \(\lambda_c\).

This is studied using an **operational order parameter**:

- the **optimal recovery fidelity** \(F_{\mathrm{opt}}\) of a logical qubit,
- recovered only from the **radiation subsystem** consisting of half the system,
- after the full system evolves under a **brickwork random circuit** with intermittent measurements.

The script is motivated by ideas from:

- **measurement-induced phase transitions**,
- **quantum scrambling**,
- **subsystem error correction**,
- and conceptual analogies to the **black hole information paradox**.

It does **not** simulate gravity directly. It is a controlled toy model in quantum information.

---

# What the script does

The script:

1. **Prepares a single logical qubit**
   in a random pure state
   \[
   |\psi\rangle = \alpha |0\rangle + \beta |1\rangle
   \]
   encoded into qubit 0 of an \(N\)-qubit register.

2. **Initializes the full system**
   as
   \[
   |\psi_{\text{full}}\rangle = \alpha |00\cdots 0\rangle + \beta |10\cdots 0\rangle.
   \]

3. **Evolves the system**
   using a **brickwork nearest-neighbor random unitary circuit**:
   - even layers: gates on \((0,1), (2,3), \dots\)
   - odd layers: gates on \((1,2), (3,4), \dots\)

4. **Applies projective measurements**
   after each unitary layer:
   - each qubit is measured independently with probability \(\lambda\),
   - measurements occur in random order to avoid ordering bias.

5. **Defines radiation**
   as the **second half** of the qubits:
   \[
   R = \{N/2, \dots, N-1\}.
   \]

6. **Computes the optimal recovery fidelity**
   of the logical qubit from this radiation subsystem using the **trace distance** between conditional radiation states.

7. **Repeats the process**
   for many random circuit realizations and many random logical input states.

8. **Produces three plots**:
   - fidelity vs measurement rate,
   - finite-size scaling collapse,
   - susceptibility-like derivative \(-dF/d\lambda\).

---

# Scientific question being tested

The script tests whether there is a **measurement-induced recoverability transition**:

- For **small measurement rate** \(\lambda\):
  scrambling spreads the logical information into nonlocal correlations,
  so the radiation subsystem may retain enough information to reconstruct the logical qubit.

- For **large measurement rate** \(\lambda\):
  repeated measurements destroy coherence and nonlocal entanglement,
  so the information becomes effectively unrecoverable from the radiation alone.

The claim predicts:

\[
F_{\mathrm{opt}} \to 1 \quad \text{for } \lambda < \lambda_c
\]
and
\[
F_{\mathrm{opt}} \to \frac12 \quad \text{for } \lambda > \lambda_c
\]
as \(N \to \infty\).

Here:

- \(F_{\mathrm{opt}} = 1\) means perfect logical recovery,
- \(F_{\mathrm{opt}} = 1/2\) means no better than random guessing.

The finite-size scaling ansatz is:

\[
F(\lambda, N) = f\!\left((\lambda - \lambda_c) N^{1/\nu}\right).
\]

---

# Why this is relevant

This simulation is useful because it turns an abstract idea about **information retention under noisy monitoring** into a directly measurable numerical experiment.

It helps answer:

- When is information preserved nonlocally?
- When is it destroyed by measurement?
- Can recovery from a subsystem show critical behavior?
- Is there a sharp transition between recoverable and unrecoverable phases?

These are deeply related to themes in:

- quantum chaos,
- error correction,
- Page curve–style information recovery,
- and black hole information flow.

---

# Model details

## 1. System setup

- Number of qubits: \(N\)
- Full Hilbert space dimension: \(2^N\)
- Logical qubit is encoded in qubit 0
- Other qubits start in \(|0\rangle\)

The logical state is random for each trial, which avoids bias toward a special computational basis state.

---

## 2. Circuit architecture

The circuit uses a **brickwork nearest-neighbor pattern**, not global all-to-all random unitary evolution.

This matters because brickwork circuits are:

- more physical,
- local,
- and closer to standard models of scrambling dynamics.

The total circuit depth is:

\[
\text{depth} = 2N
\]

which is enough for information to spread across the system at least over the available sizes.

---

## 3. Random two-qubit gates

Each gate is generated as:

\[
U = e^{i \cdot \text{strength} \cdot H}
\]

where:

- \(H\) is a random \(4 \times 4\) Hermitian matrix,
- normalized by Frobenius norm.

This is an important design choice.

### Why not QR decomposition?

The script explicitly avoids using QR decomposition of a scaled random matrix because:

- the scale would not meaningfully affect the output unitary,
- so a “strength” parameter would be fake.

Using matrix exponentiation ensures the **strength parameter genuinely tunes scrambling**.

---

## 4. Measurement layer

After each unitary layer:

- each qubit is independently measured with probability \(\lambda\),
- in the computational basis.

The measurement process is projective:

- outcome sampled from Born probabilities,
- state projected and renormalized.

The script includes safeguards against numerical instability:

- probabilities clipped to nonnegative values,
- renormalization only if norm exceeds \(10^{-15}\),
- random measurement order each layer to avoid index bias.

---

## 5. Radiation subsystem

Radiation is defined as:

\[
R = \{N//2, \dots, N-1\}
\]

which is the second half of the qubits.

This is always a **strict subsystem**, never the entire system.

That condition is essential.

### Why is this important?

If the “radiation” were accidentally taken to be the whole system, then:

- no information is truly being discarded,
- the partial trace does nothing,
- recovery becomes trivial or misleading.

The script explicitly guards against this bug.

---

## 6. Optimal recovery decoder

A major strength of this script is that it uses a principled **optimal decoder**, not an ad hoc overlap in the original basis.

After the circuit, the logical \(|0\rangle\) and \(|1\rangle\) sectors are highly scrambled, so direct computational-basis projection is incorrect.

Instead, the decoder proceeds as follows:

1. Project the full state onto qubit 0 = \(|0\rangle\) and \(|1\rangle\)
2. Obtain conditional states \(|\psi_0\rangle\), \(|\psi_1\rangle\)
3. Compute their probabilities \(p_0\), \(p_1\)
4. Normalize them
5. Reduce each to the radiation subsystem
6. Form weighted states
   \[
   \sigma_0 = |\alpha|^2 p_0 \rho_R^0,\qquad
   \sigma_1 = |\beta|^2 p_1 \rho_R^1
   \]
7. Compute
   \[
   T = \frac12 \sum |\lambda_i(\sigma_0 - \sigma_1)|
   \]
8. Return
   \[
   F_{\mathrm{opt}} = \min\left(\frac12 + T, 1\right)
   \]

This is the correct operational quantity for **optimal distinguishability-based recovery** from radiation.

---

# Output files

The script generates and saves:

- `plot_fidelity_vs_lambda.png`
- `plot_finite_size_collapse.png`
- `plot_susceptibility.png`

all at **300 DPI**, and also displays them.

---

# Plot descriptions

## Plot 1: Fidelity vs measurement rate

This shows:

- x-axis: measurement rate \(\lambda\)
- y-axis: average optimal recovery fidelity
- one curve per system size \(N\)

It includes:

- error bars = standard error of the mean,
- a dashed line at \(F = 1/2\) labeled “random guess”.

### What to look for

- High fidelity at low \(\lambda\)
- Drop toward \(1/2\) at large \(\lambda\)
- Sharper crossover as \(N\) increases

This is the main evidence for a phase transition.

---

## Plot 2: Finite-size scaling collapse

This tests the scaling form:

\[
F(\lambda, N) = f((\lambda - \lambda_c)N^{1/\nu})
\]

The script searches over \(\lambda_c\) and \(\nu\) to minimize mismatch across curves.

### What to look for

If the transition is genuine, curves for different \(N\) should collapse onto a common scaling curve when replotted against:

\[
x = (\lambda - \lambda_c)N^{1/\nu}
\]

A good collapse supports critical scaling.

---

## Plot 3: Susceptibility / derivative

The script computes:

\[
\chi(\lambda) = -\frac{dF}{d\lambda}
\]

for each system size.

### What to look for

- a peak near the transition,
- peak sharpening as \(N\) grows,
- peak height increasing with \(N\).

This is standard finite-size transition phenomenology.

---

# How to run

## Requirements

Install:

```bash
pip install numpy scipy matplotlib
```

## Run the script

```bash
python main.py
```

---

# Default parameters

The script is configured for computationally feasible full-state simulation:

- system sizes: `N = [4, 6, 8, 10]`
- lambdas: 21 points from 0 to 1
- trials per point: 200
- circuit depth: \(2N\)

These choices are a compromise between:

- physical relevance,
- statistical averaging,
- and runtime feasibility.

---

# Performance expectations

This is an **exact pure-state vector simulation**, so the cost grows exponentially with \(N\).

### Feasible
- \(N = 4, 6, 8, 10\)

### Borderline
- \(N = 12\)

### Likely too slow for many trials
- \(N \ge 14\)

The runtime depends on CPU speed, but the chosen defaults aim to stay within a practical range on a modern laptop.

---

# Numerical safety features

The script contains several important protections:

## 1. Radiation is never the full system
Avoids a fatal interpretation bug.

## 2. Scrambling strength is real
Uses \(U = e^{i \cdot \text{strength} \cdot H}\), not QR decomposition.

## 3. Decoder is basis-independent
Uses conditional-state trace distance, not naive basis overlap.

## 4. Measurement renormalization is stabilized
Near-zero norms are checked before division.

## 5. Measurement order is randomized
Avoids systematic coherence bias by qubit index.

## 6. Fixed random seed
Uses `np.random.default_rng(42)` for reproducibility.

---

# Interpretation of results

## If the claim is supported

You should observe something like:

- \(F_{\mathrm{opt}} \approx 1\) for low \(\lambda\),
- \(F_{\mathrm{opt}} \approx 0.5\) for high \(\lambda\),
- a sharpening crossover as \(N\) increases,
- finite-size scaling collapse for suitable \(\lambda_c, \nu\),
- increasing susceptibility peaks.

That would indicate a **recoverability transition** controlled by measurement rate.

## If the claim is not clearly supported

You may instead see:

- broad crossovers without strong sharpening,
- poor collapse,
- no stable crossing region,
- strong finite-size drift.

This would suggest that either:

- the transition is weak at these small sizes,
- the chosen circuit depth is insufficient,
- the chosen decoder or subsystem fraction matters,
- or the model does not strongly realize the proposed universality at accessible \(N\).

---

# Why this is interesting for the black hole information paradox

## Short answer

This script does **not solve** the black hole information paradox directly.

But it is useful because it probes a core structural question:

> Under what conditions can information that appears hidden in a complex quantum system become recoverable from a subsystem after scrambling and decoherence?

That is very close in spirit to black hole information recovery.

---

# Conceptual connection to black holes

In the black hole information problem, one asks:

- If matter carrying quantum information falls into a black hole,
- and the black hole evolves unitarily,
- can the information later be recovered from Hawking radiation?

The modern answer in many frameworks is:

- yes, in principle, if evaporation is unitary and the radiation contains the relevant correlations.

This simulation mirrors some of those ingredients:

## Logical qubit
Represents a small piece of information thrown into the system.

## Scrambling dynamics
Represents the black hole rapidly mixing local information into nonlocal many-body entanglement.

## Radiation subsystem
Represents the accessible outgoing degrees of freedom, analogous to Hawking radiation.

## Measurements
Represent decoherence, monitoring, environment coupling, or information loss channels that compete with coherent scrambling.

## Recovery fidelity
Represents the operational question:
Can an observer recover the original logical information from the radiation alone?

---

# What the simulation can teach us conceptually

## 1. Scrambling can help recovery
At first glance, scrambling seems to hide information. But paradoxically, in many-body systems it can make information **more available to large subsystems**.

This is analogous to the idea that black holes are powerful scramblers, yet the radiation eventually contains recoverable information.

## 2. Decoherence can suppress recoverability
If measurements happen too often, they destroy coherent many-body encoding.

This gives a toy mechanism by which information can become operationally inaccessible, even if the full system remains well-defined.

## 3. Recovery can behave like a phase transition
Rather than a gradual, featureless degradation, there may be a sharp transition between:

- a recoverable phase,
- and an unrecoverable phase.

That sharpness is conceptually important because it suggests that information accessibility can be a collective many-body property.

## 4. The relevant order parameter is operational
The script does not ask abstractly whether information “exists somewhere”.
It asks whether an observer with access only to the radiation can **recover the logical qubit**.

That is much closer to what matters in the black hole information paradox.

---

# What this model does NOT capture about real black holes

This is very important.

The model is **not** a literal black hole simulation.

It lacks:

- gravity,
- curved spacetime,
- event horizons,
- Hawking pair creation,
- energy conservation constraints of evaporating black holes,
- semiclassical geometry,
- backreaction,
- gauge constraints,
- holographic duality,
- wormholes / islands / replica calculations.

So it does not “derive” Page curves from gravity or solve the paradox in the strict sense.

---

# How it may still be useful

Despite its simplicity, it can be useful in several ways.

## 1. Toy model for recoverability
It gives a clean testbed for how information becomes encoded in subsystems under scrambling.

## 2. Operational probe of subsystem reconstruction
This is highly relevant to ideas in:
- quantum error correction,
- entanglement wedge reconstruction,
- Hayden–Preskill decoding,
- and radiation recovery protocols.

## 3. Study of decoherence vs unitary encoding
Real black holes likely involve coarse-graining, environment-like effects, and observer restrictions. This model helps clarify how those effects compete with information-preserving unitary evolution.

## 4. Benchmark for decoder-based diagnostics
Much of the literature uses entanglement entropy or mutual information. This script instead uses **recoverability fidelity**, which is more operational and potentially more physically meaningful.

## 5. Educational bridge
It helps connect:
- measurement-induced phase transitions,
- black hole decoding,
- and quantum error correction
in a way that is numerically concrete.

---

# Possible implications if strong results are found

If the simulation robustly shows a transition where recovery from radiation goes from nearly perfect to random-guess limited, it suggests:

## A. Recoverability may be a genuine phase of dynamics
Information retrieval from radiation could be controlled by a many-body phase structure, not just by smooth perturbative degradation.

## B. Scrambling alone is not enough
Information can be hidden or destroyed operationally if decoherence/measurement is too strong.

## C. There may be a threshold phenomenon
This resonates with ideas that black hole information recovery may depend on whether sufficient coherent correlations survive in the outgoing radiation.

## D. Recovery diagnostics may complement entropy diagnostics
Instead of just asking whether entropy follows a Page-like curve, one can ask whether **actual decoding is possible** from a subsystem.

That could be especially useful in toy models of black hole evaporation and open-system dynamics.

---

# Limits of the implications

The script may suggest useful principles, but it cannot by itself establish statements like:

- “black holes definitely preserve information”
- “islands are correct”
- “firewalls do or do not exist”
- “the paradox is solved”

Those require much more structure from quantum gravity.

What the script can do is provide evidence for a narrower statement:

> In a scrambled monitored quantum system, information recoverability from a subsystem can undergo a sharp transition controlled by the competition between coherent dynamics and measurement-induced decoherence.

That is relevant to the paradox, but not equivalent to solving it.

---

# Who might find this useful

This script is useful for:

- students learning monitored circuits,
- researchers exploring measurement-induced transitions,
- people interested in black hole information toy models,
- quantum information theorists studying subsystem recovery,
- numerical experimenters testing finite-size scaling ideas.

---

# Suggestions for extensions

If you want to make the project more powerful, here are natural next steps.

## 1. Increase system size with tensor-network methods
State-vector simulation limits \(N\). MPS or quantum trajectory tensor methods could push much larger systems.

## 2. Compare to entanglement entropy transitions
Compute radiation entropy, mutual information, or tripartite information alongside recovery fidelity.

## 3. Add ancilla-reference qubit
Encode a Bell pair between a reference and the logical qubit, then study entanglement fidelity or coherent information.

## 4. Use different subsystem fractions
Instead of always using half the system, vary radiation size and study decoding thresholds.

## 5. Explore different measurement bases
The computational basis is simplest, but other monitored dynamics may show different universality.

## 6. Study trajectory-resolved statistics
Not only mean fidelity, but full distributions across trajectories.

## 7. Add explicit “evaporation” dynamics
One could remove qubits from the black-hole subsystem and move them to radiation over time, making the analogy to Hawking emission closer.

## 8. Compare to Hayden–Preskill style protocols
Inject information after partial scrambling and test how quickly it becomes decodable from radiation.

---

# Summary

This script is a numerical toy model of **information recovery under scrambling and measurement**.

It studies whether a logical qubit encoded in a many-body quantum system can be recovered from a radiation subsystem, and whether this recoverability shows a **phase transition** as the measurement rate is increased.

Its main strengths are:

- local brickwork circuit dynamics,
- physically meaningful scrambling-strength control,
- numerically careful projective measurement,
- strict subsystem radiation definition,
- and an operationally correct **optimal decoder** based on trace distance.

## In relation to the black hole information paradox

It does **not** solve the paradox directly, because it contains no gravity.

But it is highly useful as a **conceptual and computational laboratory** for testing ideas central to the paradox:

- scrambling,
- subsystem information recovery,
- decoherence vs unitary evolution,
- and the emergence or loss of recoverable information in radiation.

So the most honest conclusion is:

> This script is not a black hole model in the gravitational sense, but it is a meaningful quantum-information toy model for one of the paradox’s central operational questions: when and how information becomes recoverable from radiation-like subsystems.

# Monitored Random Unitary Circuits and Information Recovery

A single-file Python simulation of **measurement-induced information recovery transitions** in a **monitored brickwork random unitary circuit**.

This project tests whether a logical qubit, initially localized in one qubit of an \(N\)-qubit system, can still be **optimally recovered from a radiation subsystem** after scrambling dynamics and stochastic measurements. The operational quantity of interest is the **optimal recovery fidelity** from half of the system.

The model is motivated by ideas from:

- **quantum information scrambling**
- **measurement-induced phase transitions**
- **quantum error correction**
- **black hole information recovery / Hayden-Preskill-type intuition**
- the broader **black hole information paradox**

It is **not** a gravity simulation and does **not** model horizons, Hawking radiation, or spacetime geometry directly. It is a controlled quantum information toy model.

---

## What this code does

The script simulates:

- an \(N\)-qubit pure state,
- a single logical qubit initially encoded in qubit 0,
- a **brickwork nearest-neighbor random unitary circuit**,
- repeated **probabilistic computational-basis measurements** after each layer,
- and an **optimal decoder** that asks:

> How well can an observer recover the original logical qubit using only the radiation subsystem consisting of the second half of the qubits?

The key hypothesis is that the average optimal recovery fidelity behaves like an order parameter:

- **low measurement rate** \(\lambda < \lambda_c\): information survives scrambling and can be recovered from radiation, so \(F_{\mathrm{opt}} \to 1\)
- **high measurement rate** \(\lambda > \lambda_c\): measurements destroy coherence and decoding becomes impossible, so \(F_{\mathrm{opt}} \to 1/2\)

Here \(1/2\) is the fidelity of a **random guess** for a qubit.

The script also attempts a **finite-size scaling collapse** using

\[
F(\lambda, N) = f\left((\lambda - \lambda_c) N^{1/\nu}\right)
\]

and computes a susceptibility-like diagnostic

\[
-dF/d\lambda
\]

to estimate the transition region.

---

## Why this is interesting

This project explores a concrete operational question:

> When does scrambling spread information into accessible subsystems, and when do measurements destroy it faster than it can be recovered?

That question is central to several modern themes:

- **scrambling vs decoherence**
- **error correction in many-body systems**
- **emergent code subspaces**
- **recoverability of information from subsystems**
- **toy models of black hole information flow**

In black hole language, one may loosely interpret:

- the full system as a toy “black hole + environment” quantum system,
- the radiation subsystem as “accessible outgoing degrees of freedom,”
- random unitary dynamics as a toy model for **chaotic scrambling**,
- measurements as an effective model of **information loss, decoherence, or environment-induced collapse**,
- recovery fidelity as an operational measure of whether the information about an infallen qubit is still encoded in accessible radiation.

---

## Repository structure

This repository is intentionally minimal.

```text
.
├── main.py
├── README.md
```

Generated outputs after running:

```text
plot_fidelity_vs_lambda.png
plot_finite_size_collapse.png
plot_susceptibility.png
```

---

## Requirements

- Python 3.9+
- `numpy`
- `scipy`
- `matplotlib`

Install dependencies with:

```bash
pip install numpy scipy matplotlib
```

---

## Running

Run the simulation with:

```bash
python main.py
```

The script will:

1. simulate fidelity data over a grid of system sizes and measurement rates,
2. print progress during Monte Carlo sampling,
3. estimate finite-size scaling parameters \(\lambda_c\) and \(\nu\),
4. generate and save 3 plots.

---

## Default simulation settings

The default script uses:

- system sizes:
  ```python
  N = [4, 6, 8, 10]
  ```

- measurement rates:
  ```python
  lambda in [0, 1] with 21 points
  ```

- Monte Carlo samples:
  ```python
  200 trials per (N, lambda)
  ```

- circuit depth:
  ```python
  2*N layers
  ```

- radiation subsystem:
  ```python
  qubits N//2 through N-1
  ```

These values were chosen to remain feasible for a **full state-vector simulation**.

---

## Model details

## 1. Initial state

A logical qubit is drawn at random:

\[
|\psi\rangle = \alpha |0\rangle + \beta |1\rangle
\]

with Haar-random coefficients \(\alpha,\beta\).

It is embedded into the full system as:

\[
|\psi_{\text{full}}\rangle = \alpha |00\cdots0\rangle + \beta |10\cdots0\rangle
\]

So:

- qubit 0 carries the logical information,
- all other qubits start in \(|0\rangle\).

---

## 2. Brickwork circuit

Each layer applies nearest-neighbor two-qubit unitaries in a brickwork pattern:

- even layers:
  \[
  (0,1), (2,3), (4,5), \dots
  \]

- odd layers:
  \[
  (1,2), (3,4), (5,6), \dots
  \]

Each two-qubit gate is generated as

\[
U = e^{i \, s H}
\]

where:

- \(s\) is the gate strength,
- \(H\) is a random Hermitian \(4\times4\) matrix,
- \(H\) is normalized by Frobenius norm.

This matters because it makes the **strength parameter physically meaningful**. A QR decomposition of a scaled random matrix would not have done that.

---

## 3. Measurements

After every unitary layer:

- each qubit is measured independently with probability \(\lambda\),
- measurement is in the computational basis,
- measured qubits are processed in **random order** to avoid systematic bias.

Measurements project the state and renormalize it.

The code includes numerical safeguards:

- clips probabilities to nonnegative values,
- checks for near-zero norms before normalization.

---

## 4. Radiation subsystem

The radiation subsystem is defined as the **second half** of the qubits:

\[
R = \{N/2, N/2+1, \dots, N-1\}
\]

This is always a **strict subset** of the full system.

That restriction is crucial. If one incorrectly used the whole system as the “radiation” subsystem, the reduced state would just be the full pure state, which would trivialize the decoding problem.

---

## 5. Optimal decoder

The script does **not** use an ad hoc projection decoder.

Instead, it uses the correct operational quantity for distinguishing the encoded logical alternatives after scrambling.

### Decoder recipe

Given the final pure state:

1. Project onto logical qubit 0 = \(|0\rangle\) and \(|1\rangle\), obtaining conditional states \(|\psi_0\rangle\), \(|\psi_1\rangle\).
2. Compute their weights:
   \[
   p_0 = \|\psi_0\|^2, \quad p_1 = \|\psi_1\|^2
   \]
3. Normalize them.
4. Compute reduced density matrices on radiation:
   \[
   \rho_R^0,\ \rho_R^1
   \]
5. Form:
   \[
   \sigma_0 = |\alpha|^2 p_0 \rho_R^0,\quad
   \sigma_1 = |\beta|^2 p_1 \rho_R^1
   \]
6. Compute trace distance:
   \[
   T = \frac12 \sum_i |\lambda_i(\sigma_0 - \sigma_1)|
   \]
7. Define:
   \[
   F_{\mathrm{opt}} = \min(0.5 + T,\ 1.0)
   \]

If either branch has nearly zero weight, the code returns:

\[
F_{\mathrm{opt}} = 0.5
\]

which corresponds to no recoverable information.

---

## Outputs

## Plot 1: Fidelity vs measurement rate

Shows:

- \(F_{\mathrm{opt}}\) vs \(\lambda\),
- one curve per system size,
- error bars = standard error of the mean,
- dashed line at \(F=0.5\) labeled “random guess”.

Expected qualitative signature:

- lower \(\lambda\): higher fidelity
- higher \(\lambda\): fidelity approaches 0.5
- larger \(N\): sharper crossover if a transition is present

---

## Plot 2: Finite-size scaling collapse

Attempts to fit the scaling form

\[
F(\lambda, N)=f((\lambda-\lambda_c)N^{1/\nu})
\]

by scanning over candidate values of:

- critical measurement rate \(\lambda_c\)
- exponent \(\nu\)

A good collapse is suggestive, though not definitive, evidence of critical scaling.

---

## Plot 3: Susceptibility / derivative

Plots:

\[
-dF/d\lambda
\]

The peak location gives a finite-size estimate of the transition point. If the system exhibits a phase transition, one expects:

- the peak to sharpen with increasing \(N\),
- the peak height to grow,
- the peak location to drift toward \(\lambda_c\).

---

## Scientific interpretation

## What a positive result would mean

If the simulation shows:

- \(F_{\mathrm{opt}}\approx 1\) for small \(\lambda\),
- \(F_{\mathrm{opt}}\approx 0.5\) for large \(\lambda\),
- size-dependent sharpening,
- and reasonable finite-size collapse,

then the model supports the claim that there is a **measurement-induced transition in recoverability**.

Operationally, this means there is a threshold between two regimes:

### 1. Scrambling-dominated regime
Unitary dynamics spread the logical information nonlocally into the system, allowing it to be reconstructed from a subsystem.

### 2. Measurement-dominated regime
Measurements destroy the coherence needed for that information to remain encoded in accessible correlations.

This is a concrete example of a transition not in ordinary thermodynamic order, but in **quantum information recoverability**.

---

## Connection to the black hole information paradox

## Short version

This code does **not solve the black hole information paradox**.

But it is useful because it probes a question that lies close to its modern quantum-information reformulations:

> Under what conditions does information that was initially localized become recoverable from an outgoing subsystem after chaotic evolution?

That is conceptually related to black hole evaporation and information recovery.

---

## Why black hole people care about models like this

Modern understanding of the black hole information problem increasingly uses tools from:

- quantum information theory,
- entanglement structure,
- random circuit models,
- quantum error correction,
- subsystem recovery,
- replica/statistical mechanics mappings,
- Page-curve reasoning,
- Hayden-Preskill decoding ideas.

In those frameworks, black holes are not treated merely as thermodynamic absorbers, but as **scramblers** that may encode information in highly nonlocal ways.

This code studies exactly the tension between:

- **scrambling**, which hides information globally but can make it recoverable from sufficiently large subsystems,
- and **measurement/decoherence**, which can erase the very correlations needed for recovery.

That tension is highly relevant to black hole information discussions because any real recovery mechanism must survive interactions with environments, coarse graining, and partial access.

---

## More concrete conceptual analogies

This model loosely mirrors the following black-hole-inspired narrative:

### Logical qubit
Represents a small quantum message thrown into a chaotic system.

### Brickwork random circuit
Represents chaotic scrambling dynamics.

### Radiation subsystem
Represents the accessible outgoing degrees of freedom available to an external observer.

### Optimal decoder
Represents the best physically allowed recovery protocol from what is accessible.

### Measurements
Represent decohering processes, leakage to inaccessible sectors, or effective collapse channels that compete with coherent scrambling.

This does not reproduce actual semiclassical evaporation, but it tests whether **recoverability itself** can undergo a sharp transition.

---

## Possible implications for black hole information research

## 1. Recoverability can be an operational order parameter

One of the most useful aspects of this model is that it uses a directly operational quantity:

\[
F_{\mathrm{opt}}
\]

This is better than qualitative claims like “information is somewhere in the wavefunction.” It asks whether information is **actually recoverable from a subsystem**.

For black hole discussions, this is valuable because the paradox is not just about whether global evolution is unitary in principle, but whether information is encoded in an accessible way in Hawking radiation.

---

## 2. Scrambling alone is not enough

A system may scramble very efficiently, but if sufficiently strong measurements or decoherence continually collapse important correlations, subsystem recovery can fail.

This is relevant because many naive statements about black holes say:

> “The information is scrambled, therefore it is preserved.”

But in operational settings, preservation is not the same as **recoverability from accessible radiation**.

This model highlights that subtle distinction.

---

## 3. Phase-transition language may help organize the problem

Measurement-induced transitions have taught the field that quantum information properties can exhibit sharp, universal changes as competing processes are tuned.

If analogous mechanisms matter in gravitational contexts, then information recovery may be better framed not as a binary yes/no issue, but as a competition between:

- scrambling,
- leakage,
- measurement,
- coarse graining,
- observer access,
- and subsystem size.

That perspective could be useful even if the exact circuit model is not fundamental.

---

## 4. Supports quantum error correction intuition

The low-measurement regime resembles a dynamically generated code:

- logical information is hidden nonlocally,
- local disturbances do not fully erase it,
- a sufficiently large subsystem can decode it.

That is very much in the spirit of the view that black hole interiors and bulk reconstruction may be understood through **quantum error-correcting structure**.

---

## 5. Helps test robustness of “Page-curve-like” intuition

The modern resolution of the black hole information paradox often involves entanglement wedge ideas, islands, and Page-curve reasoning.

This code does not model those directly. However, it does study whether information becomes accessible in a subsystem as a function of system size and decohering effects. That is conceptually adjacent to the same family of questions:

- when does radiation contain the information?
- how robust is that statement to noise and monitoring?
- what destroys decoding?

---

## What this model cannot tell us

It is important not to overclaim.

This project does **not** include:

- gravity,
- curved spacetime,
- Hawking pair creation,
- energy conservation constraints of evaporation,
- horizon dynamics,
- semiclassical backreaction,
- wormholes,
- islands,
- replica wormhole physics,
- gauge constraints,
- locality in spacetime geometry,
- AdS/CFT correspondence.

So it cannot by itself establish any definitive claim about real astrophysical or semiclassical black holes.

It is a **quantum information toy model**.

---

## Best way to describe its relevance

A fair statement is:

> This simulation does not solve the black hole information paradox, but it provides a clean operational laboratory for studying when scrambled quantum information remains recoverable from accessible subsystems in the presence of monitoring or decoherence. That makes it useful as a toy model for aspects of black-hole-inspired information recovery.

That is the right level of ambition.

---

## Known limitations

## 1. Small system sizes

This is a full state-vector simulation, so the Hilbert space grows as:

\[
2^N
\]

Practical default sizes are:

- 4
- 6
- 8
- 10

Maybe 12 with patience.

This is enough to see trends, but not enough for high-precision critical exponent extraction.

---

## 2. Finite-size scaling is only indicative

The code performs a simple grid search for \(\lambda_c\) and \(\nu\). This is useful for visualization, but should not be mistaken for a precision critical scaling analysis.

For publishable exponent estimates, one would want:

- larger systems,
- more samples,
- better collapse objective functions,
- possibly crossing analysis,
- bootstrap uncertainties,
- and comparison to known monitored-circuit universality classes.

---

## 3. Measurements are projective and basis-fixed

This model uses computational-basis projective measurements. Other monitoring schemes may produce different behavior.

Possible variations include:

- weak measurements
- random measurement bases
- continuously monitored dynamics
- ancilla-assisted measurements
- erasure channels instead of projections

---

## 4. One-dimensional nearest-neighbor geometry

The circuit is 1D and local. Black-hole-inspired scrambling is often modeled by all-to-all or fast-scrambling circuits.

This 1D brickwork choice is deliberate because it is physically structured and computationally manageable, but it is not a fast scrambler in the black-hole sense.

---

## 5. Radiation is chosen geometrically, not dynamically emitted

The “radiation subsystem” here is just the second half of the chain. It is not dynamically emitted radiation.

A more black-hole-like model might involve:

- explicit emission of ancilla qubits,
- moving qubits into an external register,
- or coupling a central system to outgoing modes.

---

## Potential extensions

This repository is a good starting point for many follow-up projects.

## Physics extensions

- add **ancilla radiation qubits** that are emitted over time
- study **Page-time-like decoding** rather than a fixed half-system partition
- compare **local 1D circuits** with **all-to-all fast scramblers**
- replace projective measurements with **dephasing**, **erasure**, or **amplitude damping**
- vary circuit depth and gate strength
- compare to **entanglement entropy** and **mutual information**
- study **tripartite information** or **decoupling diagnostics**
- test different subsystem fractions besides \(N/2\)
- analyze disorder averaging more systematically

## Numerical extensions

- parallelize trials across CPU cores
- cache gate ensembles
- use JIT tools or tensor-network methods
- support command-line flags for:
  - system sizes
  - number of lambda points
  - trial count
  - gate strength
  - random seed

## Decoder extensions

- compare optimal trace-distance decoding to:
  - pretty good measurement
  - learned decoders
  - tomography-based decoders
  - channel reconstruction methods

---

## Reproducibility

The script uses:

```python
np.random.default_rng(42)
```

for reproducibility.

Because the simulation is stochastic, exact numerical values depend on:

- system size,
- lambda grid,
- number of trials,
- gate strength,
- and random seed.

---

## Expected qualitative behavior

If the claim is supported, you should typically see:

- high fidelity at low measurement rates,
- a crossover region at intermediate \(\lambda\),
- fidelity trending toward \(0.5\) at high measurement rates,
- larger systems showing a steeper drop,
- susceptibility peaks sharpening with \(N\),
- a nontrivial approximate data collapse.

If the curves do not show this, possible reasons include:

- too few trials,
- gate strength not strong enough,
- system sizes too small,
- finite-size effects dominating,
- or the chosen model not exhibiting a clean transition in this observable.

---

## Example interpretation of outcomes

### If \(F_{\mathrm{opt}}\) stays high for all \(\lambda\)
Then measurements are not strong enough, the chosen circuit depth is too short, or the observable is too forgiving.

### If \(F_{\mathrm{opt}}\) quickly drops to 0.5 for all \(\lambda > 0\)
Then measurements dominate too strongly, or the system never forms a robust scrambled code.

### If curves cross and sharpen with \(N\)
That is the most interesting regime and consistent with a genuine finite-size precursor of a transition.

---

## Suggested citation language

If you use this repository in a report, a fair summary is:

> We numerically studied a monitored brickwork random unitary circuit in which a logical qubit is initially localized and later decoded from a half-system radiation subsystem using an optimal trace-distance decoder. The average recovery fidelity serves as an operational order parameter for the competition between scrambling and measurement, and exhibits finite-size behavior consistent with a measurement-induced transition.

---

## Disclaimer

This repository is for research exploration, pedagogy, and toy-model experimentation.

It should **not** be described as:

- a proof about real black holes,
- a direct simulation of Hawking radiation,
- or a resolution of the black hole information paradox.

It **is** useful as:

- a recoverability experiment,
- a monitored quantum dynamics toy model,
- and a bridge between quantum information language and black-hole-inspired questions.

---

## Bottom line

This project is useful because it asks a sharp, operational question:

> Can information that has been scrambled through chaotic dynamics still be decoded from an accessible subsystem when measurements continuously compete with coherence?

That question sits very close to the conceptual heart of modern information-theoretic approaches to the black hole information paradox.

So while this code does not solve the paradox, it can help illuminate:

- when subsystem recovery works,
- when it fails,
- how decoherence competes with scrambling,
- and why quantum error correction and recoverability are so central to the black-hole information story.

---

## License suggestion

If you want to add a license, MIT is a simple default:

```text
MIT License
```

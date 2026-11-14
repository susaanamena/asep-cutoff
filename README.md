# Cutoff in the ASEP: simulating the bounds on finite segments

Code from my B.Sc. thesis in Mathematics at the **University of Bonn** (2025).

The asymmetric simple exclusion process (ASEP) on a segment of length $N$ shows a **cutoff phenomenon**. The total variation distance to equilibrium drops abruptly from almost 1 to almost 0 within a window of order $N^{1/3}$. Bufetov and Nejjar proved that the limiting cutoff profile is $1 - F_{\mathrm{GUE}}$, the Tracy–Widom distribution, and gave explicit upper and lower bounds for the distance to equilibrium.

Those results are asymptotic ($N \to \infty$). This project asks how well they describe **finite systems of realistic size**. I simulate both bounds for $N$ between 50 and 2000 with up to 1,000,000 independent runs per configuration and compare them with the Tracy–Widom profile.

## Main results

All simulations use particle density $\alpha = 0.5$ and right-jump probability $p = 0.75$.

**The upper bound is already very accurate for small systems.** The upper bound comes from the hitting time of the maximal configuration. Its mean squared error against $1 - F_{\mathrm{GUE}}$ decreases steadily with $N$, and from $N = 1000$ the two curves are almost indistinguishable.

| $N$ | 50 | 100 | 200 | 500 | 1000 | 2000 |
|---|---|---|---|---|---|---|
| MSE | 5.0 · 10⁻⁴ | 2.7 · 10⁻⁴ | 1.6 · 10⁻⁴ | 8.5 · 10⁻⁵ | 5.3 · 10⁻⁵ | 4.3 · 10⁻⁵ |

**The lower bound is rougher but captures the transition.** The lower bound depends on a parameter $l_0$ (how far the leftmost particle is from its final position). The best choice depends on the system size. $l_0 = 1$ works best for $N \le 100$, while $l_0 = 2$ becomes better for $N = 2000$. This suggests that larger values of $l_0$ become optimal as $N$ grows.

**The stationary position of the leftmost particle follows Liggett's theorem.** Its empirical distribution matches $\rho(x) = c\,(q/p)^x / (1 + c\,(q/p)^x)$, with a fitted constant $c = 1.12$.

## How the simulation works

The process runs in continuous time. At each step the simulation waits an $\mathrm{Exp}(N)$ time, picks one of the $N-1$ bonds uniformly at random and attempts a jump to the right with probability $p$ or to the left with probability $q = 1-p$, respecting the exclusion rule. The position of the leftmost particle is tracked throughout, since both bounds depend only on it.

- **Upper bound.** Start from the minimal configuration (all particles on the left) and record the time until the maximal configuration is reached. After rescaling, the tail probabilities $P(Z > c)$ give the upper bound.
- **Stationary distribution.** Start from the maximal configuration and run the process until twice the cutoff time, then record the gap of the leftmost particle.
- **Lower bound.** Start from the minimal configuration and record whether the event $A_N(l_0)$ holds at times $g(k,c)$ for $c \in [-12, 4]$ in steps of 0.1. Subtracting the stationary probability gives the lower bound.

The simulations are written in C++ because Python was too slow for more than about 100,000 trials. The analysis and the plots were done in Wolfram Mathematica.

## Repository structure

```
src/
  ASEPSimulation.hpp   class declaration
  ASEPSimulation.cpp   simulation methods for both bounds
  main.cpp             command-line program that runs all simulations
figures/               plots of the main results
```

## Build and run

You only need a C++17 compiler.

```bash
g++ -O3 -std=c++17 src/main.cpp src/ASEPSimulation.cpp -o asep
./asep -N 100 -a 0.5 --p 0.75 -t 10000
```

| Option | Meaning | Default |
|---|---|---|
| `-N` | segment length | 100 |
| `-a`, `--alpha` | particle density, $k = \lceil N\alpha \rceil$ | 0.5 |
| `--p` | right-jump probability | 0.75 |
| `-t`, `--trials` | number of independent runs | 10000 |

The program writes its results as `.dat` files in the current folder.

| File | Content |
|---|---|
| `hittingtimes{N}.dat` | raw hitting times |
| `ht_prob{N}.dat` | upper bound $P(Z > c)$ for each value of $c$ |
| `stationary{N}.dat` | gap of the leftmost particle at twice the cutoff time |
| `prob{N}_L{l0}.dat` | occurrence of $A_N(l_0)$ (0/1) for each run and each $c$ |
| `rel.probs{N}_L{l0}.dat` | empirical probability of $A_N(l_0)$ for each $c$ |

The raw data files are opened in append mode, so running the program again adds new trials to the existing ones. Delete the `.dat` files to start from scratch. Large values of $N$ with many trials can take several hours.

## References

1. A. Bufetov and P. Nejjar, *Cutoff profile of ASEP on a segment*, Probab. Theory Related Fields 183 (2022), 229–253.
2. C. Labbé and H. Lacoin, *Cutoff phenomenon for the asymmetric simple exclusion process and the biased card shuffling*, Ann. Probab. 47 (2019), 1541–1586.
3. D. A. Levin and Y. Peres, *Markov chains and mixing times*, 2nd ed., AMS, 2017.
4. T. M. Liggett, *Stochastic interacting systems: contact, voter and exclusion processes*, Springer, 1999.

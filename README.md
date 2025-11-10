# R&D / AI Assignment – Curve Parameter Recovery

This repository solves for the unknown parameters \( \theta, M, X \) in:
\[
x = t\cos\theta - e^{M|t|}\sin(0.3t)\sin\theta + X,\quad
y = 42 + t\sin\theta + e^{M|t|}\sin(0.3t)\cos\theta,
\]
given sample points \((x_i,y_i)\) for \(6<t<60\).

## How to use (super simple)
1. Put your data file `xy_data.csv` (columns: `x,y`) in the same folder as `fit_params.py`.
2. Run:
   ```bash
   python fit_params.py
   ```
3. The script prints:
   - Estimated \(\theta\) (degrees), \(M\), and \(X\)
   - An **L1 score** on a uniform \(t\)-grid (as per the rubric)
   - A **ready-to-paste** submission string for your README/Desmos.

> **Optional:** `python fit_params.py --plot` will show the rotated data vs the fitted curve in \((t,s)\)-space.

## Method (what you can explain in your report)
1. **Undo the transform:** Translate points by \((-X,-42)\) and rotate by \(-\theta\). In the correct frame, the cloud becomes \([t,\ s]\) where \(s = e^{M|t|}\sin(0.3t)\).
2. **Recover \(M\) analytically:** From \( |s/\sin(0.3t)| = e^{M|t|} \Rightarrow \log|s/\sin(0.3t)| = M|t| \), we get a simple linear fit for \(M\) once \((t,s)\) are known.
3. **Search over \(\theta, X\):** We do a coarse grid in the allowed ranges \((0^\circ,50^\circ)\) and \((0,100)\), compute \((t,s)\) for each, fit \(M\), and score the residual \( \| s - e^{M|t|}\sin(0.3t)\|^2 \). Then we refine locally.
4. **L1 scoring:** We build a uniform grid in \(t\) within both the observed and required ranges, interpolate the “expected” \(s(t)\) from data in the rotated frame, and compare against the model’s \(s(t)\). Both are mapped back to \((x,y)\) and we report the mean \(L1 = \mathbb{E}[|Δx|+|Δy|]\).

## Submission
=== Estimated Parameters ===
theta (degrees): 29.996445
M               : 0.030012
X               : 54.995622
Example (format only):
```
(t*cos(0.523537) - exp(0.030012*abs(t))*sin(0.3*t)*sin(0.523537) + 54.995622, 42 + t*sin(0.523537) + exp(0.030012*abs(t))*sin(0.3*t)*cos(0.523537))
```

## Notes
- Constraints are enforced: \(0^\circ<\theta<50^\circ\), \(-0.05<M<0.05\), \(0<X<100\).
- CSV headers must be `x,y`. If your first two columns are the data but named differently, the script will auto-rename them.
- If you want raw residuals or more metrics, say the word and we’ll add them.

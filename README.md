# R&D / AI Assignment – Curve Parameter Recovery

This repository solves for the unknown parameters \( \theta, M, X \) in:
\[
x = t\cos\theta - e^{M|t|}\sin(0.3t)\sin\theta + X,\quad
y = 42 + t\sin\theta + e^{M|t|}\sin(0.3t)\cos\theta,
\]
given sample points \((x_i,y_i)\) for \(6<t<60\).


## Method 
1. **Undo the transform:** Translate points by \((-X,-42)\) and rotate by \(-\theta\). In the correct frame, the cloud becomes \([t,\ s]\) where \(s = e^{M|t|}\sin(0.3t)\).
2. **Recover \(M\) analytically:** From \( |s/\sin(0.3t)| = e^{M|t|} \Rightarrow \log|s/\sin(0.3t)| = M|t| \), we get a simple linear fit for \(M\) once \((t,s)\) are known.
3. **Search over \(\theta, X\):** We do a coarse grid in the allowed ranges \((0^\circ,50^\circ)\) and \((0,100)\), compute \((t,s)\) for each, fit \(M\), and score the residual \( \| s - e^{M|t|}\sin(0.3t)\|^2 \). Then we refine locally.
4. **L1 scoring:** We build a uniform grid in \(t\) within both the observed and required ranges, interpolate the “expected” \(s(t)\) from data in the rotated frame, and compare against the model’s \(s(t)\). Both are mapped back to \((x,y)\) and we report the mean \(L1 = \mathbb{E}[|Δx|+|Δy|]\).

## Submission
=== Estimated Parameters ===

theta (degrees): 29.996445

M               : 0.030012

X               : 54.995622
```
(t*cos(0.523537) - exp(0.030012*abs(t))*sin(0.3*t)*sin(0.523537) + 54.995622, 42 + t*sin(0.523537) + exp(0.030012*abs(t))*sin(0.3*t)*cos(0.523537))
```

Desmos link
```
https://www.desmos.com/calculator/b81saq0vux
```

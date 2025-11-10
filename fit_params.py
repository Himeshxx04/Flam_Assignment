#!/usr/bin/env python3
import argparse, math, sys
from pathlib import Path
import numpy as np

import pandas as pd

def rotate_minus_theta(points, theta):
    # R(-θ) = [[cosθ, sinθ], [-sinθ, cosθ]]
    c, s = math.cos(theta), math.sin(theta)
    R = np.array([[c, s], [-s, c]])
    return points @ R.T

def estimate_M_from_ts_si(t, s):
    sinp = np.sin(0.3 * t)
    mask = (np.abs(sinp) > 1e-6) & (np.isfinite(s)) & (np.abs(s) > 1e-12)
    if mask.sum() < 5:
        return None, np.inf
    a = np.abs(s[mask] / sinp[mask])
    m2 = np.isfinite(a) & (a > 0) & np.isfinite(np.abs(t[mask]))
    if m2.sum() < 5:
        return None, np.inf
    T = np.abs(t[mask][m2])
    Y = np.log(a[m2])
    denom = float(np.dot(T, T))
    if denom <= 0:
        return None, np.inf
    M_est = float(np.dot(T, Y) / denom)
    rmse = float(np.sqrt(np.mean((Y - M_est*T)**2)))
    return M_est, rmse

def objective(theta_deg, X, xy):
    theta = math.radians(theta_deg)
    shifted = xy - np.array([X, 42.0])
    ts = rotate_minus_theta(shifted, theta)
    t = ts[:,0]; s = ts[:,1]
    order = np.argsort(t)
    t = t[order]; s = s[order]
    M_est, _ = estimate_M_from_ts_si(t, s)
    if M_est is None or not (-0.05 < M_est < 0.05):
        return np.inf, None, None, None
    s_pred = np.exp(M_est * np.abs(t)) * np.sin(0.3 * t)
    mse = float(np.mean((s - s_pred)**2))
    return mse, M_est, t, s

def fit_params(xy, seed=0):
    thetas = np.linspace(0.5, 49.5, 180)
    Xs = np.linspace(0.5, 99.5, 180)
    best = {'mse': np.inf, 'theta': None, 'X': None, 'M': None}
    for th in thetas:
        for X in Xs:
            mse, M_est, _, _ = objective(th, X, xy)
            if mse < best['mse']:
                best.update({'mse':mse, 'theta':th, 'X':X, 'M':M_est})
    # refine
    rng = np.random.default_rng(seed)
    theta = best['theta']; X = best['X']; M = best['M']
    mse = best['mse']
    for _ in range(3000):
        th_try = float(np.clip(theta + rng.normal(scale=0.05), 0.01, 49.99))
        X_try  = float(np.clip(X + rng.normal(scale=0.10), 0.01, 99.99))
        mse_try, M_try, _, _ = objective(th_try, X_try, xy)
        if mse_try < mse and (M_try is not None):
            theta, X, M, mse = th_try, X_try, M_try, mse_try
    return theta, M, X, mse

def interpolate_lin(x_src, y_src, x_tgt):
    # simple linear interpolation within range
    return np.interp(x_tgt, x_src, y_src)

def l1_score(theta_deg, M, X, xy, t_min=6.0, t_max=60.0, n=500):
    """Compute L1 distance on a uniform t-grid by comparing the expected curve
    reconstructed from rotated data (interpolated s(t)) vs model prediction.
    """
    theta = math.radians(theta_deg)
    shifted = xy - np.array([X, 42.0])
    ts = rotate_minus_theta(shifted, theta)
    t_data = ts[:,0]; s_data = ts[:,1]
    order = np.argsort(t_data)
    t_data = t_data[order]; s_data = s_data[order]

    # Restrict to valid observed t-range intersected with [t_min, t_max]
    t_lo = max(t_min, float(np.min(t_data)))
    t_hi = min(t_max, float(np.max(t_data)))
    if t_hi <= t_lo:
        return None

    t_grid = np.linspace(t_lo, t_hi, n)
    # expected curve in (x,y): from data s via interpolation
    s_exp = interpolate_lin(t_data, s_data, t_grid)
    # model predicted s
    s_pred = np.exp(M * np.abs(t_grid)) * np.sin(0.3 * t_grid)

    # Map back to (x,y) for both expected and predicted
    th = theta
    c, s = math.cos(th), math.sin(th)
    # forward rotation R(θ) = [[c, -s],[s, c]] and translation (X,42)
    x_exp = t_grid * c - s_exp * s + X
    y_exp = 42.0 + t_grid * s + s_exp * c
    x_pred = t_grid * c - s_pred * s + X
    y_pred = 42.0 + t_grid * s + s_pred * c

    l1 = np.mean(np.abs(x_exp - x_pred) + np.abs(y_exp - y_pred))
    return float(l1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="xy_data.csv", help="Path to xy_data.csv (columns: x,y)")
    ap.add_argument("--plot", action="store_true", help="Show plot (requires matplotlib)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        print(f"ERROR: CSV not found at {csv_path}. Place your file there or pass --csv path.")
        sys.exit(2)

    df = pd.read_csv(csv_path)
    if not {'x','y'}.issubset(df.columns):
        cols = list(df.columns)
        if len(cols) < 2:
            print("CSV must have at least two columns (x,y).")
            sys.exit(3)
        df = df.rename(columns={cols[0]:'x', cols[1]:'y'})
    xy = df[['x','y']].to_numpy(dtype=float)

    theta_deg, M, X, mse = fit_params(xy, seed=args.seed)
    print("=== Estimated Parameters ===")
    print(f"theta (degrees): {theta_deg:.6f}")
    print(f"M               : {M:.6f}")
    print(f"X               : {X:.6f}")
    print(f"MSE(s-residual) : {mse:.6e}")

    # L1 score (uniform t grid)
    L1 = l1_score(theta_deg, M, X, xy)
    if L1 is not None:
        print(f"L1 score (uniform t grid): {L1:.6f}")
    else:
        print("L1 score: not computed (insufficient overlapping t-range)")

    # Submission string
    th_rad = math.radians(theta_deg)
    sub = (
        f"(t*cos({th_rad:.6f}) - exp({M:.6f}*abs(t))*sin(0.3*t)*sin({th_rad:.6f}) + {X:.6f}, "
        f"42 + t*sin({th_rad:.6f}) + exp({M:.6f}*abs(t))*sin(0.3*t)*cos({th_rad:.6f}))"
    )
    print("\nSubmission string (paste this in README/Desmos):")
    print(sub)

    if args.plot:
        try:
            import matplotlib.pyplot as plt
            # Rotate to (t,s) for visualization
            shifted = xy - np.array([X, 42.0])
            ts = rotate_minus_theta(shifted, math.radians(theta_deg))
            t = ts[:,0]; s = ts[:,1]
            order = np.argsort(t); t = t[order]; s = s[order]
            s_fit = np.exp(M * np.abs(t)) * np.sin(0.3 * t)

            plt.figure()
            plt.plot(t, s, '.', label='data (rotated)')
            plt.plot(t, s_fit, '-', label='model fit')
            plt.xlabel('t'); plt.ylabel('s')
            plt.legend()
            plt.title('Rotated data vs model in (t, s)')
            plt.show()
        except Exception as e:
            print("Plotting failed (matplotlib missing or other issue):", e)

if __name__ == "__main__":
    main()

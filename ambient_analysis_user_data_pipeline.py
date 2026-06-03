#!/usr/bin/env python3
"""
ambient_analysis_user_data_pipeline.py

Purpose
-------
Modified version of the MATLAB/Python ambient vibration pipeline for USER DATA.

Main change from the synthetic version:
    - The synthetic data generation block is removed.
    - Y_setup is read from Excel or CSV files.
    - The rest of the workflow remains similar:
        PSD -> singular value spectrum -> local modal identification -> global mode shape assembly.

Expected wide format per setup
------------------------------
Each setup must have one sheet/file with columns:
    time_s, acc_DOF1_mps2, acc_DOF2_mps2, ...

Example:
    Setup1: time_s, acc_DOF1_mps2, acc_DOF2_mps2, acc_DOF3_mps2, acc_DOF4_mps2
    Setup2: time_s, acc_DOF3_mps2, acc_DOF4_mps2, acc_DOF5_mps2, acc_DOF6_mps2
    Setup3: time_s, acc_DOF2_mps2, acc_DOF4_mps2, acc_DOF5_mps2, acc_DOF6_mps2

Install
-------
    pip install numpy scipy matplotlib pandas openpyxl

CSV-only mode does not need pandas/openpyxl.

Important
---------
For real OMA/BAYOMA analysis, use long enough data. The original script used 600 s at 200 Hz.
The example template is intentionally short, only to show the format.
"""

from __future__ import annotations

from pathlib import Path
import csv
import math
import warnings

import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import minimize
import matplotlib.pyplot as plt


# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Option 1: "excel"
# Option 2: "csv"
DATA_MODE = "excel"

EXCEL_FILE = "ambient_vibration_data_template.xlsx"
SHEET_NAMES = ["Setup1", "Setup2", "Setup3"]

CSV_FILES = [
    "ambient_setup1_example.csv",
    "ambient_setup2_example.csv",
    "ambient_setup3_example.csv",
]

# DOFs measured by each setup, in the same order as columns after time_s.
SETUP_DOFS = [
    [1, 2, 3, 4],
    [3, 4, 5, 6],
    [2, 4, 5, 6],
]

N_DOF_GLOBAL = 6
REF_DOF = 4

# Analysis frequency band
F_L = 0.5
F_U = 20.0
MF = 400  # number of intervals in the target frequency band

# Target mode for local identification.
# For real data, use the dominant peak from the SV spectrum as a guide.
TARGET_FN_GUESS_HZ = 2.0
TARGET_BAND_LOWER_HZ = 1.6
TARGET_BAND_UPPER_HZ = 2.4

OUTPUT_DIR = "ambient_analysis_output"


# =============================================================================
# DATA LOADING
# =============================================================================

def _estimate_fs_from_time(time_s: np.ndarray) -> float:
    dt = np.diff(time_s)
    dt = dt[np.isfinite(dt)]
    if len(dt) == 0:
        raise ValueError("Cannot estimate sampling frequency because time_s has too few samples.")
    dt_med = float(np.median(dt))
    if dt_med <= 0:
        raise ValueError("time_s must be strictly increasing.")
    # Check if nearly uniform
    if np.std(dt) / dt_med > 0.01:
        warnings.warn(
            "The time step is not perfectly uniform. PSD/FFT assumes uniform sampling. "
            "Consider resampling your data before analysis."
        )
    return 1.0 / dt_med


def load_excel_wide(path: str | Path, sheet_names: list[str], setup_dofs: list[list[int]]):
    import pandas as pd

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Excel file not found: {path}")

    Y_setup = []
    time_setup = []

    for sheet_name, dofs in zip(sheet_names, setup_dofs):
        df = pd.read_excel(path, sheet_name=sheet_name)
        expected_cols = ["time_s"] + [f"acc_DOF{d}_mps2" for d in dofs]

        missing = [c for c in expected_cols if c not in df.columns]
        if missing:
            raise ValueError(
                f"Sheet {sheet_name} is missing columns: {missing}. "
                f"Expected columns are: {expected_cols}"
            )

        df = df[expected_cols].dropna()
        time_s = df["time_s"].to_numpy(dtype=float)
        Y = df[[f"acc_DOF{d}_mps2" for d in dofs]].to_numpy(dtype=float).T

        time_setup.append(time_s)
        Y_setup.append(Y)

    fs = _estimate_fs_from_time(time_setup[0])
    return Y_setup, time_setup, fs


def load_csv_wide(csv_files: list[str | Path], setup_dofs: list[list[int]]):
    Y_setup = []
    time_setup = []

    for csv_file, dofs in zip(csv_files, setup_dofs):
        csv_file = Path(csv_file)
        if not csv_file.exists():
            raise FileNotFoundError(f"CSV file not found: {csv_file}")

        data = np.genfromtxt(csv_file, delimiter=",", names=True, dtype=float, encoding="utf-8")
        names = list(data.dtype.names)
        expected_cols = ["time_s"] + [f"acc_DOF{d}_mps2" for d in dofs]

        missing = [c for c in expected_cols if c not in names]
        if missing:
            raise ValueError(
                f"CSV {csv_file} is missing columns: {missing}. "
                f"Expected columns are: {expected_cols}"
            )

        time_s = np.asarray(data["time_s"], dtype=float)
        Y = np.vstack([np.asarray(data[f"acc_DOF{d}_mps2"], dtype=float) for d in dofs])

        time_setup.append(time_s)
        Y_setup.append(Y)

    fs = _estimate_fs_from_time(time_setup[0])
    return Y_setup, time_setup, fs


# =============================================================================
# PSD AND MODAL IDENTIFICATION FUNCTIONS
# =============================================================================

def compute_psd_matrices(
    Y_setup: list[np.ndarray],
    fs: float,
    f_l: float,
    f_u: float,
    mf: int,
):
    dt = 1.0 / fs
    n0 = min(Y.shape[1] for Y in Y_setup)

    # Use common length across all setups.
    Y_setup = [Y[:, :n0] for Y in Y_setup]

    df_target = (f_u - f_l) / mf
    M = int(np.floor(n0 * dt * df_target))

    if M < 1:
        M = 1
        warnings.warn(
            "Data duration is short relative to MF and the frequency band. "
            "Using M=1 segment. For real OMA, use longer data or reduce MF."
        )

    Nw = int(np.floor(n0 / M))
    N = Nw * M

    f_seg = np.arange(Nw) / (Nw * dt)
    k_band = np.where((f_seg >= f_l) & (f_seg <= f_u))[0]
    f_band = f_seg[k_band]

    E_avg_setup = []
    for Y in Y_setup:
        Y = Y - np.mean(Y, axis=1, keepdims=True)
        Y_trunc = Y[:, :N]
        nr = Y.shape[0]
        E_acc = np.zeros((nr, nr, len(k_band)), dtype=np.complex128)

        for seg in range(M):
            idx0 = seg * Nw
            idx1 = idx0 + Nw
            Y_seg = Y_trunc[:, idx0:idx1]

            F_all = np.sqrt(dt / Nw) * np.fft.fft(Y_seg, n=Nw, axis=1)
            F_band = F_all[:, k_band]

            for ki in range(len(k_band)):
                Fk = F_band[:, ki].reshape(-1, 1)
                E_acc[:, :, ki] += Fk @ Fk.conj().T

        E_avg_setup.append(E_acc / M)

    return E_avg_setup, f_band, M, Nw


def compute_sv_spectrum(E_avg: np.ndarray):
    nr, _, nk = E_avg.shape
    SV = np.zeros((nr, nk))
    EV = np.zeros((nr, nr, nk))

    for ki in range(nk):
        ReE = np.real(E_avg[:, :, ki])
        ReE = 0.5 * (ReE + ReE.T)
        ev, V = np.linalg.eigh(ReE)
        idx = np.argsort(ev)[::-1]
        SV[:, ki] = ev[idx]
        EV[:, :, ki] = V[:, idx]

    return SV, EV


def dynamic_amplification(f_vec: np.ndarray, fn: float, zeta: float, q: int = 0):
    bk = f_vec / fn
    return (2 * np.pi * f_vec) ** (-2 * q) / ((1 - bk ** 2) ** 2 + (2 * zeta * bk) ** 2)


def nllf_reduced(theta, f_id, E_band, nr, nf, phi_u):
    fn = theta[0]
    zeta = abs(theta[1])
    S = np.exp(theta[2])
    Se = np.exp(theta[3])

    if fn <= 0 or zeta <= 0 or zeta >= 1 or not np.isfinite(S) or not np.isfinite(Se):
        return 1e15

    Dk = dynamic_amplification(f_id, fn, zeta, q=0)
    L = 0.0

    for ki in range(nf):
        SDk = S * Dk[ki] + Se
        if SDk <= 0 or Se <= 0:
            return 1e15

        ReEk = np.real(E_band[:, :, ki])
        tr_Ek = np.real(np.trace(ReEk))
        uAu = float(np.real(phi_u.T @ ReEk @ phi_u))
        fac = 1.0 / (1.0 + Se / (S * Dk[ki]))

        L += np.log(SDk)
        L += (1.0 / Se) * (tr_Ek - fac * uAu)

    L += (nr - 1) * nf * np.log(Se)
    return float(L)


def identify_local_modes(
    E_avg_setup: list[np.ndarray],
    f_band: np.ndarray,
    setup_dofs: list[list[int]],
    target_lower: float,
    target_upper: float,
):
    k_id = np.where((f_band >= target_lower) & (f_band <= target_upper))[0]
    if len(k_id) < 5:
        raise ValueError(
            "Too few frequency points in target identification band. "
            "Widen TARGET_BAND_LOWER/UPPER or use longer data."
        )

    f_id = f_band[k_id]

    local_fn = []
    local_zeta = []
    local_S = []
    local_Se = []
    local_phi = []

    for r, E_avg in enumerate(E_avg_setup):
        E_band = E_avg[:, :, k_id]
        nr = E_band.shape[0]
        nf = len(k_id)

        A0 = np.zeros((nr, nr))
        d_tot = 0.0
        for ki in range(nf):
            ReEk = np.real(E_band[:, :, ki])
            A0 += ReEk
            d_tot += np.real(np.trace(ReEk))

        ev0, V0 = np.linalg.eigh(A0)
        imax = int(np.argmax(ev0))
        phi0 = np.real(V0[:, imax])
        phi0 = phi0 / np.linalg.norm(phi0)

        # Initial fn from largest first singular/eigenvalue in the target band for this setup.
        SV_r, _ = compute_sv_spectrum(E_avg[:, :, k_id])
        fn0 = float(f_id[int(np.argmax(SV_r[0, :]))])

        zeta0 = 0.01
        Dk0 = dynamic_amplification(f_id, fn0, zeta0, q=0)
        k0_hat = float(phi0.T @ A0 @ phi0)

        d_hat = 0.0
        for ki in range(nf):
            d_hat += float(phi0.T @ np.real(E_band[:, :, ki]) @ phi0)

        Se0 = max((d_tot - k0_hat) / max((nr - 1) * nf, 1), 1e-20)
        S0 = max(np.sum(d_hat / Dk0) / nf, 1e-20)

        theta0 = np.array([fn0, zeta0, np.log(S0), np.log(Se0)])

        res = minimize(
            nllf_reduced,
            theta0,
            args=(f_id, E_band, nr, nf, phi0),
            method="Nelder-Mead",
            options={"maxiter": 5000, "maxfev": 5000, "xatol": 1e-10, "fatol": 1e-10},
        )

        theta = res.x
        fn_hat = float(theta[0])
        zeta_hat = float(abs(theta[1]))
        S_hat = float(np.exp(theta[2]))
        Se_hat = float(np.exp(theta[3]))

        Dk_hat = dynamic_amplification(f_id, fn_hat, zeta_hat, q=0)
        A_hat = np.zeros((nr, nr))
        for ki in range(nf):
            factor = 1.0 / (1.0 + Se_hat / (S_hat * Dk_hat[ki]))
            A_hat += factor * np.real(E_band[:, :, ki])

        ev_hat, V_hat = np.linalg.eigh(A_hat)
        imax_hat = int(np.argmax(ev_hat))
        phi_hat = np.real(V_hat[:, imax_hat])
        phi_hat = phi_hat / np.linalg.norm(phi_hat)

        local_fn.append(fn_hat)
        local_zeta.append(zeta_hat)
        local_S.append(S_hat)
        local_Se.append(Se_hat)
        local_phi.append(phi_hat)

        print(f"Setup {r+1}: fn={fn_hat:.4f} Hz, zeta={zeta_hat*100:.3f}%, "
              f"S={S_hat:.3e}, Se={Se_hat:.3e}")

    return np.array(local_fn), np.array(local_zeta), np.array(local_S), np.array(local_Se), local_phi


def assemble_global_mode_shape(local_phi, setup_dofs, n_dof_global):
    n_setups = len(setup_dofs)
    L_cells = []

    for dofs in setup_dofs:
        Lr = np.zeros((len(dofs), n_dof_global))
        for j, dof in enumerate(dofs):
            Lr[j, dof - 1] = 1.0
        L_cells.append(Lr)

    weights = np.ones(n_setups) / n_setups

    A0 = np.zeros((n_dof_global, n_dof_global))
    for r in range(n_setups):
        Lr = L_cells[r]
        tt = local_phi[r]
        nr = len(tt)
        A0 += weights[r] * (Lr.T @ (np.eye(nr) - np.outer(tt, tt)) @ Lr)

    ev, V = np.linalg.eigh(A0)
    phi_global = np.real(V[:, int(np.argmin(ev))])
    phi_global = phi_global / np.linalg.norm(phi_global)

    max_iter = 300
    tol = 1e-12

    for iteration in range(max_iter):
        phi_old = phi_global.copy()

        c_r = np.zeros(n_setups)
        lambda_r = np.zeros(n_setups)

        for r in range(n_setups):
            tr = L_cells[r] @ phi_global
            tt = local_phi[r]
            c_r[r] = np.linalg.norm(tr) ** 2
            lambda_r[r] = weights[r] * (1 - abs(tt.T @ tr) / max(np.linalg.norm(tr), 1e-15))

        A_mat = np.zeros((n_dof_global, n_dof_global))
        b_vec = np.zeros(n_dof_global)

        for r in range(n_setups):
            Lr = L_cells[r]
            tt = local_phi[r]
            A_mat += (weights[r] + lambda_r[r]) * (Lr.T @ Lr)
            b_vec += -weights[r] * math.sqrt(max(c_r[r], 0)) * (Lr.T @ tt)

        D_mat = np.block([
            [A_mat, np.outer(b_vec, b_vec)],
            [np.eye(n_dof_global), A_mat],
        ])

        evD, VD = np.linalg.eig(D_mat)
        idx = int(np.argmin(np.real(evD)))
        phi_global = np.real(VD[:n_dof_global, idx])
        phi_global = phi_global / np.linalg.norm(phi_global)

        if np.linalg.norm(phi_global - phi_old) < tol:
            break

    # Sign convention: make the reference DOF positive if possible.
    if phi_global[REF_DOF - 1] < 0:
        phi_global = -phi_global

    return phi_global


def save_vector_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


# =============================================================================
# MAIN
# =============================================================================

def main():
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)

    if DATA_MODE.lower() == "excel":
        Y_setup, time_setup, fs = load_excel_wide(EXCEL_FILE, SHEET_NAMES, SETUP_DOFS)
    elif DATA_MODE.lower() == "csv":
        Y_setup, time_setup, fs = load_csv_wide(CSV_FILES, SETUP_DOFS)
    else:
        raise ValueError("DATA_MODE must be either 'excel' or 'csv'.")

    print("Loaded user data:")
    print(f"  sampling frequency = {fs:.4f} Hz")
    for i, Y in enumerate(Y_setup):
        print(f"  Setup {i+1}: channels={Y.shape[0]}, samples={Y.shape[1]}, DOFs={SETUP_DOFS[i]}")

    E_avg_setup, f_band, M, Nw = compute_psd_matrices(Y_setup, fs, F_L, F_U, MF)
    print(f"\nPSD settings: M={M} segments, Nw={Nw} samples/segment, frequency points={len(f_band)}")

    # SV spectrum for setup 1
    SV_s1, EV_s1 = compute_sv_spectrum(E_avg_setup[0])
    df_eff = float(np.median(np.diff(f_band))) if len(f_band) > 1 else 0.0
    distance = max(1, int(round(0.5 / df_eff))) if df_eff > 0 else 1
    peaks, props = find_peaks(np.log(np.maximum(SV_s1[0, :], 1e-300)), prominence=0.5, distance=distance)
    peak_freqs = f_band[peaks]

    print("\nDetected resonance peaks from Setup 1 / SV1:")
    if len(peak_freqs) == 0:
        print("  No strong peak detected. Check frequency band, sensor placement, signal quality, or peak prominence.")
    else:
        print("  " + ", ".join(f"{p:.3f} Hz" for p in peak_freqs))

    # Plot root PSD for setup 1
    root_psd_s1 = np.zeros((len(SETUP_DOFS[0]), len(f_band)))
    for ki in range(len(f_band)):
        d = np.diag(np.real(E_avg_setup[0][:, :, ki]))
        root_psd_s1[:, ki] = np.sqrt(np.abs(d))

    plt.figure(figsize=(10, 4))
    for i, dof in enumerate(SETUP_DOFS[0]):
        plt.semilogy(f_band, root_psd_s1[i, :], label=f"DOF {dof}")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Root PSD [(m/s²)/√Hz]")
    plt.title("Root PSD Spectrum — Setup 1")
    plt.grid(True)
    plt.xlim([F_L, F_U])
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "root_psd_setup1.png", dpi=300)

    # Plot root SV for setup 1
    plt.figure(figsize=(10, 4))
    for i in range(SV_s1.shape[0]):
        plt.semilogy(f_band, np.sqrt(np.abs(SV_s1[i, :])), label=f"SV{i+1}")
    for pf in peak_freqs:
        plt.axvline(pf, linestyle=":", linewidth=0.8)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Root SV")
    plt.title("Singular Value Spectrum — Setup 1")
    plt.grid(True)
    plt.xlim([F_L, F_U])
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "sv_spectrum_setup1.png", dpi=300)

    # Local modal identification and global mode-shape assembly
    print("\nLocal modal identification:")
    local_fn, local_zeta, local_S, local_Se, local_phi = identify_local_modes(
        E_avg_setup,
        f_band,
        SETUP_DOFS,
        TARGET_BAND_LOWER_HZ,
        TARGET_BAND_UPPER_HZ,
    )

    phi_global = assemble_global_mode_shape(local_phi, SETUP_DOFS, N_DOF_GLOBAL)

    print("\nGlobal relative mode shape:")
    for dof, val in enumerate(phi_global, start=1):
        print(f"  DOF {dof}: {val:+.6f}")

    save_vector_csv(
        output_dir / "modal_results_summary.csv",
        ["setup_id", "fn_hz", "zeta_percent", "S", "Se"],
        [
            [i + 1, local_fn[i], local_zeta[i] * 100, local_S[i], local_Se[i]]
            for i in range(len(local_fn))
        ],
    )

    save_vector_csv(
        output_dir / "global_mode_shape.csv",
        ["dof", "relative_mode_shape"],
        [[i + 1, phi_global[i]] for i in range(len(phi_global))],
    )

    plt.figure(figsize=(5, 7))
    dofs = np.arange(1, N_DOF_GLOBAL + 1)
    plt.plot(np.r_[0, phi_global], np.r_[0, dofs], "o-", linewidth=2)
    plt.axvline(0, linewidth=0.8)
    plt.xlabel("Relative mode-shape component")
    plt.ylabel("DOF")
    plt.title(f"Global Mode Shape, mean fn={np.mean(local_fn):.3f} Hz")
    plt.grid(True)
    plt.yticks(dofs)
    plt.tight_layout()
    plt.savefig(output_dir / "global_mode_shape.png", dpi=300)

    print(f"\nFiles saved in: {output_dir.resolve()}")


if __name__ == "__main__":
    main()

"""
Ambient Vibration Data Analysis — Full 5-Phase Pipeline
Python translation of the provided MATLAB script.

Dependencies:
    numpy scipy matplotlib

Run:
    python ambient_analysis_full_pipeline_python.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.signal import find_peaks


def dynamic_amplification(f_vec, fn, zeta, q=0):
    """Dynamic amplification function.

    Parameters
    ----------
    f_vec : array-like
        Frequency vector in Hz.
    fn : float
        Natural frequency in Hz.
    zeta : float
        Damping ratio.
    q : int
        0 = acceleration, 1 = velocity, 2 = displacement.
    """
    f_vec = np.asarray(f_vec, dtype=float)
    bk = f_vec / fn
    return (2 * np.pi * f_vec) ** (-2 * q) / ((1 - bk**2) ** 2 + (2 * zeta * bk) ** 2)


def nllf_reduced(theta, f_id, e_band, nr, nf, phi_u):
    """Reduced negative log-likelihood function.

    theta = [fn, zeta, log(S), log(Se)]
    e_band shape: (nf, nr, nr)
    """
    fn = theta[0]
    zeta = abs(theta[1])
    S = np.exp(theta[2])
    Se = np.exp(theta[3])

    if fn <= 0 or zeta <= 0 or zeta >= 1:
        return 1e15

    Dk = dynamic_amplification(f_id, fn, zeta, q=0)

    L = 0.0
    for ki in range(nf):
        SDk = S * Dk[ki] + Se
        if SDk <= 0:
            return 1e15

        ReEk = np.real(e_band[ki])
        tr_Ek = np.real(np.trace(ReEk))
        uAu = np.real(phi_u.T @ ReEk @ phi_u)
        fac = 1.0 / (1.0 + Se / (S * Dk[ki]))

        L += np.log(SDk)
        L += (1.0 / Se) * (tr_Ek - fac * uAu)

    L += (nr - 1) * nf * np.log(Se)
    return float(L)


def simulate_sdof(fn, zeta, S_p, Dt, N, rng):
    """Newmark-beta SDOF simulation; returns acceleration response."""
    omega_n = 2 * np.pi * fn
    c_coef = 2 * zeta * omega_n
    k_stif = omega_n**2

    beta = 0.25
    gamma = 0.50
    a1 = 1 / (beta * Dt**2)
    a2 = 1 / (beta * Dt)
    a3 = 1 / (2 * beta) - 1
    a4 = gamma / (beta * Dt)
    a5 = 1 - gamma / beta
    a6 = Dt * (1 - gamma / (2 * beta))
    k_eff = k_stif + a1 + a4 * c_coef

    sigma_p = np.sqrt(S_p / Dt)
    p = sigma_p * rng.standard_normal(N)

    x = np.zeros(N)
    xdot = np.zeros(N)
    xddot = np.zeros(N)

    for j in range(1, N):
        p_eff = (
            p[j]
            + (a1 * x[j - 1] + a2 * xdot[j - 1] + a3 * xddot[j - 1])
            + c_coef * (a4 * x[j - 1] + a5 * xdot[j - 1] + a6 * xddot[j - 1])
        )
        x[j] = p_eff / k_eff
        xddot[j] = a1 * (x[j] - x[j - 1]) - a2 * xdot[j - 1] - a3 * xddot[j - 1]
        xdot[j] = xdot[j - 1] + Dt * ((1 - gamma) * xddot[j - 1] + gamma * xddot[j])

    return xddot


def mac(phi_a, phi_b):
    """Modal Assurance Criterion."""
    phi_a = np.asarray(phi_a, dtype=float)
    phi_b = np.asarray(phi_b, dtype=float)
    return (abs(phi_a.T @ phi_b) ** 2) / ((phi_a.T @ phi_a) * (phi_b.T @ phi_b))


def main(save_figures=True, show_figures=False):
    rng = np.random.default_rng(42)

    print("==========================================================")
    print("  FULL AMBIENT VIBRATION ANALYSIS PIPELINE")
    print("  Python version of the MATLAB 5-phase pipeline")
    print("==========================================================\n")

    # ------------------------------------------------------------------
    # PHASE 1 — Planning and synthetic data generation
    # ------------------------------------------------------------------
    print("PHASE 1: Planning & Data Generation")
    print("-------------------------------------")

    n_dof = 6
    fn_true = np.array([2.0, 5.8, 9.2, 12.1, 14.5, 16.2])
    zeta_true = np.array([0.02, 0.02, 0.02, 0.02, 0.02, 0.02])

    Phi_true = np.zeros((n_dof, n_dof))
    for i in range(n_dof):
        for j in range(n_dof):
            # MATLAB: sin(j*i*pi/(2*(n+1))) with j,i = 1..n
            Phi_true[j, i] = np.sin((j + 1) * (i + 1) * np.pi / (2 * (n_dof + 1)))
        Phi_true[:, i] /= np.linalg.norm(Phi_true[:, i])

    mode_id = 1  # user-facing mode number; Python index is mode_idx
    mode_idx = mode_id - 1
    phi_true = Phi_true[:, mode_idx]
    fn_true_1 = fn_true[mode_idx]
    zeta_true_1 = zeta_true[mode_idx]

    # DOF numbers are user-facing 1..6; convert to Python indices for arrays.
    setup_dofs_user = [np.array([1, 2, 3, 4]), np.array([3, 4, 5, 6]), np.array([2, 4, 5, 6])]
    setup_dofs = [d - 1 for d in setup_dofs_user]
    n_setups = len(setup_dofs)
    ref_dof_user = 4

    fs = 200.0
    Dt = 1 / fs
    Ttot = 600.0
    N0 = int(Ttot * fs)

    S_modal = np.array([1.0, 0.6, 0.4, 0.3, 0.25, 0.2]) * 1e-6
    Se_true = 5e-9

    print(f"  Structure: {n_dof}-DOF shear frame")
    print(f"  Target: Mode {mode_id}, fn = {fn_true_1:.2f} Hz, zeta = {zeta_true_1*100:.1f}%")
    print("  Sensor setups:")
    for r, dofs in enumerate(setup_dofs_user, start=1):
        print(f"    Setup {r}: DOFs {dofs.tolist()} (ref DOF: {ref_dof_user})")
    print(f"  Acquisition: fs = {fs:g} Hz, T = {Ttot:g} s per setup")

    Y_setup = []
    for r, dofs in enumerate(setup_dofs, start=1):
        nr = len(dofs)
        Y_r = np.zeros((nr, N0))
        for i in range(n_dof):
            phi_i = Phi_true[dofs, i]
            g_i = simulate_sdof(fn_true[i], zeta_true[i], S_modal[i], Dt, N0, rng)
            Y_r += np.outer(phi_i, g_i)

        sigma_e = np.sqrt(Se_true / Dt)
        Y_r += sigma_e * rng.standard_normal((nr, N0))
        Y_setup.append(Y_r)
        print(f"    Setup {r}: {nr} channels generated")
    print("  Done.\n")

    # ------------------------------------------------------------------
    # PHASE 2 — PSD calculation
    # ------------------------------------------------------------------
    print("PHASE 2: PSD Calculation")
    print("-------------------------")

    fL = 0.5
    fU = 20.0
    Mf = 400
    df = (fU - fL) / Mf
    M = int(np.floor(N0 * Dt * df))
    Nw = int(np.floor(N0 / M))
    N = Nw * M
    f_seg = np.arange(Nw) / (Nw * Dt)

    k_band = np.where((f_seg >= fL) & (f_seg <= fU))[0]
    nk = len(k_band)
    f_band = f_seg[k_band]

    print(f"  Band: [{fL:.1f}, {fU:.1f}] Hz, Mf = {Mf}")
    print(f"  df = {df:.4f} Hz, M = {M} segments, Nw = {Nw} samples/segment ({Nw*Dt:.1f} s)")

    E_avg_setup = []
    for Y_r in Y_setup:
        Y_r = Y_r - np.mean(Y_r, axis=1, keepdims=True)
        nr = Y_r.shape[0]
        Y_trunc = Y_r[:, :N]
        E_acc = np.zeros((nk, nr, nr), dtype=complex)

        for seg in range(M):
            idx0 = seg * Nw
            idx1 = idx0 + Nw
            Y_seg = Y_trunc[:, idx0:idx1]
            F_all = np.sqrt(Dt / Nw) * np.fft.fft(Y_seg, n=Nw, axis=1)
            F_band = F_all[:, k_band]
            for ki in range(nk):
                Fk = F_band[:, ki]
                E_acc[ki] += np.outer(Fk, np.conjugate(Fk))
        E_avg_setup.append(E_acc / M)

    print(f"  PSD matrices computed for all {n_setups} setups.\n")

    rootPSD_s1 = np.sqrt(np.abs(np.real(np.array([np.diag(E_avg_setup[0][ki]) for ki in range(nk)])))).T

    plt.figure(figsize=(9, 3))
    for ch in range(rootPSD_s1.shape[0]):
        plt.semilogy(f_band, rootPSD_s1[ch], label=f"DOF {setup_dofs_user[0][ch]}")
    for f in fn_true:
        if fL <= f <= fU:
            plt.axvline(f, linestyle=":", linewidth=0.8)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("sqrt(PSD) [(m/s²)/sqrt(Hz)]")
    plt.title(f"Phase 2: Root PSD Spectrum — Setup 1, M={M} segments")
    plt.grid(True)
    plt.xlim(fL, fU)
    plt.legend()
    plt.tight_layout()
    if save_figures:
        plt.savefig("phase2_root_psd.png", dpi=300)

    # ------------------------------------------------------------------
    # PHASE 3 — Singular value spectrum
    # ------------------------------------------------------------------
    print("PHASE 3: Singular Value Spectrum")
    print("---------------------------------")

    nr1 = len(setup_dofs[0])
    SV_s1 = np.zeros((nr1, nk))
    EV_s1 = np.zeros((nk, nr1, nr1))

    for ki in range(nk):
        ReE = np.real(E_avg_setup[0][ki])
        ReE = 0.5 * (ReE + ReE.T)
        ev, V = np.linalg.eigh(ReE)
        idx_s = np.argsort(ev)[::-1]
        SV_s1[:, ki] = ev[idx_s]
        EV_s1[ki] = V[:, idx_s]

    plt.figure(figsize=(9, 3))
    for i in range(nr1):
        plt.semilogy(f_band, np.sqrt(np.abs(SV_s1[i])), label=f"SV{i+1}")
    for f in fn_true:
        if fL <= f <= fU:
            plt.axvline(f, linestyle=":", linewidth=0.8)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("sqrt(SV) [(m/s²)/sqrt(Hz)]")
    plt.title(f"Phase 3: Root SV Spectrum — Setup 1 ({nr1} channels)")
    plt.grid(True)
    plt.xlim(fL, fU)
    plt.legend()
    plt.tight_layout()
    if save_figures:
        plt.savefig("phase3_sv_spectrum.png", dpi=300)

    peaks, properties = find_peaks(np.log(np.maximum(SV_s1[0], 1e-300)), prominence=0.5, distance=round(0.5 / df))
    peak_freqs = f_band[peaks]
    print("  Resonance peaks detected (SV1):", " ".join(f"{f:.2f}" for f in peak_freqs), "Hz")
    print("  True natural frequencies:       ", " ".join(f"{f:.2f}" for f in fn_true), "Hz\n")

    # ------------------------------------------------------------------
    # PHASE 4 — Local modal identification
    # ------------------------------------------------------------------
    print("PHASE 4: Local Modal Identification")
    print("-------------------------------------")

    bw_factor = 5
    bw = bw_factor * fn_true_1 * zeta_true_1
    fL_id = fn_true_1 - bw
    fU_id = fn_true_1 + bw
    print(f"  Mode {mode_id} band: [{fL_id:.3f}, {fU_id:.3f}] Hz (±{bw_factor}× half-power BW)")

    k_id = np.where((f_band >= fL_id) & (f_band <= fU_id))[0]
    f_id = f_band[k_id]
    Nf_id = len(k_id)
    print(f"  Nf (FFT points in band) = {Nf_id}\n")

    local_fn = np.zeros(n_setups)
    local_zeta = np.zeros(n_setups)
    local_S = np.zeros(n_setups)
    local_Se = np.zeros(n_setups)
    local_phi = []

    for r in range(n_setups):
        print(f"  --- Setup {r+1} ---")
        E_band_r = E_avg_setup[r][k_id]
        nr = E_band_r.shape[1]

        A0_r = np.sum(np.real(E_band_r), axis=0)
        d_tot = sum(np.real(np.trace(E_band_r[ki])) for ki in range(Nf_id))

        ev0r, V0r = np.linalg.eigh(A0_r)
        imax = np.argmax(ev0r)
        phi0_r = np.real(V0r[:, imax])
        phi0_r = phi0_r / np.linalg.norm(phi0_r)

        SV1_band = SV_s1[0, k_id]
        kpeak = int(np.argmax(SV1_band))
        fn0_r = f_id[kpeak]

        Dk0 = dynamic_amplification(f_id, fn0_r, 0.01, q=0)
        k0_hat = phi0_r.T @ A0_r @ phi0_r
        d_hat = sum(phi0_r.T @ np.real(E_band_r[ki]) @ phi0_r for ki in range(Nf_id))
        Se0_r = max((d_tot - k0_hat) / ((nr - 1) * Nf_id), 1e-20)
        S0_r = max(np.sum(d_hat / Dk0) / Nf_id, 1e-20)
        zeta0_r = 0.01

        print(f"    Initial: fn={fn0_r:.3f} Hz, zeta={zeta0_r*100:.2f}%, S={S0_r:.3e}, Se={Se0_r:.3e}")

        theta0 = np.array([fn0_r, zeta0_r, np.log(S0_r), np.log(Se0_r)])
        result = minimize(
            nllf_reduced,
            theta0,
            args=(f_id, E_band_r, nr, Nf_id, phi0_r),
            method="Nelder-Mead",
            options={"maxfev": 5000, "maxiter": 5000, "xatol": 1e-10, "fatol": 1e-10, "disp": False},
        )

        theta_hat = result.x
        fn_hat = theta_hat[0]
        zeta_hat = abs(theta_hat[1])
        S_hat = np.exp(theta_hat[2])
        Se_hat = np.exp(theta_hat[3])

        Dk_hat = dynamic_amplification(f_id, fn_hat, zeta_hat, q=0)
        A_hat = np.zeros((nr, nr))
        for ki in range(Nf_id):
            factor_k = 1.0 / (1.0 + Se_hat / (S_hat * Dk_hat[ki]))
            A_hat += factor_k * np.real(E_band_r[ki])

        ev_hat, V_hat = np.linalg.eigh(A_hat)
        phi_hat = np.real(V_hat[:, np.argmax(ev_hat)])
        phi_hat = phi_hat / np.linalg.norm(phi_hat)

        local_fn[r] = fn_hat
        local_zeta[r] = zeta_hat
        local_S[r] = S_hat
        local_Se[r] = Se_hat
        local_phi.append(phi_hat)

        phi_true_r = phi_true[setup_dofs[r]]
        phi_true_r = phi_true_r / np.linalg.norm(phi_true_r)
        mac_r = mac(phi_true_r, phi_hat)

        print(f"    MPV:     fn={fn_hat:.3f} Hz, zeta={zeta_hat*100:.2f}%, S={S_hat:.3e}, Se={Se_hat:.3e}")
        print(f"    True:    fn={fn_true_1:.3f} Hz, zeta={zeta_true_1*100:.2f}%")
        print(f"    Local mode shape MAC = {mac_r:.4f}\n")

    print(f"  Local identification summary (Mode {mode_id}):")
    print("  Setup    fn (Hz)      zeta (%)    MAC (true)")
    for r in range(n_setups):
        phi_true_r = phi_true[setup_dofs[r]]
        phi_true_r = phi_true_r / np.linalg.norm(phi_true_r)
        mac_r = mac(phi_true_r, local_phi[r])
        print(f"  {r+1:<8d} {local_fn[r]:<12.4f} {local_zeta[r]*100:<10.3f} {mac_r:<12.4f}")
    print()

    # ------------------------------------------------------------------
    # PHASE 5 — Global mode shape assembly
    # ------------------------------------------------------------------
    print("PHASE 5: Global Mode Shape Assembly")
    print("-------------------------------------")

    L_cells = []
    for dofs in setup_dofs:
        nr = len(dofs)
        Lr = np.zeros((nr, n_dof))
        for j, dof_idx in enumerate(dofs):
            Lr[j, dof_idx] = 1.0
        L_cells.append(Lr)

    w_weights = np.ones(n_setups) / n_setups

    A0_gls = np.zeros((n_dof, n_dof))
    for r in range(n_setups):
        Lr = L_cells[r]
        nr = Lr.shape[0]
        tt = local_phi[r].reshape(-1, 1)
        A0_gls += w_weights[r] * (Lr.T @ (np.eye(nr) - tt @ tt.T) @ Lr)

    ev_gls, V_gls = np.linalg.eigh(A0_gls)
    phi_global = np.real(V_gls[:, np.argmin(ev_gls)])
    phi_global = phi_global / np.linalg.norm(phi_global)

    print("  Initial global guess:", np.array2string(phi_global, precision=4))

    max_iter = 300
    tol = 1e-12
    for it in range(1, max_iter + 1):
        phi_old = phi_global.copy()

        c_r = np.zeros(n_setups)
        lambda_r = np.zeros(n_setups)
        for r in range(n_setups):
            Lr = L_cells[r]
            tr = Lr @ phi_global
            tt = local_phi[r]
            c_r[r] = np.linalg.norm(tr) ** 2
            lambda_r[r] = w_weights[r] * (1 - abs(tt.T @ tr) / max(np.linalg.norm(tr), 1e-15))

        A_mat = np.zeros((n_dof, n_dof))
        b_vec = np.zeros(n_dof)
        for r in range(n_setups):
            Lr = L_cells[r]
            tt = local_phi[r]
            A_mat += (w_weights[r] + lambda_r[r]) * (Lr.T @ Lr)
            b_vec += -w_weights[r] * np.sqrt(c_r[r]) * (Lr.T @ tt)

        D_mat = np.block([[A_mat, np.outer(b_vec, b_vec)], [np.eye(n_dof), A_mat]])
        ev_D, V_D = np.linalg.eig(D_mat)
        idx_min = np.argmin(np.real(ev_D))
        v_full = np.real(V_D[:, idx_min])
        phi_global = v_full[:n_dof]
        phi_global = phi_global / np.linalg.norm(phi_global)

        if np.linalg.norm(phi_global - phi_old) < tol:
            break

    print(f"  GLS converged in {it} iterations.")

    if phi_true.T @ phi_global < 0:
        phi_global = -phi_global
    phi_true_norm = phi_true / np.linalg.norm(phi_true)

    mac_global = mac(phi_true_norm, phi_global)
    rmse_global = np.sqrt(np.mean((phi_global - phi_true_norm) ** 2))

    print("\n  Global mode shape (unit norm):")
    print("    Identified:", np.array2string(phi_global, precision=4))
    print("    True:      ", np.array2string(phi_true_norm, precision=4))
    print(f"    MAC  = {mac_global:.6f} (1.0 = perfect)")
    print(f"    RMSE = {rmse_global:.6f}\n")

    # Final plots
    dofs_plot_user = np.arange(1, n_dof + 1)
    plt.figure(figsize=(11, 7))

    plt.subplot(2, 3, (1, 4))
    plt.plot(np.r_[0, phi_true_norm], np.r_[0, dofs_plot_user], "--o", linewidth=2, markersize=8, label="True")
    plt.plot(np.r_[0, phi_global], np.r_[0, dofs_plot_user], "-o", linewidth=2, markersize=8, label=f"GLS (MAC={mac_global:.4f})")
    for r in range(n_setups):
        dofs_user = setup_dofs_user[r]
        phi_local_r = L_cells[r] @ phi_global
        scale_r = np.sign(phi_local_r.T @ local_phi[r]) / max(np.linalg.norm(phi_local_r), 1e-15)
        phi_local_scaled = local_phi[r] * scale_r
        plt.plot(phi_local_scaled, dofs_user, "s--", linewidth=1.2, markersize=6, label=f"Local S{r+1}")
    plt.axvline(0, linewidth=0.5)
    plt.axhline(ref_dof_user, linestyle=":", linewidth=0.8)
    plt.xlabel("Mode shape component")
    plt.ylabel("DOF (floor)")
    plt.title(f"Mode {mode_id} Global Assembly\nfn={np.mean(local_fn):.3f} Hz, MAC={mac_global:.4f}")
    plt.yticks(dofs_plot_user)
    plt.ylim(0.5, n_dof + 0.5)
    plt.grid(True)
    plt.legend()

    plt.subplot(2, 3, 2)
    for i in range(nr1):
        plt.semilogy(f_id, np.sqrt(np.abs(SV_s1[i, k_id])), linewidth=1.2, label=f"SV{i+1}")
    plt.axvline(fn_true_1, linestyle=":", linewidth=1, label="fn true")
    plt.axvline(np.mean(local_fn), linestyle="--", linewidth=1, label="fn MPV")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("sqrt(SV)")
    plt.title("SV Spectrum (Setup 1, ID band)")
    plt.grid(True)
    plt.xlim(fL_id, fU_id)
    plt.legend()

    plt.subplot(2, 3, 3)
    for r in range(n_setups):
        dofs_user = setup_dofs_user[r]
        phi_local_r = L_cells[r] @ phi_global
        scale_r = np.sign(phi_local_r.T @ local_phi[r]) / max(np.linalg.norm(phi_local_r), 1e-15)
        phi_local_scaled = local_phi[r] * scale_r
        plt.plot(phi_local_scaled, dofs_user, "s-", linewidth=1.5, markersize=7, label=f"Local S{r+1}")
    plt.plot(phi_true_norm, dofs_plot_user, "--", linewidth=2, label="True (full)")
    plt.axvline(0, linewidth=0.5)
    plt.xlabel("Mode shape component")
    plt.ylabel("DOF")
    plt.title("Local Mode Shapes vs. True")
    plt.yticks(dofs_plot_user)
    plt.ylim(0.5, n_dof + 0.5)
    plt.grid(True)
    plt.legend()

    plt.subplot(2, 3, 5)
    x = np.arange(n_setups)
    width = 0.35
    plt.bar(x - width / 2, local_fn, width, label="Identified")
    plt.bar(x + width / 2, np.full(n_setups, fn_true_1), width, label="True")
    plt.xticks(x, [f"S{i+1}" for i in range(n_setups)])
    plt.ylabel("Natural Frequency (Hz)")
    plt.title("fn: Identified vs. True")
    plt.grid(True)
    plt.legend()

    plt.subplot(2, 3, 6)
    plt.bar(x - width / 2, local_zeta * 100, width, label="Identified")
    plt.bar(x + width / 2, np.full(n_setups, zeta_true_1 * 100), width, label="True")
    plt.xticks(x, [f"S{i+1}" for i in range(n_setups)])
    plt.ylabel("Damping Ratio (%)")
    plt.title("zeta: Identified vs. True")
    plt.grid(True)
    plt.legend()

    plt.suptitle(f"5-Phase Ambient Vibration Analysis — Mode {mode_id} (fn={fn_true_1:.2f} Hz)")
    plt.tight_layout()
    if save_figures:
        plt.savefig("phase5_global_mode_shape.png", dpi=300)

    print("==========================================================")
    print(f"  FINAL RESULTS SUMMARY — MODE {mode_id}")
    print("==========================================================")
    print(f"  True fn        = {fn_true_1:.4f} Hz")
    print(f"  Mean MPV fn    = {np.mean(local_fn):.4f} Hz (error: {abs(np.mean(local_fn)-fn_true_1)/fn_true_1*100:.2f}%)")
    print(f"  True zeta      = {zeta_true_1*100:.2f}%")
    print(f"  Mean MPV zeta  = {np.mean(local_zeta)*100:.2f}% (error: {abs(np.mean(local_zeta)-zeta_true_1)/zeta_true_1*100:.2f}%)")
    print(f"  Global MAC     = {mac_global:.6f}")
    print(f"  Global RMSE    = {rmse_global:.6f}")
    print("==========================================================")
    print("Done.")

    if show_figures:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main(save_figures=True, show_figures=False)

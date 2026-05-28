import argparse
import os
from datetime import datetime

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from arch import arch_model
from statsmodels.stats.diagnostic import acorr_ljungbox, het_arch

GARCH_SCALE_FACTOR = 100.0


def parse_args():
    parser = argparse.ArgumentParser(description="Monte Carlo Simulation with GARCH")
    parser.add_argument("ticker", help="Ticker symbol")
    parser.add_argument("--data-dir", type=str, default="data",
                        help="Data directory (default: data)")
    parser.add_argument("--start-price", type=float, default=None,
                        help="Starting price for simulation (default: last close)")
    parser.add_argument("--days", type=int, default=20,
                        help="Simulation days (default: 20)")
    parser.add_argument("--paths", type=int, default=5000,
                        help="Number of simulation paths (default: 5000)")
    parser.add_argument("--target", type=float, default=None,
                        help="Target price for probability calc (default: start + 7%%)")
    parser.add_argument("--stop-loss", type=str, default="2xatr",
                        help="Stop-loss: ATR-multiple (e.g. '2xatr'), percentage (e.g. 0.03 = 3%%) or absolute price (e.g. 379). Default: 2xatr")
    parser.add_argument("--trailing-stop", type=str, default="2xatr",
                        help="Trailing stop: ATR-multiple (e.g. '2xatr') or percentage (e.g. 0.03). 0 to disable. Default: 2xatr")
    parser.add_argument("--no-trailing", action="store_true",
                        help="Disable trailing stop")
    parser.add_argument("--short", action="store_true",
                        help="Simulate a short position (inverted stops and target)")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for reproducibility (default: none)")
    parser.add_argument("--lookback", type=int, default=None,
                        help="Antal senaste handelsdagar att fitta GARCH på (default: hela historiken)")
    return parser.parse_args()


def load_data(ticker, data_dir, lookback=None):
    if not ticker.endswith(".csv"):
        csv_path = os.path.join(data_dir, f"{ticker}.csv")
    else:
        csv_path = os.path.join(data_dir, ticker)
        ticker = ticker.replace(".csv", "")

    df = pd.read_csv(csv_path)
    df["Date"] = pd.to_datetime(
        df["Date"].astype(str).str.split(" ").str[0], errors="coerce"
    )
    df = df.dropna(subset=["Date", "Close"]).sort_values("Date").set_index("Date")

    # Log returns: ln(P_t / P_{t-1})
    returns = np.log(df["Close"] / df["Close"].shift(1)).dropna()
    if lookback is not None and lookback > 0:
        returns = returns.tail(lookback)

    last_close = float(df["Close"].iloc[-1])
    last_date = df.index[-1]

    # ATR(14) as fraction of last close — used for ATR-multiple stop sizing
    atr_pct = None
    if {"High", "Low", "Close"}.issubset(df.columns) and len(df) >= 15:
        high = df["High"].astype(float)
        low = df["Low"].astype(float)
        close = df["Close"].astype(float)
        prev_close = close.shift(1)
        tr = pd.concat([
            high - low,
            (high - prev_close).abs(),
            (low - prev_close).abs(),
        ], axis=1).max(axis=1)
        atr = tr.rolling(14).mean().iloc[-1]
        if pd.notna(atr) and last_close > 0:
            atr_pct = float(atr / last_close)

    if len(returns) < 30:
        raise ValueError(f"Bara {len(returns)} datapunkter — minst 30 krävs för GARCH-fit")
    if len(returns) < 60:
        print(f"  Varning: Bara {len(returns)} datapunkter — GARCH-fit kan bli opålitlig")

    return ticker, returns, last_close, last_date, atr_pct


def fit_garch(returns):
    scaled = returns * GARCH_SCALE_FACTOR
    garch = arch_model(scaled, vol="Garch", p=1, o=1, q=1, mean="Constant", dist="t")
    fit = garch.fit(disp="off")

    try:
        std_resid = fit.std_resid.dropna().values
        if len(std_resid) < 10:
            std_resid = np.random.normal(0, 1, 1000)
    except Exception:
        std_resid = np.random.normal(0, 1, 1000)

    try:
        params = fit.params
        omega = float(params.get("omega", 1e-6)) / (GARCH_SCALE_FACTOR ** 2)
        alpha = float(params.get("alpha[1]", 0.1))
        gamma = float(params.get("gamma[1]", 0.0))
        beta = float(params.get("beta[1]", 0.8))
        mu = float(params.get("mu", 0.0)) / GARCH_SCALE_FACTOR
    except Exception:
        omega, alpha, gamma, beta, mu = 1e-6, 0.1, 0.05, 0.8, 0.0

    try:
        last_var = float(fit.conditional_volatility.iloc[-1] ** 2) / (GARCH_SCALE_FACTOR ** 2)
    except Exception:
        last_var = 1e-6

    try:
        last_ret = float(returns.iloc[-1])
    except Exception:
        last_ret = 0.0

    return omega, alpha, gamma, beta, mu, std_resid, last_var, last_ret


def diagnose_garch(std_resid, alpha, gamma, beta):
    persistence = alpha + gamma / 2 + beta

    lb_result = acorr_ljungbox(std_resid**2, lags=[10, 20], return_df=True)
    lb_pvalue = float(lb_result["lb_pvalue"].iloc[-1])

    arch_stat, arch_pvalue, _, _ = het_arch(std_resid, nlags=10)
    arch_pvalue = float(arch_pvalue)

    warnings = []
    if persistence >= 1.0:
        warnings.append(f"Persistence α+γ/2+β={persistence:.3f} ≥ 1 — icke-stationär variansprocess")
    if lb_pvalue < 0.05:
        warnings.append(f"Ljung-Box p={lb_pvalue:.3f} — kvarvarande ARCH-effekter")
    if arch_pvalue < 0.05:
        warnings.append(f"ARCH-LM p={arch_pvalue:.3f} — kvarvarande heteroskedasticitet")

    return {
        "passed": len(warnings) == 0,
        "warnings": warnings,
        "persistence": persistence,
        "lb_pvalue": lb_pvalue,
        "arch_pvalue": arch_pvalue,
    }


def simulate(start_price, num_paths, sim_days, stop_loss_pct, trailing_stop,
             trailing_stop_pct, omega, alpha, gamma, beta, mu, std_resid, last_var, last_ret,
             short=False):
    prices = np.zeros((sim_days + 1, num_paths))
    prices[0, :] = start_price

    stopped = np.zeros(num_paths, dtype=bool)
    stop_day = np.full(num_paths, -1, dtype=int)
    stop_type = np.zeros(num_paths, dtype=int)  # 0=not stopped, 1=trailing, 2=fixed
    current_extreme = np.full(num_paths, start_price)
    trail_levels = np.zeros_like(prices)

    if short:
        stop_price_fixed = start_price * (1 + stop_loss_pct)
        trail_levels[0, :] = start_price * (1 + trailing_stop_pct) if trailing_stop else np.nan
    else:
        stop_price_fixed = start_price * (1 - stop_loss_pct)
        trail_levels[0, :] = start_price * (1 - trailing_stop_pct) if trailing_stop else np.nan

    var = np.full(num_paths, last_var)
    ret = np.full(num_paths, last_ret)

    for t in range(1, sim_days + 1):
        active = ~stopped
        if not active.any():
            prices[t] = prices[t - 1]
            trail_levels[t] = trail_levels[t - 1]
            continue

        shock = (ret[active] - mu) ** 2
        leverage = np.where(ret[active] - mu < 0, 1.0, 0.0)
        var[active] = omega + (alpha + gamma * leverage) * shock + beta * var[active]
        var = np.maximum(var, 1e-9)

        z = np.random.choice(std_resid, size=num_paths)
        daily_log_ret = mu + np.sqrt(var) * z

        # Apply log returns via exp: P_t = P_{t-1} * exp(r_t)
        prices[t, active] = prices[t - 1, active] * np.exp(daily_log_ret[active])
        prices[t, ~active] = prices[t - 1, ~active]

        if short:
            current_extreme[active] = np.minimum(current_extreme[active], prices[t, active])
        else:
            current_extreme[active] = np.maximum(current_extreme[active], prices[t, active])

        if trailing_stop:
            if short:
                trail_level = current_extreme * (1 + trailing_stop_pct)
            else:
                trail_level = current_extreme * (1 - trailing_stop_pct)
            trail_levels[t, active] = trail_level[active]
            trail_levels[t, ~active] = trail_levels[t - 1, ~active]

            if short:
                trail_hit = active & (prices[t] >= trail_level)
            else:
                trail_hit = active & (prices[t] <= trail_level)
            if trail_hit.any():
                trail_levels[t, trail_hit] = trail_level[trail_hit]
                stopped[trail_hit] = True
                stop_day[trail_hit] = t
                # Om priset aldrig rört sig gynnsamt = fixed stop-loss, annars trailing
                if short:
                    favorable = current_extreme[trail_hit] < start_price
                else:
                    favorable = current_extreme[trail_hit] > start_price
                stop_type[trail_hit] = np.where(favorable, 1, 2)
                active = active & ~trail_hit
        else:
            trail_levels[t] = np.nan

        if short:
            fixed_hit = active & (prices[t] >= stop_price_fixed)
        else:
            fixed_hit = active & (prices[t] <= stop_price_fixed)
        if fixed_hit.any():
            if trailing_stop:
                trail_levels[t, fixed_hit] = stop_price_fixed
            stopped[fixed_hit] = True
            stop_day[fixed_hit] = t
            stop_type[fixed_hit] = 2

        ret[~stopped] = daily_log_ret[~stopped]

    return prices, stopped, stop_day, stop_type, trail_levels, stop_price_fixed


def calc_statistics(prices, start_price, target_price, num_paths, stopped, stop_day, stop_type, short=False):
    final_prices = np.array(prices[-1, :], dtype=np.float64)
    final_prices = np.nan_to_num(final_prices, nan=start_price, posinf=start_price, neginf=start_price)

    pnl = (start_price - final_prices) if short else (final_prices - start_price)

    win_mask = pnl > 0
    loss_mask = pnl < 0
    num_wins = int(np.sum(win_mask))
    num_losses = int(np.sum(loss_mask))

    win_rate = float(num_wins / num_paths)
    loss_rate = float(num_losses / num_paths)

    avg_win = float(np.median(pnl[win_mask])) if num_wins > 0 else 0.0
    avg_loss = float(np.median(pnl[loss_mask])) if num_losses > 0 else 0.0

    payoff_ratio = None
    be_win_rate = 0.0

    if num_losses > 0 and avg_loss != 0:
        payoff_ratio = float(abs(avg_win / avg_loss))
        be_win_rate = (1 / (1 + payoff_ratio)) * 100

    # EV = medel-PnL över alla paths. Medel (inte median som payoff) eftersom
    # de få stora vinsterna i högersvansen är just det som ger edgen för
    # trendföljande strategier — medianen klipper bort dem.
    ev = float(np.mean(pnl))
    ev_pct = float((ev / start_price) * 100) if start_price != 0 else 0.0

    if payoff_ratio and payoff_ratio > 0 and avg_win > 0:
        kelly = max((win_rate * payoff_ratio - loss_rate) / payoff_ratio, 0.0)
    else:
        kelly = 0.0
    half_kelly = float(kelly * 50)

    stop_rate = float(np.mean(stopped) * 100)
    avg_stop_day = float(np.mean(stop_day[stopped])) if np.any(stopped) else 0.0

    trail_mask = stop_type == 1
    fixed_mask = stop_type == 2
    trail_rate = float(trail_mask.sum() / num_paths * 100)
    fixed_rate = float(fixed_mask.sum() / num_paths * 100)
    trail_win_rate = float((trail_mask & win_mask).sum() / trail_mask.sum() * 100) if trail_mask.any() else 0.0
    trail_loss_rate = float((trail_mask & loss_mask).sum() / trail_mask.sum() * 100) if trail_mask.any() else 0.0

    if short:
        prob_target = float((final_prices < target_price).mean() * 100)
    else:
        prob_target = float((final_prices > target_price).mean() * 100)

    median_win = float(np.median(final_prices[win_mask])) if num_wins > 0 else start_price
    median_loss = float(np.median(final_prices[loss_mask])) if num_losses > 0 else start_price

    return {
        "final_prices": final_prices,
        "win_rate": win_rate,
        "loss_rate": loss_rate,
        "avg_win": avg_win,
        "avg_loss": avg_loss,
        "payoff_ratio": payoff_ratio,
        "be_win_rate": be_win_rate,
        "ev": ev,
        "ev_pct": ev_pct,
        "kelly": kelly,
        "half_kelly": half_kelly,
        "stop_rate": stop_rate,
        "avg_stop_day": avg_stop_day,
        "trail_rate": trail_rate,
        "fixed_rate": fixed_rate,
        "trail_win_rate": trail_win_rate,
        "trail_loss_rate": trail_loss_rate,
        "prob_target": prob_target,
        "median_win": median_win,
        "median_loss": median_loss,
    }


def save_report(ticker, stock_name, last_close, start_price, target_price,
                sim_days, num_paths, stop_price_fixed, stop_loss_pct,
                trailing_stop, trailing_stop_pct, stats, diagnostics, fit_window, lookback_used, short=False):
    output_file = f"montecarlo_{ticker}.txt"
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("=" * 80 + "\n")
        f.write("MONTE CARLO SIMULATION (BOOTSTRAPPED GJR-GARCH)\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Rapport skapad: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"Ticker: {stock_name}\n")
        f.write(f"Position: {'SHORT' if short else 'LONG'}\n")
        f.write(f"Method: Filtered Historical Simulation (GJR-GARCH, Student-t, Log Returns)\n")
        f.write(f"Last Close: {last_close:.2f}\n")
        f.write(f"Start Price: {start_price:.2f}\n")
        f.write(f"Target Price: {target_price:.2f}\n")
        f.write(f"Simulation Days: {sim_days}\n")
        f.write(f"Paths: {num_paths:,}\n")
        f.write(f"Stop-Loss: {stop_price_fixed:.2f} ({stop_loss_pct:.0%})\n")
        f.write(f"Trailing Stop: {'Yes' if trailing_stop else 'No'} ({trailing_stop_pct:.0%})\n")
        fit_label = f"{fit_window} dagar" + (" (--lookback)" if lookback_used else " (full historik)")
        f.write(f"GARCH Fit Window: {fit_label}\n\n")

        f.write("-" * 80 + "\n")
        f.write("STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Mean Final Price: {stats['final_prices'].mean():.2f}\n")
        f.write(f"Median Exit (winners): {stats['median_win']:.2f}\n")
        f.write(f"Median Exit (losers):  {stats['median_loss']:.2f}\n")
        target_op = "<" if short else ">"
        f.write(f"Probability {target_op} Target ({target_price:.0f}): {stats['prob_target']:.2f}%\n")
        f.write(f"Win Rate: {stats['win_rate'] * 100:.1f}%\n")
        pr = stats['payoff_ratio']
        f.write(f"Payoff Ratio: {f'{pr:.2f}' if pr is not None else 'N/A'} (Break-even win rate: {stats['be_win_rate']:.1f}%)\n")
        f.write(f"Expected Value (EV): {stats['ev']:+.2f} ({stats['ev_pct']:+.2f}%)\n")
        f.write(f"Half Kelly: {stats['half_kelly']:.2f}%\n")
        f.write(f"Stopped paths: {stats['stop_rate']:.1f}% (avg day {stats['avg_stop_day']:.1f})\n")
        f.write(f"  Fixed stop-loss: {stats['fixed_rate']:.1f}% (förlust)\n")
        if stats['trail_rate'] > 0:
            f.write(f"  Trailing stop:   {stats['trail_rate']:.1f}% ({stats['trail_win_rate']:.1f}% vinst, {stats['trail_loss_rate']:.1f}% förlust)\n")
        f.write("\n")

        ok = "✓" if diagnostics["passed"] else "✗"
        lb_ok = "✓" if diagnostics["lb_pvalue"] >= 0.05 else "✗"
        arch_ok = "✓" if diagnostics["arch_pvalue"] >= 0.05 else "✗"
        pers_ok = "✓" if diagnostics["persistence"] < 1.0 else "✗"

        f.write("-" * 80 + "\n")
        f.write("MODEL DIAGNOSTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Persistence (α+γ/2+β):  {diagnostics['persistence']:.3f}  {pers_ok}\n")
        f.write(f"Ljung-Box(20) p:    {diagnostics['lb_pvalue']:.3f}  {lb_ok}\n")
        f.write(f"ARCH-LM(10) p:      {diagnostics['arch_pvalue']:.3f}  {arch_ok}\n")
        f.write(f"Model status:       {'OK' if diagnostics['passed'] else 'VARNING — se detaljer'}\n")
        if diagnostics["warnings"]:
            for w in diagnostics["warnings"]:
                f.write(f"  ⚠ {w}\n")
        f.write("\n")
        f.write("=" * 80 + "\n")

    return output_file


def plot_results(ticker, stock_name, prices, trail_levels, stats, stop_day, stop_type,
                 start_price, target_price, stop_price_fixed,
                 sim_days, num_paths, last_date,
                 trailing_stop, trailing_stop_pct, fit_window, lookback_used, short=False):
    date_range = pd.bdate_range(start=last_date, periods=sim_days + 1)

    # Active mask: path i is "alive" at day t if not stopped or stopped after t.
    # Cone/trail percentiles use only active paths to avoid frozen stop-levels
    # biasing the distribution.
    t_idx = np.arange(sim_days + 1)[:, None]
    sd = stop_day[None, :]
    active_mask = (sd == -1) | (t_idx <= sd)

    masked_prices = np.where(active_mask, prices, np.nan)
    with np.errstate(all="ignore"):
        p95 = np.nanpercentile(masked_prices, 95, axis=1)
        p75 = np.nanpercentile(masked_prices, 75, axis=1)
        p50 = np.nanpercentile(masked_prices, 50, axis=1)
        p25 = np.nanpercentile(masked_prices, 25, axis=1)
        p05 = np.nanpercentile(masked_prices, 5, axis=1)
        mean_price = np.nanmean(masked_prices, axis=1)

    active_frac = active_mask.mean(axis=1)
    # Blank cone only when very few paths remain (need ~30 for stable percentiles).
    min_paths = max(30, int(num_paths * 0.005))
    sparse = active_mask.sum(axis=1) < min_paths
    for arr in (p95, p75, p50, p25, p05, mean_price):
        arr[sparse] = np.nan

    if trailing_stop:
        masked_trail = np.where(active_mask, trail_levels, np.nan)
        with np.errstate(all="ignore"):
            median_trail = np.nanmedian(masked_trail, axis=1)
        median_trail[sparse] = np.nan
        if short:
            median_trail[0] = start_price * (1 + trailing_stop_pct)
        else:
            median_trail[0] = start_price * (1 - trailing_stop_pct)
    else:
        median_trail = np.full(len(date_range), np.nan)

    COLORS = {
        "mean": "#1565C0",
        "median": "#6A1B9A",
        "start": "#2E7D32",
        "target": "#E65100",
        "stop": "#C62828",
        "trail": "#AD1457",
        "conf_90": "#BBDEFB",
        "conf_50": "#C8E6C9",
        "hist_bull": "#4CAF50",
        "hist_bear": "#EF5350",
    }

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 12,
        "axes.labelsize": 13,
        "axes.titlesize": 18,
        "axes.titleweight": "bold",
        "axes.grid": True,
        "grid.alpha": 0.4,
        "grid.linestyle": "-",
        "grid.linewidth": 0.5,
        "lines.linewidth": 2.5,
        "legend.fontsize": 11,
        "legend.framealpha": 0.95,
        "figure.figsize": (18, 9),
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
    })

    fig = plt.figure(figsize=(18, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=[3, 1], height_ratios=[4, 1],
                          wspace=0.02, hspace=0.08)
    ax_main = fig.add_subplot(gs[0, 0])
    ax_hist = fig.add_subplot(gs[0, 1])
    ax_stop = fig.add_subplot(gs[1, 0], sharex=ax_main)

    ax_main.fill_between(date_range, p05, p95, color=COLORS["conf_90"], alpha=0.6, label="90% Konfidensintervall")
    ax_main.fill_between(date_range, p25, p75, color=COLORS["conf_50"], alpha=0.7, label="50% Konfidensintervall")

    # Sample-banor: visa upp till 15 stoppade + 15 överlevande som tunna linjer.
    # Stoppade banor klipps vid stop_day så det syns var/när de träffar SL.
    stopped_idx = np.where(stop_day != -1)[0]
    survived_idx = np.where(stop_day == -1)[0]
    rng = np.random.default_rng(0)
    n_stop_show = min(15, len(stopped_idx))
    n_surv_show = min(15, len(survived_idx))
    if n_stop_show > 0:
        sel = rng.choice(stopped_idx, size=n_stop_show, replace=False)
        for i in sel:
            d = stop_day[i]
            ax_main.plot(date_range[: d + 1], prices[: d + 1, i],
                         color="#616161", linewidth=0.6, alpha=0.35, zorder=1)
    if n_surv_show > 0:
        sel = rng.choice(survived_idx, size=n_surv_show, replace=False)
        for i in sel:
            ax_main.plot(date_range, prices[:, i],
                         color="#616161", linewidth=0.6, alpha=0.35, zorder=1)

    ax_main.plot(date_range, mean_price, color=COLORS["mean"], label="Medelvärde", linewidth=3, linestyle="--")
    ax_main.plot(date_range, p50, color=COLORS["median"], label="Median", linewidth=3)

    ax_main.axhline(start_price, color=COLORS["start"], linewidth=2.5, linestyle="-", label=f"Start: {start_price:.0f}")
    ax_main.axhline(target_price, color=COLORS["target"], linewidth=2.5, linestyle="--", label=f"Mål: {target_price:.0f}")
    ax_main.axhline(stop_price_fixed, color=COLORS["stop"], linewidth=2.5, linestyle=":", label=f"Stop-loss: {stop_price_fixed:.0f}")

    if trailing_stop:
        ax_main.plot(date_range, median_trail, color=COLORS["trail"], linewidth=2, alpha=0.9, linestyle="-.", label="Trailing stop (median)")

    ax_main.scatter([date_range[0]], [start_price], color=COLORS["start"], s=100, zorder=5, edgecolors="white", linewidths=2)

    ax_main.set_title(f"{stock_name} — {sim_days}-dagars Monte Carlo (GARCH, Log Returns)\n{num_paths:,} simuleringar", pad=15, fontsize=16)
    ax_main.set_ylabel("Pris", fontsize=13, fontweight="bold")
    ax_main.set_xlabel("Datum", fontsize=13, fontweight="bold")

    ax_main.xaxis.set_major_locator(mdates.DayLocator(interval=3))
    ax_main.xaxis.set_minor_locator(mdates.DayLocator())
    ax_main.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    fig.autofmt_xdate(rotation=45, ha="right")

    ax_main.yaxis.set_major_locator(mticker.MaxNLocator(nbins=12))
    ax_main.yaxis.set_major_formatter(mticker.StrMethodFormatter("{x:,.0f}"))

    y_min_data = min(np.nanmin(p05), stop_price_fixed, target_price) * 0.98
    y_max_limit = max(np.nanmax(p95), stop_price_fixed, target_price) * 1.08
    ax_main.set_ylim(y_min_data, y_max_limit)

    # Knockout-panel: kumulativ stop-out per dag, uppdelat på fixed/trailing.
    t_days = np.arange(sim_days + 1)
    fixed_per_day = np.zeros(sim_days + 1)
    trail_per_day = np.zeros(sim_days + 1)
    fixed_days = stop_day[stop_type == 2]
    trail_days = stop_day[stop_type == 1]
    for d in fixed_days:
        if 0 <= d <= sim_days:
            fixed_per_day[d] += 1
    for d in trail_days:
        if 0 <= d <= sim_days:
            trail_per_day[d] += 1
    cum_fixed = np.cumsum(fixed_per_day) / num_paths * 100
    cum_trail = np.cumsum(trail_per_day) / num_paths * 100

    ax_stop.fill_between(date_range, 0, cum_fixed,
                         color=COLORS["stop"], alpha=0.7, linewidth=0,
                         label=f"Fixed SL ({cum_fixed[-1]:.0f}%)")
    if trailing_stop and cum_trail[-1] > 0:
        ax_stop.fill_between(date_range, cum_fixed, cum_fixed + cum_trail,
                             color="#FB8C00", alpha=0.7, linewidth=0,
                             label=f"Trailing ({cum_trail[-1]:.0f}%)")
    total_stop = cum_fixed + cum_trail
    ax_stop.plot(date_range, total_stop, color="#424242", linewidth=1.2)
    ax_stop.set_ylim(0, max(100, total_stop.max() * 1.05) if total_stop.max() > 0 else 100)
    ax_stop.set_ylabel("Stoppade %", fontsize=11, fontweight="bold")
    ax_stop.set_xlabel("Datum", fontsize=13, fontweight="bold")
    ax_stop.yaxis.set_major_formatter(mticker.PercentFormatter(decimals=0))
    ax_stop.legend(loc="upper left", fontsize=10, framealpha=0.9)
    ax_stop.grid(True, alpha=0.4)

    ax_main.set_xlabel("")
    plt.setp(ax_main.get_xticklabels(), visible=False)

    ax_main.legend(loc="upper right", ncol=2, fontsize=10, framealpha=0.9)

    final_prices = stats["final_prices"]
    target_label = f"P(≤ {target_price:.0f})" if short else f"P(≥ {target_price:.0f})"
    fit_label = f"{fit_window}d" + ("*" if lookback_used else "")
    stats_lines = [
        f"  Fit window:   {fit_label}",
        f"{'─' * 22}",
        f"  EV:           {stats['ev']:+.2f} ({stats['ev_pct']:+.2f}%)",
        f"  Win Rate:     {stats['win_rate'] * 100:.1f}%",
        f"  Payoff Ratio: {stats['payoff_ratio']:.2f}" if stats['payoff_ratio'] is not None else "  Payoff Ratio: N/A",
        f"  Break-even:   {stats['be_win_rate']:.1f}%",
        f"  Half Kelly:   {stats['half_kelly']:.1f}%",
        f"{'─' * 22}",
        f"  {target_label}:    {stats['prob_target']:.1f}%",
        f"  Stopped:      {stats['stop_rate']:.1f}%",
        f"    Fixed SL:   {stats['fixed_rate']:.1f}%",
        *(
            [f"    Trailing:   {stats['trail_rate']:.1f}% ({stats['trail_win_rate']:.0f}%V/{stats['trail_loss_rate']:.0f}%F)"]
            if stats['trail_rate'] > 0 else []
        ),
    ]
    stats_text = "\n".join(stats_lines)
    props = dict(boxstyle="square,pad=0.5", facecolor="white", edgecolor="0.8", linewidth=0.8, alpha=0.9)
    ax_main.text(0.02, 0.98, stats_text, transform=ax_main.transAxes, fontsize=10,
                 verticalalignment="top", horizontalalignment="left", bbox=props)

    bins = np.linspace(y_min_data, y_max_limit, 51)
    if short:
        bull_mask = final_prices <= start_price
        bear_mask = final_prices > start_price
    else:
        bull_mask = final_prices >= start_price
        bear_mask = final_prices < start_price

    # Klippa extrema svansar till bin-rangen så de hamnar i kant-bins istället för att försvinna
    bull_clipped = np.clip(final_prices[bull_mask], y_min_data, y_max_limit)
    bear_clipped = np.clip(final_prices[bear_mask], y_min_data, y_max_limit)

    ax_hist.hist(bull_clipped, bins=bins, orientation="horizontal",
                 color=COLORS["hist_bull"], alpha=0.8, edgecolor="white", linewidth=0.5, label="Vinst")
    ax_hist.hist(bear_clipped, bins=bins, orientation="horizontal",
                 color=COLORS["hist_bear"], alpha=0.8, edgecolor="white", linewidth=0.5, label="Förlust")

    ax_hist.axhline(start_price, color=COLORS["start"], linewidth=2, linestyle="-")
    ax_hist.axhline(target_price, color=COLORS["target"], linewidth=2, linestyle="--")
    ax_hist.axhline(stop_price_fixed, color=COLORS["stop"], linewidth=2, linestyle=":")

    ax_hist.set_xlabel("Antal simuleringar", fontsize=11, fontweight="bold")
    ax_hist.set_title(f"Slutprisfördelning\n(Dag {sim_days})", fontsize=12, fontweight="bold")
    ax_hist.tick_params(axis="y", labelleft=False)
    ax_hist.legend(loc="upper right", fontsize=10)
    ax_hist.set_ylim(y_min_data, y_max_limit)

    for ax in [ax_main, ax_hist]:
        for spine in ax.spines.values():
            spine.set_color("#333333")
            spine.set_linewidth(1)

    graph_file = f"montecarlo_{ticker}.png"
    plt.savefig(graph_file, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)

    return graph_file


def main():
    args = parse_args()

    ticker = args.ticker
    num_paths = args.paths
    sim_days = args.days

    try:
        ticker, returns, last_close, last_date, atr_pct = load_data(ticker, args.data_dir, args.lookback)
    except FileNotFoundError:
        print(f"Error: File not found for {ticker}")
        return
    except Exception as e:
        print(f"Data error: {e}")
        return

    short = args.short
    stock_name = ticker.upper()
    start_price = args.start_price if args.start_price is not None else last_close
    default_target = start_price * 0.93 if short else start_price * 1.07
    target_price = args.target if args.target else default_target

    def parse_stop(value, atr_pct, start_price, short, kind):
        """Parse stop arg: 'NxATR' multiplier, float percent/price, or None."""
        s = str(value).strip().lower()
        if s.endswith("xatr"):
            if atr_pct is None:
                print(f"  Varning: ATR saknas (OHLC?) — {kind} faller tillbaka till 3%")
                return 0.03
            mult = float(s[:-4]) if s[:-4] else 1.0
            return mult * atr_pct
        val = float(s)
        if val > 1:  # Absolute price
            if short:
                pct = (val - start_price) / start_price
            else:
                pct = (start_price - val) / start_price
            return max(pct, 0.001)
        return val

    stop_loss_pct = parse_stop(args.stop_loss, atr_pct, start_price, short, "stop-loss")
    trailing_raw = parse_stop(args.trailing_stop, atr_pct, start_price, short, "trailing")
    trailing_stop = not args.no_trailing and trailing_raw > 0
    trailing_stop_pct = trailing_raw if trailing_stop else 0.0

    print(f"Analyzing {stock_name}...")
    print(f"Last close: {last_close:.2f}, Start price: {start_price:.2f}, Target: {target_price:.2f}")

    omega, alpha, gamma, beta, mu, std_resid, last_var, last_ret = fit_garch(returns)

    diagnostics = diagnose_garch(std_resid, alpha, gamma, beta)
    if diagnostics["warnings"]:
        print("  GARCH diagnostik:")
        for w in diagnostics["warnings"]:
            print(f"    ⚠ {w}")
    else:
        print("  GARCH diagnostik: OK")

    if args.seed is not None:
        np.random.seed(args.seed)

    prices, stopped, stop_day, stop_type, trail_levels, stop_price_fixed = simulate(
        start_price, num_paths, sim_days, stop_loss_pct, trailing_stop,
        trailing_stop_pct, omega, alpha, gamma, beta, mu, std_resid, last_var, last_ret,
        short=short)

    stats = calc_statistics(prices, start_price, target_price, num_paths, stopped, stop_day, stop_type, short=short)

    fit_window = len(returns)
    lookback_used = args.lookback is not None and args.lookback > 0

    output_file = save_report(
        ticker, stock_name, last_close, start_price, target_price,
        sim_days, num_paths, stop_price_fixed, stop_loss_pct,
        trailing_stop, trailing_stop_pct, stats, diagnostics,
        fit_window, lookback_used, short=short)
    print(f"Report saved: {output_file}")

    graph_file = plot_results(
        ticker, stock_name, prices, trail_levels, stats, stop_day, stop_type,
        start_price, target_price, stop_price_fixed,
        sim_days, num_paths, last_date,
        trailing_stop, trailing_stop_pct, fit_window, lookback_used, short=short)
    print(f"Graph saved: {graph_file}")


if __name__ == "__main__":
    main()

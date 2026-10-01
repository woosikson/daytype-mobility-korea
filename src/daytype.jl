# src/daytype.jl — holds only what changes when the workday pipeline is carried over to **weekend · holiday**.
#   Matrix convention as in `src/model.jl`: A[j,i], row = origin j (night-time, residence), column = destination i (daytime).
#
# ── Definitions of the three day types ────────────────
#   weekday  — Mon–Fri days that are not public holidays
#   weekend  — runs of 1–2 consecutive non-working days (Sat·Sun, **and isolated midweek public holidays**)
#   holiday  — runs of 3 or more consecutive non-working days (Lunar New Year, Chuseok, etc.)
#   ⚠️ 2022-09 has no isolated midweek public holiday, so weekend ≡ the 6 Sat·Sun days · holiday ≡ Chuseok 9/9–12.
#      In other months, midweek public holidays such as 2022-03-01 (Tue) are classified as weekend.
#
# ── Only four things change from the workday version ───────────────────────────────────
#   ① routine seed   : scale **inter**-district commutes by λ and move the removed share to the diagonal.
#   ② routine seed   : then subtract `φ·free` from that diagonal (`routine_axis_daytype` in this file).
#   ③ non-commute    : scale the participation rate by ω and add what was subtracted in ② to the same row.
#   ④ activity rate  : scale the commute budget by λ (f^{dt} = λκ — supplied by the caller).
#   The IPF and axis-split equations are unchanged — `src/model.jl` is used as is.

# ═════════════════════════════════════════════════════════════════════════════
# φ — does someone who stops commuting necessarily become non-move?
# ═════════════════════════════════════════════════════════════════════════════
# The seed is built by **adding** the non-routine axis (row sum `π P`) to the routine axis (row sum `P`), and IPF
# brings each row back to `P`. So the realized NC share is almost exactly
#
#     NC/P ≈ π/(1+π)  <  0.5      (even as π → ∞)
#
# and **no matter how large ω is, it cannot exceed 50%** (verified numerically). The non-move level required by
# the Time Use Survey (TUS) is lower than that, so it **cannot be reached by level alone.**
#
# The deeper reason is that `routine_axis_daytype` **preserves row sums** — it scales inter-district commutes
# by λ and sends the removed share to that row's diagonal, and the budget `f=λκ` (see Methods) splits that
# diagonal again, sending most of it to N̂M. That is, **a person who commutes on weekdays necessarily becomes
# non-move upon stopping commuting.** This is imposed by the model, not stated by the data.
#
# φ = **"the fraction of people who commute (work or school) on weekdays who, on that day type, stop commuting
#        and make a non-commute trip instead"**. The commuters released in each row are
#
#     free^g_j = P^g_j κ^g_j (1 − λ^{g,dt})
#
# (κ = census commuting (work or school) rate, λ = commute reduction factor); the φ share of these is subtracted
# from the routine-axis diagonal and **added to the same row of the non-commute axis** (destination shape = KTDB `V` as is).
#
#     Σ_i R^dt_ji = P^g_j − φ·free^g_j ,    Σ_i NC^dt_ji = π^{g,dt} P^g_j + φ·free^g_j
#
# The seed row total stays `P(1+π)`, so the IPF row constraint is unchanged. The realized NC share becomes
# `(π + φκ(1−λ))/(1+π)`, lifting the 50% ceiling. **With φ = 0 it is exactly identical to the baseline (no-φ) model.**

"""
    routine_axis_daytype(R, λ, free_phi, N) -> (Rdt, moved)

Carries the workday routine axis `R = C + NM` over to a day-type version.

1. Scale the off-diagonal (inter-district commutes) by `λ` and **add the removed share
   `(1-λ)·Σ_{i≠j}R_{ji}` to that row's diagonal** — up to here row sums are preserved.
2. **Subtract** `free_phi[j]` (= φ·free^g_j) from that diagonal. The diagonal is clamped so it does not go
   negative, and the **amount actually subtracted** is returned as `moved` — the caller must add exactly that to NC for totals to match.

With `free_phi .= 0` it is identical to the baseline `routine_axis_sb`; with `λ = 1` it equals the workday version.
The inner commute : non-move split within the diagonal is not done here — the active-budget split (see Methods) re-splits it with `f^{dt}=λκ`,
so the **same rule** applies to `mobile` (diagonal arrives only as a sum) and `radiation` (c_in·nm arrive separately).
"""
function routine_axis_daytype(R, λ::Float64, free_phi::Vector{Float64}, N)
    Rdt = λ .* R
    moved = zeros(N)
    for j in 1:N
        off = sum(@view R[j, :]) - R[j, j]
        dia = R[j, j] + (1 - λ) * off
        m = min(free_phi[j], dia)
        moved[j] = m
        Rdt[j, j] = dia - m
    end
    (Rdt, moved)
end

"""
    noncommute_guess_g(V, P, πg, moved, g, N)

**Single-age-group** version of `noncommute_guess` in `src/model.jl`. The destination shape (KTDB `V`) is
used as is; only the row total is set to `πg·P^g_j + moved[j]`. Rows of `V` that are entirely 0 are skipped,
as in `noncommute_guess`. Age groups are fully independent, so this is split out for the caller to run with `Threads.@threads`.
"""
function noncommute_guess_g(V, P, πg::Float64, moved::Vector{Float64}, g::String, N)
    A = zeros(N, N)
    for j in 1:N
        row = @view V[g][j, :]
        s = sum(row)
        s > 0 || continue
        A[j, :] .= ((πg * P[g][j] + moved[j]) / s) .* row
    end
    A
end

"""
    clean_zero_routine!(sp, Rdt, N)

For age groups with `λ^g = 0` the routine-axis off-diagonal is **exactly 0**, so `Ĉ_{ji} (i≠j)` must
also be 0. However, the axis split in `active_budget_split` computes `N̂C = M̂·NC/(R+NC)` → `R̂ = M̂ − N̂C`
in that order, so floating-point residuals (±10⁻¹²) remain in `R̂` at cells with `R=0`.

This dust is removed by **moving it to `N̂C`** — `Ĉ + N̂C` is preserved per cell, so the row identity
`Ĉ + N̂C + N̂M = u` still holds. Returns the number of cells moved and the total amount (for validation output).
"""
function clean_zero_routine!(sp, Rdt, N)
    n = 0; tot = 0.0
    for j in 1:N, i in 1:N
        i == j && continue
        if Rdt[j, i] == 0.0 && sp.Chat[j, i] != 0.0
            tot += abs(sp.Chat[j, i]); n += 1
            sp.NChat[j, i] += sp.Chat[j, i]
            sp.Chat[j, i] = 0.0
        end
    end
    (; n, tot)
end

"""
    zero_breakdown(A, N)

Classifies off-diagonal cells by their value.
The OD csv contains only cells with `A[j,i] > 0`, so "zero cells" appear in the file only as **missing rows**.

  * `n_neg`  — negative (floating-point residual)
  * `n_zero` — **exactly 0** (structural zero) — no row in the csv
  * `n_tiny` — `0 < A ≤ 5e-4` — row present but fewer than 0.0005 persons per day
  * `n_pos`  — `A > 5e-4`
"""
function zero_breakdown(A, N)
    n_neg = n_zero = n_tiny = n_pos = 0
    vmin = 0.0; tiny_max = 0.0
    for j in 1:N, i in 1:N
        i == j && continue
        v = A[j, i]
        if v < 0
            n_neg += 1; v < vmin && (vmin = v)
        elseif v == 0
            n_zero += 1
        elseif v <= 5e-4
            n_tiny += 1; v > tiny_max && (tiny_max = v)
        else
            n_pos += 1
        end
    end
    (; n_neg, n_zero, n_tiny, n_pos, vmin, tiny_max)
end

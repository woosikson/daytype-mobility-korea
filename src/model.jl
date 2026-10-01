# src/model.jl — assemble initial guess → IPF → active-budget split.
#   Matrix convention as in src/inputs.jl: A[j,i], row = origin j, column = destination i.

"""
    noncommute_guess(V, P, π)

Non-commute initial guess. The destination distribution takes the shape of the KTDB trip OD `V` as is,
and the row total is set solely by the survey participation rate `π^g`.

    NC^g_{ij} = P^g_j · π^g · V^g_{ij} / Σ_{i'} V^g_{i'j}   →   Σ_i NC^g_{ij} = π^g P^g_j
"""
function noncommute_guess(V, P, π, N)
    NC = Dict(g => zeros(N, N) for g in AG)
    for g in AG, j in 1:N
        row = @view V[g][j, :]
        s = sum(row)
        s > 0 || continue
        NC[g][j, :] .= (π[g] * P[g][j] / s) .* row
    end
    NC
end

"""
    routine_axis(seed, g, j-scope...)

Routine axis `C^g_{ij} + NM^g_{ij}` (commute + non-move) of the initial guess.
The two seeds differ only in how the diagonal is set:

  * `mobile`     — observed KT OD as is (diagonal = inner commute + non-move, not separable)
  * `radiation` — off-diagonal = c_out × radiation shape, diagonal = c_in + nm (census decomposition)
"""
function routine_axis(seed::String, g::String, CA, comp, G, NC, N)
    if seed == "mobile"
        return copy(CA[g])
    elseif seed == "radiation"
        R = comp.c_out[g] .* G
        for j in 1:N; R[j, j] = comp.c_in[g][j] + comp.nm[g][j]; end
        return R
    end
    error("unknown seed: $seed")
end

"""
    ipf_fit(guess, u, v)

Biproportional fitting of the initial guess to the two marginals (`ProportionalFitting.jl`).
Margins are passed as proportions; after convergence the total is scaled to Σ_j P^g_j (= Σ u).
"""
function ipf_fit(guess, u, v; maxiter = 5000, tol = 1e-9)
    fac = ipf(guess, [u ./ sum(u), v ./ sum(v)]; maxiter = maxiter, tol = tol)
    M = Array(fac) .* guess
    M .*= sum(u) / sum(M)
    M
end

"""
    active_budget_split(M, R, NC, P, f)

Splits the IPF result into three components.

  1. axis split (preserves initial-guess proportions):  R̂ = M̂·R/(R+NC),  N̂C = M̂·NC/(R+NC)
  2. diagonal split (active budget):
        N̂M_j = clamp(Σ_i R̂_ij − P_j f_j, 0, R̂_jj),   Ĉ_jj = R̂_jj − N̂M_j,  Ĉ_ij = R̂_ij (i≠j)

Returns: `(Chat, NChat, nm, nm_raw, n_lo, n_hi)` — `nm_raw` is the pre-clamp value (for diagnostics).
"""
function active_budget_split(M, R, NC, P, f, N)
    den   = R .+ NC
    NChat = ifelse.(den .> 0, M .* NC ./ den, 0.0)
    Rhat  = M .- NChat
    Chat  = copy(Rhat)
    nm     = zeros(N); nm_raw = zeros(N)
    n_lo = 0; n_hi = 0
    for j in 1:N
        rowR = sum(@view Rhat[j, :])
        raw  = rowR - P[j] * f[j]
        cap  = Rhat[j, j]
        raw < 0   && (n_lo += 1)
        raw > cap && (n_hi += 1)
        b = clamp(raw, 0.0, cap)
        nm[j] = b; nm_raw[j] = raw
        Chat[j, j] = cap - b
    end
    (; Chat, NChat, Rhat, nm, nm_raw, n_lo, n_hi)
end

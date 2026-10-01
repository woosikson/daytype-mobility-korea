# src/pipeline.jl — the pipeline of Eqs (2a)–(2g), written to take **one month's input bundle (`Ctx`)**.
#
#   Since we loop over 72 months, one month's inputs are packed into a single `Ctx` and passed as an argument.
#   Matrix convention: A[j,i], row = origin j, column = destination i (src/inputs.jl).
#
#   Varies by month : P · de facto marginal ṽ · census 3-way split (fixed shares × P) · κ (materially only for 00·10) ·
#                      λ^{a,e}·ω^{a,e} (Sat/Sun/weekday-holiday composition of that month's weekend) · 75+ band split
#   Fixed across months: radiation shape G · KTDB V · π^{a,w} · TUS source tables · posterior of θ
#   The `mobile` seed (validation only) is built for 2022-09 only (`MOB` argument).

struct Ctx
    y::Int
    m::Int
    dts::Vector{String}                                 # day types present in that month
    comp_days::Dict{String,NamedTuple}                  # day-type composition (n, n_sat, n_sun, n_wkh, kind)
    P::Dict{String,Vector{Float64}}                     # P^a_j = u^a_j (origin, shared by all three day types)
    SUMP::Dict{String,Float64}
    PTOT::Float64
    VT::Dict{String,Dict{String,Vector{Float64}}}       # dt => a => ṽ^{a,t}_i
    comp::NamedTuple                                    # (c_in, nm, c_out)
    κ::Dict{String,Vector{Float64}}
    RWD::Dict{String,Dict{String,Matrix{Float64}}}      # seed => a => weekday routine axis
    LAM::Dict{String,Dict{String,Float64}}              # weekday 1 · weekend λ^{a,e} · holiday = baseline to be multiplied by χ
    OMB::Dict{String,Dict{String,Float64}}              # weekday 1 · weekend ω^{a,e} · holiday = baseline to be multiplied by ψ
    W7579::Float64                                      # 75+ → 70~79 / 80+ split (that month's registered population)
    PBAND::Dict{String,Float64}
    P10::Float64
    q_tus::Dict{String,NamedTuple}                      # band => (weekday, weekend) — weighted by that month's composition
    T1_OBS::Float64
    T2_OBS::Vector{Float64}
    nadd::Int                                           # number of cells with additive correction
end

ymlabel(c::Ctx) = ymstr(c.y, c.m)

# ── TUS day-of-week weighting — weekday public holidays are treated like Sundays (ref run.jl §1) ──────
tusw(w::NamedTuple, v) = (w.n_sat * v[2] + (w.n_sun + w.n_wkh) * v[3]) / w.n
tusw_wd(v) = v[1]
"Holiday baseline — **ordinary weekend** = Sat : Sun 1 : 1 (the denominator of ψ is the MOLIT \"ordinary weekend\")."
const STD_WEEKEND = (n = 2, n_sat = 1, n_sun = 1, n_wkh = 0, kind = "")

"""
    lam_omega(w, W7579, tus)

Build λ^{a,e}·ω^{a,e} of Eq. (3) from the weekend composition `w` (same as run.jl §3).
"""
function lam_omega(w::NamedTuple, W7579::Float64, tus)
    bw(g) = g == "75" ? [W7579, 1 - W7579] : [1.0]
    rate(g, code, wd::Bool) = sum(x * (wd ? tusw_wd(tus[(b, code)]) : tusw(w, tus[(b, code)]))
                                  for (b, x) in zip(BAND[g], bw(g)))
    LAM = Dict{String,Float64}(); OMB = Dict{String,Float64}()
    rows = NamedTuple[]
    for g in AG
        c_wd, s_wd = rate(g, "921", true), rate(g, "930", true)
        c_e, s_e = rate(g, "921", false), rate(g, "930", false)
        LAM[g] = g in ("00", "10") ? 0.0 : c_e / (c_wd + s_wd)
        OMB[g] = uni(vcat([rate(g, x, false) for x in NC_NARROW], s_e)) /
                 uni([rate(g, x, true) for x in NC_NARROW])
        push!(rows, (age_group = g, band = join(BAND[g], "+"), commute_wd = c_wd, school_wd = s_wd,
                     commute_weekend = c_e, school_weekend = s_e,
                     lambda_weekend = LAM[g], omega_base = OMB[g]))
    end
    (LAM, OMB, rows)
end

"Values by age group → sums over the 8 TUS bands (`00` is excluded)."
function to_band(x::Dict{String,Float64}, W7579::Float64)
    d = Dict(b => 0.0 for b in TUS_BANDS)
    for (g, bs) in G2BAND
        ws = g == "75" ? [W7579, 1 - W7579] : [1.0]
        for (bb, w) in zip(bs, ws); d[bb] += w * x[g]; end
    end
    d
end
to_band(c::Ctx, x::Dict{String,Float64}) = to_band(x, c.W7579)

# ═════════════════════════════════════════════════════════════════════════════
# Functions that solve at a single θ
# ═════════════════════════════════════════════════════════════════════════════
"The two seed axes of age group g for a given day type and θ (Methods, steps ①–⑤)."
function seed_axes(c::Ctx, dt::String, seed::String, g::String, φ::Float64, ψ::Float64, χ::Float64)
    λg = dt == "holiday" ? min(χ * c.LAM["holiday"][g], 1.0) : min(c.LAM[dt][g], 1.0)
    ωg = dt == "holiday" ? ψ * c.OMB["holiday"][g] : c.OMB[dt][g]
    πg = min(ωg * PI_WD[g], 0.999)
    R  = c.RWD[seed][g]
    free = dt == "weekday" ? zeros(N) : φ .* c.P[g] .* c.κ[g] .* (1 - λg)
    Rd, moved = dt == "weekday" ? (R, zeros(N)) : routine_axis_daytype(R, λg, free, N)
    NC = noncommute_guess_g(V, c.P, πg, moved, g, N)
    (; Rd, NC, λg, πg, moved)
end

"Solve a single age group g at one θ end to end — (2f) IPF · (2g) split."
function solve_group(c::Ctx, dt::String, seed::String, g::String, θ)
    ax = seed_axes(c, dt, seed, g, dt == "weekday" ? 0.0 : Float64(θ[1]), Float64(θ[2]), Float64(θ[3]))
    M  = ipf_fit(ax.Rd .+ ax.NC, c.P[g], c.VT[dt][g])
    sp = active_budget_split(M, ax.Rd, ax.NC, c.P[g], ax.λg .* c.κ[g], N)
    clean_zero_routine!(sp, ax.Rd, N)
    (; sp, M, λg = ax.λg, πg = ax.πg)
end

"National **totals** of the 3 components by age group at one θ only (for the grid and draws)."
function solve(c::Ctx, dt::String, seed::String; φ = 0.0, ψ = 1.0, χ = 0.5)
    res = Vector{NTuple{3,Float64}}(undef, length(AG))
    Threads.@threads for k in eachindex(AG)
        sp = solve_group(c, dt, seed, AG[k], (φ, ψ, χ)).sp
        res[k] = (sum(sp.Chat), sum(sp.NChat), sum(sp.nm))
    end
    Dict(AG[k] => res[k] for k in eachindex(AG))
end

"3-component Dict → (10+ commute participation rate, travel participation rate q by band, national C·NC·NM)."
function summarise(c::Ctx, r::Dict{String,NTuple{3,Float64}})
    nmb = to_band(c, Dict(g => r[g][3] for g in AG if g != "00"))
    (; c10 = sum(r[g][1] for g in AG if g != "00") / c.P10,
       q  = Dict(b => 1 - nmb[b] / c.PBAND[b] for b in TUS_BANDS),
       c  = sum(r[g][1] for g in AG), nc = sum(r[g][2] for g in AG), nm = sum(r[g][3] for g in AG))
end

"Same as `solve` but returns full matrices and adm2-level values (for outputs)."
function solve_full(c::Ctx, dt::String, seed::String, θ)
    out = Vector{Any}(undef, length(AG))
    Threads.@threads for k in eachindex(AG)
        g = AG[k]
        r = solve_group(c, dt, seed, g, (θ.φ, θ.ψ, θ.χ))
        sp = r.sp
        dcol = vec(sum(r.M, dims = 1))
        out[k] = (; g, sp, λg = r.λg, πg = r.πg, dcol,
                  dm = maximum(abs.(dcol .- c.VT[dt][g])),
                  id = maximum(abs.(vec(sum(sp.Chat, dims = 2)) .+ vec(sum(sp.NChat, dims = 2)) .+
                                    sp.nm .- c.P[g])))
    end
    out
end

# ═════════════════════════════════════════════════════════════════════════════
# Build one month's input bundle
# ═════════════════════════════════════════════════════════════════════════════
"""
    make_ctx(y, m; REG, RAW, CMT, CAL, MOB, tus, tt)

Collect that month's marginals · census 3-way split · weekday routine axes of both seeds · λ·ω · targets.
`MOB` is that month's mobile-phone corrected OD `Dict(g => A)` (if absent, no `mobile` seed is built).
"""
function make_ctx(y, m; REG, RAW, CMT, CAL, MOB, tus, tt)
    cd = month_composition(CAL, y, m)
    dts = [dt for dt in ("weekday", "weekend", "holiday") if haskey(cd, dt)]
    P = REG.P[(y, m)]
    SUMP = Dict(g => sum(P[g]) for g in AG); PTOT = sum(values(SUMP))
    VT, nadd = daytype_marginals(RAW[(y, m)], P, S, dts)
    comp = month_components(CMT, y, m, P, REG.p03[(y, m)], REG.p1214[(y, m)], S)
    κ = census_commute_rate(comp, P, N)
    RWD = Dict{String,Dict{String,Matrix{Float64}}}()
    RWD["radiation"] = Dict(g => max.(routine_axis("radiation", g, nothing, comp, Gshape, nothing, N), 0.0)
                            for g in AG)
    MOB === nothing || (RWD["mobile"] = Dict(g => max.(MOB[g], 0.0) for g in AG))
    W7579 = REG.p7579[(y, m)] / (REG.p7579[(y, m)] + REG.p80[(y, m)])
    lw, ow, _ = lam_omega(cd["weekend"], W7579, tus)
    lh, oh, _ = lam_omega(STD_WEEKEND, W7579, tus)
    LAM = Dict("weekday" => Dict(g => 1.0 for g in AG), "weekend" => lw, "holiday" => lh)
    OMB = Dict("weekday" => Dict(g => 1.0 for g in AG), "weekend" => ow, "holiday" => oh)
    PBAND = to_band(Dict(g => SUMP[g] for g in AG if g != "00"), W7579)
    P10 = sum(values(PBAND))
    q_tus = Dict(b => (weekday = tt[b][1] / 100, weekend = tusw(cd["weekend"], tt[b]) / 100)
                 for b in TUS_BANDS)
    T1 = tusw(cd["weekend"], (TUS_COMMUTE_10P.weekday, TUS_COMMUTE_10P.saturday, TUS_COMMUTE_10P.sunday))
    T2 = [q_tus[b].weekend - q_tus[b].weekday for b in TUS_BANDS]
    Ctx(y, m, dts, cd, P, SUMP, PTOT, VT, comp, κ, RWD, LAM, OMB, W7579, PBAND, P10, q_tus, T1, T2, nadd)
end

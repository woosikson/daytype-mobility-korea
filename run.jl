# run.jl — day-type-resolved district OD matrices, 2018-01 … 2023-12
#
# Monthly **weekday · weekend · holiday** population movement between districts, 2018-01 … 2023-12 (adm2 250 × age group 15).
#
#   * Eqs. (1), (2a)–(2g), (3) and the calibration: see the Methods of the Data Descriptor. Main seed `radiation`; `mobile` only for 2022-09 validation.
#   * **θ = (φ, ψ, χ): the posterior median θ̂ of the 2022-09 calibration is applied unchanged to all 72 months.**
#   * All inputs are read from `input/` and `contents/`.
#
# Notation: the paper writes A_{ij} (i = destination, j = origin); the code stores `X[j, i]`.
# No plotting here — `plot.R` reads only `data/*.csv`.
#
# Layout
#   Part A  inputs    §1 calendar · §2 monthly inputs (2022-09 reproduction check) · §3 TUS
#   Part B  2022-09   §4 weekday · §5 grid · §6 posterior · §7 outputs · §8 validation · §9 figure summaries
#   Part C  monthly   §10 72-month OD · §11 weekday three-component time series
using CSV, DataFrames, ProportionalFitting, Printf, Statistics, Logging, Distributions, SHA, Dates

global_logger(ConsoleLogger(stderr, Logging.Warn))     # silence IPF convergence @info

const HERE  = @__DIR__
const INP   = joinpath(HERE, "input")
const CONT  = joinpath(HERE, "contents")
const DATA  = joinpath(HERE, "data");  isdir(DATA) || mkpath(DATA)
const ODDIR = joinpath(DATA, "od");    isdir(ODDIR) || mkpath(ODDIR)

com(x) = replace(@sprintf("%d", round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")
"Storage precision — 3 decimals if ≥ 1 person, else 6 significant digits."
mv(x) = abs(x) >= 1 ? round(x, digits = 3) : round(x, sigdigits = 6)
sha16(p) = bytes2hex(open(sha256, p))[1:16]
hr(t) = (println(); println("-"^78); println(t); flush(stdout))

include(joinpath(HERE, "src", "inputs.jl"))    # AG, AG3, load_*, census_commute_rate …
include(joinpath(HERE, "src", "model.jl"))     # routine_axis, ipf_fit, active_budget_split
include(joinpath(HERE, "src", "daytype.jl"))   # routine_axis_daytype(φ), noncommute_guess_g, clean_zero_routine!
include(joinpath(HERE, "src", "tus.jl"))       # TUS PDF parsing
include(joinpath(HERE, "src", "molit.jl"))     # MOLIT ψ bracket
include(joinpath(HERE, "src", "visitor.jl"))   # Korea Tourism Data Lab visitor data comparison
include(joinpath(HERE, "src", "monthly.jl"))   # monthly inputs (calendar · registered population · de facto · census)
include(joinpath(HERE, "src", "pipeline.jl"))  # Ctx · seed_axes · solve · solve_full · make_ctx

# ── Fixed settings ──────────────────────────────────────────────────────────────
const DTS     = ("weekday", "weekend", "holiday")
const SEEDS   = ("radiation", "mobile")               # first element is the **main** seed
const CALSEED = "radiation"
const Y0, M0  = 2022, 9                               # calibration month
# TUS (2024) Vol. 1-1, Table 1-2: participation rate of activity 921 (commuting to work), age 10+ — input to T1.
const TUS_COMMUTE_10P = (weekday = 0.507, saturday = 0.225, sunday = 0.147)

println("="^78)
println("Day-type mobility matrices — 2018-01 … 2023-12, θ = 2022-09 posterior median")
println("  threads = $(Threads.nthreads())")
println("="^78)
const T_START = time()

# ═════════════════════════════════════════════════════════════════════════════
# Part A — inputs
# ═════════════════════════════════════════════════════════════════════════════
# ── §1 calendar ───────────────────────────────────────────────────────────────
hr("1. Calendar — 2018-01 … 2023-12")
const CAL = load_calendar(joinpath(INP, "daytype_2018-2023.csv"))
@assert nrow(CAL) == 2191
let rows = NamedTuple[]
    for (y, m) in YM_ALL, (dt, w) in month_composition(CAL, y, m)
        push!(rows, (year = y, month = m, ym = ymstr(y, m), daytype = dt,
                     n_days = w.n, n_sat = w.n_sat, n_sun = w.n_sun, n_wkh = w.n_wkh, kind = w.kind))
    end
    df = sort(DataFrame(rows), [:year, :month, :daytype])
    CSV.write(joinpath(DATA, "month_daytype_composition.csv"), df)
    h = df[df.daytype .== "holiday", :]
    @printf("  months with a holiday: %d · months whose weekend includes a stand-alone weekday public holiday: %d\n",
            nrow(h), count(df.daytype .== "weekend" .&& df.n_wkh .> 0))
    for k in ("seollal", "chuseok", "other")
        @printf("    %-8s %s\n", k, join((r.ym * "(" * string(r.n_days) * ")" for r in eachrow(h) if occursin(k, r.kind)), " "))
    end
end
let w = month_composition(CAL, Y0, M0)
    @assert (w["weekday"].n, w["weekend"].n, w["holiday"].n) == (20, 6, 4)
    @assert (w["weekend"].n_sat, w["weekend"].n_sun, w["weekend"].n_wkh) == (3, 3, 0)
    println("  2022-09: weekday 20 · weekend 6 (Sat 3 · Sun 3) · holiday 4 (Chuseok)")
end

# ── §2 monthly inputs ──────────────────────────────────────────────────────────
hr("2. Monthly inputs — 2022-09 checked against the 2022-09 input files")

# adm2 250 — ascending code order
const S = sort(unique(CSV.read(joinpath(INP, "marginals_202209.csv"), DataFrame;
                               select = [:sgg_cd], types = Dict(:sgg_cd => String)).sgg_cd))
const N = length(S); const IDX = Dict(s => k for (k, s) in enumerate(S))
@assert N == 250

t0 = time()
const REG = load_registered(joinpath(INP, "kosis_registered_pop_sgg_sexage_2018-2023.csv"), S, IDX)
const SGGNM = REG.sggnm
@printf("  registered population, 72 months (%.0f s) — 2018-01 %s · 2022-09 %s · 2023-12 %s persons\n", time() - t0,
        com(sum(sum(REG.P[(2018, 1)][g]) for g in AG)), com(sum(sum(REG.P[(2022, 9)][g]) for g in AG)),
        com(sum(sum(REG.P[(2023, 12)][g]) for g in AG)))

VERIN = NamedTuple[]            # input reproduction check (2022-09: 72-month rules vs 2022-09 input files)
let mg = CSV.read(joinpath(INP, "marginals_202209.csv"), DataFrame; types = Dict(:sgg_cd => String, :age_group => String))
    d = maximum(abs(r.registered - REG.P[(2022, 9)][r.age_group][IDX[r.sgg_cd]]) for r in eachrow(mg))
    push!(VERIN, (input = "registered P (marginals_202209)", max_abs_diff = d, max_rel_diff = NaN))
    @printf("  [repro] registered 2022-09 vs marginals_202209.csv : max diff %.3g persons\n", d)
end

const DFF = joinpath(INP, "defacto_daytype_2018-2023.csv")
const RAW = load_defacto_daytype(DFF)
@assert length(RAW) == 72 && all(length(v) == N * length(AG) for v in values(RAW))
# Day counts — do they match the calendar? Means are taken over **days on which that tz exists**. tz15 must be complete;
# tz03 is used only for weekday (Eq. 1b), so weekday tz03 must be complete too.
let dd = unique(CSV.read(DFF, DataFrame; select = [:year, :month, :daytype, :timezn, :n_days],
                         types = Dict(:timezn => String)))
    short = NamedTuple[]
    for r in eachrow(dd)
        n = month_composition(CAL, r.year, r.month)[r.daytype].n
        r.n_days == n && continue
        (r.timezn == "03" && r.daytype != "weekday" && r.n_days < n) ||
            error("de facto day count differs from calendar: $(r.year)-$(r.month) $(r.daytype) tz$(r.timezn) $(r.n_days) ≠ $n")
        push!(short, (ym = ymstr(r.year, r.month), daytype = r.daytype, timezn = r.timezn, n_days = r.n_days, calendar = n))
    end
    CSV.write(joinpath(DATA, "defacto_missing_days.csv"), DataFrame(short))
    println("  de facto: 72 months × adm2 250 × age 15 · tz15 and weekday tz03 day counts all match the calendar")
    for s in short
        @printf("    ⚠️ %s %s tz%s has %d / %d days (days missing in source) — not used by Eq. (1b), no effect on results\n",
                s.ym, s.daytype, s.timezn, s.n_days, s.calendar)
    end
end
let rf = CSV.read(joinpath(INP, "defacto_daytype_202209.csv"), DataFrame;
                  types = Dict(:sgg_cd => String, :age_group => String, :timezn => String, :daytype => String))
    rf = rf[in.(rf.timezn, Ref(["03", "15"])), :]
    d = maximum(abs(r.raw_mean - RAW[(2022, 9)][(r.sgg_cd, r.age_group)][(DTMAP[r.daytype], r.timezn)]) for r in eachrow(rf))
    push!(VERIN, (input = "de facto raw_mean (defacto_daytype_202209)", max_abs_diff = d, max_rel_diff = NaN))
    @printf("  [repro] de facto 2022-09 vs defacto_daytype_202209.csv : max diff %.3g persons (%d rows)\n", d, nrow(rf))
end

const CMT = load_census_comm(joinpath(INP, "sgg_commute_agegroup_2018-2023.csv"))
let rc = CSV.read(joinpath(INP, "regular_components_202209.csv"), DataFrame; types = Dict(:sgg_cd => String, :age_group => String))
    cp = month_components(CMT, 2022, 9, REG.P[(2022, 9)], REG.p03[(2022, 9)], REG.p1214[(2022, 9)], S)
    d = maximum(max(abs(r.gwannae - cp.c_in[r.age_group][IDX[r.sgg_cd]]), abs(r.bimove - cp.nm[r.age_group][IDX[r.sgg_cd]]),
                    abs(r.tasgg - cp.c_out[r.age_group][IDX[r.sgg_cd]])) for r in eachrow(rc))
    push!(VERIN, (input = "census components (regular_components_202209)", max_abs_diff = d, max_rel_diff = NaN))
    @printf("  [repro] census 3-way split 2022-09 vs regular_components_202209.csv : max diff %.3g persons\n", d)
end

# KTDB non-commute destination shape — **age-independent** V_ij (X[j,i], row = origin).
#   The input's age-group axis is V_ij split by age weights constant per origin, so it carries no destination-shape information
#   (it cancels in the row normalisation of Eq. (9)). Summing the 15 age groups recovers V_ij, used for every age group.
#   The age-free file (`nc_od_ageless_capped.csv`) is not used: 3-decimal rounding turns 163 zero-fill cells into 0.
const VALL = let Vg = load_od(joinpath(INP, "nc_od_agegroup.csv"), IDX, N;
                             from = :from_cd_stable, to = :to_cd_stable, age = :age_group, val = :nc_trips)
    sum(Vg[g] for g in AG)
end
@assert count(==(0.0), VALL[j, i] for j in 1:N, i in 1:N if i != j) == 0
const V = Dict(g => VALL for g in AG)            # same matrix for every age group (keeps the noncommute_guess_g interface)
const PI_WD  = load_pi(joinpath(INP, "nc_person_row_factor.csv")).π
const Gshape = load_radiation_shape(joinpath(INP, "radiation_T_ij_2020.csv"), IDX, N)

# ── §3 Time Use Survey ───────────────────────────────────────────────────────
hr("3. TUS → λ · ω · ψ bracket · χ prior")
const TUS, TUS_BANDS_READ = read_tus_travel(CONT)
const TT = read_tus_travel_total(CONT)
println("  Table 2-2 bands: ", length(TUS_BANDS_READ))
const PR_CHI  = Beta(2.0, 2.0)
const CHI_MID = median(PR_CHI)
@printf("  ψ bracket [%.3f, %.4f] · prior median %.4f · χ prior Beta(2,2) median %.2f\n",
        PSI_MIN, PSI_MAX, PSI_MID, CHI_MID)

# ═════════════════════════════════════════════════════════════════════════════
# Part B — 2022-09 : calibration and 2022-09 outputs
# ═════════════════════════════════════════════════════════════════════════════
hr("4. 2022-09 — input bundle · weekday")
# Mobile-phone night–day residence OD 2022-09 — `mobile` seed (validation only). X[j, i] : row = night-time district j, column = daytime district i
const MOB9 = let rf = CSV.read(joinpath(INP, "od_corrected_202209.csv"), DataFrame;
                               types = Dict(:ngt_sgg_cd => String, :day_sgg_cd => String, :age_group => String))
    X = Dict(g => zeros(N, N) for g in AG)
    for r in eachrow(rf); X[r.age_group][IDX[r.ngt_sgg_cd], IDX[r.day_sgg_cd]] = r.pop; end
    X
end
mkctx(y, m, mob) = make_ctx(y, m; REG, RAW, CMT, CAL, MOB = mob, tus = TUS, tt = TT)
const C9 = mkctx(Y0, M0, MOB9)
let mg = load_marginals(joinpath(INP, "marginals_202209.csv"))
    d = maximum(let v = mg.d15[g] .* (C9.SUMP[g] / sum(mg.d15[g]))
                    maximum(abs.(C9.VT["weekday"][g] .- v) ./ max.(v, 1.0)) end for g in AG)
    push!(VERIN, (input = "weekday destination marginal (marginals_202209 defacto15_corr)", max_abs_diff = NaN, max_rel_diff = d))
    @printf("  [repro] weekday ṽ vs marginals_202209.csv defacto15_corr : max rel diff %.3g (storage rounding level)\n", d)
end
CSV.write(joinpath(DATA, "verification_inputs.csv"), DataFrame(VERIN))

const P, SUMP, PTOT = C9.P, C9.SUMP, C9.PTOT
@printf("  adm2 %d · age group %d · registered %s persons · correction additions %d cells\n", N, length(AG), com(PTOT), C9.nadd)
let (_, _, rows) = lam_omega(C9.comp_days["weekend"], C9.W7579, TUS)
    CSV.write(joinpath(DATA, "lambda_omega_agegroup.csv"), DataFrame([
        merge(r, (registered = mv(SUMP[r.age_group]),
                  kappa_mean = mv(sum(C9.κ[r.age_group] .* P[r.age_group]) / SUMP[r.age_group]),
                  pi_weekday = mv(PI_WD[r.age_group]))) for r in rows]))
    wP = [SUMP[g] / PTOT for g in AG]
    @printf("  75+ split 70~79 : 80+ = %.4f : %.4f · registered-weighted λ^e %.3f · ω^e %.3f\n", C9.W7579, 1 - C9.W7579,
            sum(wP .* [C9.LAM["weekend"][g] for g in AG]), sum(wP .* [C9.OMB["weekend"][g] for g in AG]))
end
CSV.write(joinpath(DATA, "tus_travel_total.csv"), DataFrame([
    (band = b, weekday = TT[b][1], saturday = TT[b][2], sunday = TT[b][3],
     q_weekday = mv(C9.q_tus[b].weekday), q_weekend = mv(C9.q_tus[b].weekend),
     dq_weekend = mv(C9.q_tus[b].weekend - C9.q_tus[b].weekday)) for b in TUS_BANDS]))

const SWD = summarise(C9, solve(C9, "weekday", CALSEED))
@printf("  weekday (radiation): C %.2f%% · NC %.2f%% · NM %.2f%% · 10+ commute %.2f%%\n",
        100SWD.c/PTOT, 100SWD.nc/PTOT, 100SWD.nm/PTOT, 100SWD.c10)

# ── §5 grid ───────────────────────────────────────────────────────────────
hr("5. Grid (2022-09)")
const GPHI  = collect(range(0.0, 1.0, length = 41))
const GPSI  = exp.(range(log(0.80), log(1.60), length = 17))
const GPHIH = collect(range(0.0, 1.0, length = 41))        # holiday φ grid: same 0.025 spacing as weekend
const GCHI  = collect(range(0.0, 1.0, length = 13))
qcols(s) = NamedTuple{Tuple(Symbol.("q_" .* TUS_BANDS))}(Tuple(s.q[b] for b in TUS_BANDS))
row_of(dt, φ, ψ, χ, s) = merge((daytype = dt, phi = φ, psi = ψ, chi = χ, c10 = s.c10,
                                commute = s.c, noncommute = s.nc, nonmove = s.nm), qcols(s))
function build_grid()
    grid = NamedTuple[row_of("weekday", 0.0, 1.0, 1.0, SWD)]
    t0 = time(); done = 0
    tot = length(GPHI) + length(GPHIH) * length(GPSI) * length(GCHI)
    for φ in GPHI
        push!(grid, row_of("weekend", φ, 1.0, 1.0, summarise(C9, solve(C9, "weekend", CALSEED; φ = φ))))
        done += 1
    end
    for φ in GPHIH
        for ψ in GPSI, χ in GCHI
            push!(grid, row_of("holiday", φ, ψ, χ,
                               summarise(C9, solve(C9, "holiday", CALSEED; φ = φ, ψ = ψ, χ = χ))))
            done += 1
        end
        el = time() - t0
        @printf("    holiday φ = %.3f  %6d / %d  (%.0f s · %.1f min left)\n", φ, done, tot, el, (tot - done) * el / done / 60)
        flush(stdout)
    end
    @printf("  grid %s points · %.1f min\n", com(tot), (time() - t0) / 60)
    grid
end
const GRIDF = joinpath(DATA, "grid_summaries.csv")
const GRIDK = joinpath(DATA, "grid_key.txt")
const PARAMZN = "phi_scalar+psi+chi"
gridkey() = join(vcat(["PARAM", PARAMZN, "SEED", CALSEED],
                      ["GPHI"], string.(round.(GPHI, digits = 9)), ["GPSI"], string.(round.(GPSI, digits = 9)),
                      ["GPHIH"], string.(round.(GPHIH, digits = 9)), ["GCHI"], string.(round.(GCHI, digits = 9)),
                      ["LAM"], [string(round(C9.LAM["weekend"][g], digits = 9)) for g in AG],
                      ["OMB"], [string(round(C9.OMB["weekend"][g], digits = 9)) for g in AG],
                      ["VT"], [string(round(sum(C9.VT[dt][g]), digits = 3)) for dt in DTS for g in AG],
                      ["V"], [bytes2hex(sha256(Vector{UInt8}(reinterpret(UInt8, vec(VALL)))))[1:16]]), "|")
if isfile(GRIDF) && isfile(GRIDK) && read(GRIDK, String) == gridkey()
    println("  reusing data/grid_summaries.csv (grid axes · λ · ω · seed unchanged)")
else
    CSV.write(GRIDF, DataFrame(build_grid()))
    write(GRIDK, gridkey())
end

# ── §6 posterior ───────────────────────────────────────────────────────────
hr("6. posterior (2022-09)")
gs = CSV.read(GRIDF, DataFrame)
const QSYM   = Symbol.("q_" .* TUS_BANDS)
const SUMCOL = vcat([:c10, :commute, :noncommute, :nonmove], QSYM)
wdref = Dict(c => gs[gs.daytype .== "weekday", c][1] for c in SUMCOL)
const LGPSI = log.(GPSI)
TWE = Dict(c => zeros(length(GPHI)) for c in SUMCOL)
THO = Dict(c => zeros(length(GPHIH), length(GPSI), length(GCHI)) for c in SUMCOL)
let iφ = Dict(v => k for (k, v) in enumerate(GPHI)), iψ = Dict(v => k for (k, v) in enumerate(GPSI)),
    ih = Dict(v => k for (k, v) in enumerate(GPHIH)), ic = Dict(v => k for (k, v) in enumerate(GCHI))
    for r in eachrow(gs[gs.daytype .== "weekend", :]), c in SUMCOL; TWE[c][iφ[r.phi]] = r[c]; end
    for r in eachrow(gs[gs.daytype .== "holiday", :]), c in SUMCOL
        THO[c][ih[r.phi], iψ[r.psi], ic[r.chi]] = r[c]
    end
end
function locate(ax, x)
    x = clamp(x, first(ax), last(ax))
    k = clamp(searchsortedlast(ax, x), 1, length(ax) - 1)
    (k, (x - ax[k]) / (ax[k+1] - ax[k]))
end
wev(c, φ) = ((i, t) = locate(GPHI, φ); (1 - t) * TWE[c][i] + t * TWE[c][i+1])
function hov(c, φ, ψ, χ)
    (i, ti) = locate(GPHIH, φ); (j, tj) = locate(LGPSI, log(ψ)); (k, tk) = locate(GCHI, χ)
    A = THO[c]
    ((1-ti)*(1-tj)*(1-tk)*A[i,j,k] + ti*(1-tj)*(1-tk)*A[i+1,j,k] + (1-ti)*tj*(1-tk)*A[i,j+1,k] +
     ti*tj*(1-tk)*A[i+1,j+1,k] + (1-ti)*(1-tj)*tk*A[i,j,k+1] + ti*(1-tj)*tk*A[i+1,j,k+1] +
     (1-ti)*tj*tk*A[i,j+1,k+1] + ti*tj*tk*A[i+1,j+1,k+1])
end
# The only target is T2 (weekend Δq in 8 bands). T1 (weekend 10+ commute) barely moves θ
# (range 1.35e-5 over the grid) and comes from the same TUS as λ, so it is kept out of the likelihood and used for **validation** only.
const T1_OBS, T1_SD = C9.T1_OBS, 0.015
const T2_OBS, T2_SD = C9.T2_OBS, 0.05
lnorm(z) = -0.5 * log(2π) - 0.5 * z^2
@printf("  target T2 (σ %.3f): %s  · validation T1 = %.4f\n", T2_SD, join((@sprintf("%+.3f", x) for x in T2_OBS), " "), T1_OBS)
const PR_PSI = (log(PSI_MID), 0.10)
const FPHI = collect(range(0.0, 1.0, length = 81))
const FPSI = exp.(range(log(GPSI[1]), log(GPSI[end]), length = 81))
const FCHI = collect(range(0.0, 1.0, length = 41))

"Joint posterior — T2 constrains weekend only, so the likelihood depends on φ alone (ψ, χ follow their priors)."
function posterior(phi_free::Bool; chi_prior = PR_CHI)
    φax = phi_free ? FPHI : [0.0]
    nf, np, nc = length(φax), length(FPSI), length(FCHI)
    W = zeros(nf, np, nc)
    QQ = Dict(k => zeros(nf, np, nc) for k in
              (:phi, :psi, :chi, :we_c, :we_nc, :we_nm, :we_c10, :ho_c, :ho_nc, :ho_nm))
    for (i, φ) in enumerate(φax)
        ll_we = 0.0
        for (k, bd) in enumerate(QSYM)
            ll_we += lnorm((wev(bd, φ) - wdref[bd] - T2_OBS[k]) / T2_SD) - log(T2_SD)
        end
        we = (wev(:commute, φ), wev(:noncommute, φ), wev(:nonmove, φ), wev(:c10, φ))
        for (j, ψ) in enumerate(FPSI)
            lpψ = lnorm((log(ψ) - PR_PSI[1]) / PR_PSI[2]) - log(PR_PSI[2])
            for (l, χ) in enumerate(FCHI)
                lpχ = log(max(pdf(chi_prior, χ), 1e-300))
                W[i,j,l] = ll_we + lpψ + lpχ
                QQ[:phi][i,j,l] = φ; QQ[:psi][i,j,l] = ψ; QQ[:chi][i,j,l] = χ
                QQ[:we_c][i,j,l] = we[1]; QQ[:we_nc][i,j,l] = we[2]
                QQ[:we_nm][i,j,l] = we[3]; QQ[:we_c10][i,j,l] = we[4]
                QQ[:ho_c][i,j,l] = hov(:commute, φ, ψ, χ); QQ[:ho_nc][i,j,l] = hov(:noncommute, φ, ψ, χ)
                QQ[:ho_nm][i,j,l] = hov(:nonmove, φ, ψ, χ)
            end
        end
    end
    mx = maximum(W); Wn = exp.(W .- mx)
    dφ = phi_free ? 1.0 / (nf - 1) : 1.0
    dlψ = (log(FPSI[end]) - log(FPSI[1])) / (np - 1)
    dχ = 1.0 / (nc - 1)
    (; W = Wn ./ sum(Wn), Q = QQ, logZ = mx + log(sum(Wn) * dφ * dlψ * dχ))
end
function wq(v, w, p)
    o = sortperm(vec(v)); vs = vec(v)[o]; cw = cumsum(vec(w)[o]); cw ./= cw[end]
    vs[searchsortedfirst(cw, p)]
end
wmean(v, w) = sum(vec(v) .* vec(w))

# canon = free φ (main) · phi0 = φ ≡ 0 (for log BF comparison)
const VARIANTS = [("canon", true), ("phi0", false)]
const CANON = "canon"
post = Dict{String,Any}()
psum = NamedTuple[]; pmarg = NamedTuple[]; mcomp = NamedTuple[]
for (nm, pf) in VARIANTS
    r = posterior(pf)
    post[nm] = r
    @printf("  %-8s logZ = %+10.3f\n", nm, r.logZ)
    for (k, lab) in ((:phi,"phi"), (:psi,"psi"), (:chi,"chi"), (:we_c10,"weekend_c10"),
                     (:we_c,"weekend_commute"), (:we_nc,"weekend_noncommute"), (:we_nm,"weekend_nonmove"),
                     (:ho_c,"holiday_commute"), (:ho_nc,"holiday_noncommute"), (:ho_nm,"holiday_nonmove"))
        v = r.Q[k]
        push!(psum, (variant = nm, quantity = lab, mean = mv(wmean(v, r.W)), q025 = mv(wq(v, r.W, 0.025)),
                     median = mv(wq(v, r.W, 0.5)), q975 = mv(wq(v, r.W, 0.975))))
    end
    for (k, lab) in ((:phi,"phi"), (:psi,"psi"), (:chi,"chi"))
        v = vec(r.Q[k]); w = vec(r.W); lv = sort(unique(v))
        length(lv) == 1 && continue
        d = Dict(x => 0.0 for x in lv); for (x, ww) in zip(v, w); d[x] += ww; end
        for x in lv; push!(pmarg, (variant = nm, param = lab, value = x, weight = d[x])); end
    end
end
CSV.write(joinpath(DATA, "posterior_summary.csv"), DataFrame(psum))
CSV.write(joinpath(DATA, "posterior_marginals.csv"), DataFrame(pmarg))
let prow = NamedTuple[]
    for (lab, ax, d) in (("phi", FPHI, Uniform(0, 1)), ("psi", FPSI, LogNormal(PR_PSI...)), ("chi", FCHI, PR_CHI))
        Δ = [i == 1 ? ax[2] - ax[1] : i == length(ax) ? ax[end] - ax[end-1] : (ax[i+1] - ax[i-1]) / 2 for i in eachindex(ax)]
        w = pdf.(d, ax) .* Δ; w ./= sum(w)
        for (x, ww) in zip(ax, w); push!(prow, (param = lab, value = x, weight = ww)); end
    end
    CSV.write(joinpath(DATA, "prior_marginals.csv"), DataFrame(prow))
end
"Minimum of the residual sum of squares J (only φ varies)."
function mindist(phi_free::Bool)
    best = (J = Inf, φ = 0.0)
    for φ in (phi_free ? FPHI : [0.0])
        J = 0.0
        for (k, bd) in enumerate(QSYM); J += ((wev(bd, φ) - wdref[bd] - T2_OBS[k]) / T2_SD)^2; end
        J < best.J && (best = (J = J, φ = φ))
    end
    best
end
for (nm, pf) in VARIANTS
    m = mindist(pf)
    ntg = length(QSYM); npar = pf ? 1 : 0; dfree = ntg - npar
    pv = ccdf(Chisq(dfree), m.J)
    push!(mcomp, (variant = nm, n_target = ntg, n_param = npar, df = dfree, J = mv(m.J), p_value = mv(pv),
                  logZ = mv(post[nm].logZ), md_phi = mv(m.φ)))
    @printf("  %-8s J = %8.2f  df %2d  p = %.3g  | φ %.3f\n", nm, m.J, dfree, pv, m.φ)
end
CSV.write(joinpath(DATA, "model_comparison.csv"), DataFrame(mcomp))
@printf("  log BF(free φ : φ≡0) = %+.2f\n", post["canon"].logZ - post["phi0"].logZ)
med(nm, q) = only(filter(x -> x.variant == nm && x.quantity == q, psum)).median
const TH = (φ = med(CANON, "phi"), ψ = med(CANON, "psi"), χ = med(CANON, "chi"))
@printf("\n  θ̂ = (φ %.4f, ψ %.4f, χ %.4f)   ← applied to all 72 months\n", TH.φ, TH.ψ, TH.χ)
# target fit — model values at θ̂ and residuals z = (g_k − T_k)/σ_k (table only)
let tfit = NamedTuple[]
    # T1 — validation (outside the likelihood). Also record weekday 10+ commute (weekday side of H1)
    push!(tfit, (target = "T1_weekend_c10", role = "held_out", observed = mv(T1_OBS), sigma = T1_SD,
                 predicted = mv(wev(:c10, TH.φ)), z = mv((wev(:c10, TH.φ) - T1_OBS) / T1_SD)))
    push!(tfit, (target = "H1_weekday_c10", role = "held_out", observed = NaN, sigma = NaN,
                 predicted = mv(SWD.c10), z = NaN))
    for (k, bd) in enumerate(QSYM)
        dq = wev(bd, TH.φ) - wdref[bd]
        push!(tfit, (target = "T2_" * TUS_BANDS[k], role = "target", observed = mv(T2_OBS[k]), sigma = T2_SD, predicted = mv(dq),
                     z = mv((dq - T2_OBS[k]) / T2_SD)))
    end
    CSV.write(joinpath(DATA, "target_fit.csv"), DataFrame(tfit))
    @printf("  target fit (θ̂): max |z| = %.2f · validation T1 %.4f vs model %.4f\n",
            maximum(abs(r.z) for r in tfit if r.role == "target"), T1_OBS, wev(:c10, TH.φ))
end
# χ prior sensitivity — holiday components (data do not update χ, so the prior sets the output uncertainty)
let crow = NamedTuple[]
    for (lab, pr) in (("Beta(2,2)", Beta(2.0, 2.0)), ("U(0,1)", Beta(1.0, 1.0)), ("Beta(2,3)", Beta(2.0, 3.0)), ("Beta(3,2)", Beta(3.0, 2.0)))
        r = posterior(true; chi_prior = pr)
        # point estimate = holiday solved at that prior's posterior median θ̂ (same definition as the released OD).
        #   median = posterior median of each quantity itself (reference) · q025/q975 = its posterior quantiles.
        θp = (φ = wq(r.Q[:phi], r.W, .5), ψ = wq(r.Q[:psi], r.W, .5), χ = wq(r.Q[:chi], r.W, .5))
        sθ = summarise(C9, solve(C9, "holiday", CALSEED; φ = θp.φ, ψ = θp.ψ, χ = θp.χ))
        pt = Dict(:chi => θp.χ, :ho_c => sθ.c, :ho_nc => sθ.nc, :ho_nm => sθ.nm)
        for (k, q) in ((:chi, "chi"), (:ho_c, "holiday_commute"), (:ho_nc, "holiday_noncommute"), (:ho_nm, "holiday_nonmove"))
            sc = q == "chi" ? 1.0 : 100 / PTOT
            push!(crow, (prior = lab, quantity = q, point_at_theta_hat = round(sc * pt[k], digits = 4),
                         median = round(sc * wq(r.Q[k], r.W, .5), digits = 4),
                         q025 = round(sc * wq(r.Q[k], r.W, .025), digits = 4), q975 = round(sc * wq(r.Q[k], r.W, .975), digits = 4)))
        end
    end
    CSV.write(joinpath(DATA, "chi_prior_sensitivity.csv"), DataFrame(crow))
end
for lab in ("phi", "psi", "chi", "weekend_commute", "weekend_noncommute", "weekend_nonmove",
            "holiday_commute", "holiday_noncommute", "holiday_nonmove")
    r = only(filter(x -> x.variant == CANON && x.quantity == lab, psum))
    sc = lab in ("phi", "psi", "chi") ? 1.0 : 100 / PTOT
    @printf("    %-20s %9.3f  [%9.3f, %9.3f]\n", lab, sc*r.median, sc*r.q025, sc*r.q975)
end

# ── §7 2022-09 final outputs ─────────────────────────────────────────────────
hr("7. 2022-09 final outputs (θ̂, both seeds)")
inner(A) = sum(A[j, j] for j in 1:N)
share = NamedTuple[]; sggrow = NamedTuple[]; cnc = NamedTuple[]; ver = NamedTuple[]; dayr = NamedTuple[]
HM  = Dict((dt, sd, k) => zeros(N, N) for dt in DTS for sd in SEEDS for k in ("commute", "noncommute"))
HM3 = Dict((dt, a3, k) => zeros(N, N) for dt in DTS for a3 in AG3S for k in ("commute", "noncommute"))
for dt in DTS, seed in SEEDS
    r = solve_full(C9, dt, seed, TH)
    T = zeros(4)
    D3 = Dict(a => zeros(N) for a in AG3S); P3 = Dict(a => zeros(N) for a in AG3S); Dall = zeros(N)
    for x in r
        g = x.g; sp = x.sp
        c_row = vec(sum(sp.Chat, dims = 2)); n_row = vec(sum(sp.NChat, dims = 2))
        cs, ns, ms, rs = sum(c_row), sum(n_row), sum(sp.nm), SUMP[g]
        T .+= (cs, ns, ms, rs)
        push!(share, (daytype = dt, seed = seed, age_group = g, lambda = mv(x.λg), pi = mv(x.πg),
            commute_pct = round(100cs/rs, digits = 2), noncommute_pct = round(100ns/rs, digits = 2),
            nonmove_pct = round(100ms/rs, digits = 2), total_pct = round(100(cs+ns+ms)/rs, digits = 2),
            commute = mv(cs), noncommute = mv(ns), nonmove = mv(ms), registered = mv(rs)))
        push!(cnc, (daytype = dt, seed = seed, age_group = g, registered = mv(rs),
            commute = mv(cs), commute_inner = mv(inner(sp.Chat)), commute_outer = mv(cs - inner(sp.Chat)),
            noncommute = mv(ns), noncommute_inner = mv(inner(sp.NChat)),
            noncommute_outer = mv(ns - inner(sp.NChat)), nonmove = mv(ms)))
        for j in 1:N
            push!(sggrow, (daytype = dt, seed = seed, sgg_cd = S[j], age_group = g,
                commute_row = mv(c_row[j]), noncommute_row = mv(n_row[j]),
                nonmove = mv(sp.nm[j]), nonmove_raw = mv(sp.nm_raw[j]), registered = mv(P[g][j])))
        end
        HM[(dt, seed, "commute")] .+= sp.Chat; HM[(dt, seed, "noncommute")] .+= sp.NChat
        D3[AG3[g]] .+= x.dcol; P3[AG3[g]] .+= P[g]; Dall .+= x.dcol
        if seed == CALSEED
            HM3[(dt, AG3[g], "commute")] .+= sp.Chat; HM3[(dt, AG3[g], "noncommute")] .+= sp.NChat
            for i in 1:N
                push!(dayr, (daytype = dt, level = "age15", age = g, sgg_cd = S[i],
                             daytime = round(x.dcol[i], digits = 2), registered = round(P[g][i], digits = 1),
                             ratio = round(x.dcol[i] / P[g][i], digits = 5)))
            end
        end
    end
    if seed == CALSEED
        for a in AG3S, i in 1:N
            push!(dayr, (daytype = dt, level = "age3", age = a, sgg_cd = S[i], daytime = round(D3[a][i], digits = 2),
                         registered = round(P3[a][i], digits = 1), ratio = round(D3[a][i] / P3[a][i], digits = 5)))
        end
        Pall = sum(P[g] for g in AG)
        for i in 1:N
            push!(dayr, (daytype = dt, level = "all", age = "all", sgg_cd = S[i], daytime = round(Dall[i], digits = 2),
                         registered = round(Pall[i], digits = 1), ratio = round(Dall[i] / Pall[i], digits = 5)))
        end
    end
    push!(ver, (daytype = dt, seed = seed, max_dest_dev = maximum(x.dm for x in r),
                max_row_identity = maximum(x.id for x in r),
                clamp_lo = sum(x.sp.n_lo for x in r), clamp_hi = sum(x.sp.n_hi for x in r),
                commute_pct = round(100T[1]/T[4], digits = 2), noncommute_pct = round(100T[2]/T[4], digits = 2),
                nonmove_pct = round(100T[3]/T[4], digits = 2)))
    @printf("  %-8s %-10s C %6.2f · NC %6.2f · NM %6.2f\n", dt, seed, 100T[1]/T[4], 100T[2]/T[4], 100T[3]/T[4])
end
CSV.write(joinpath(DATA, "decomp_share_final.csv"), DataFrame(share))
CSV.write(joinpath(DATA, "decomp_share_long.csv"), DataFrame([
    (daytype = r.daytype, seed = r.seed, age_group = r.age_group, component = c, pct = p)
    for r in share for (c, p) in (("commute", r.commute_pct), ("noncommute", r.noncommute_pct), ("nonmove", r.nonmove_pct))]))
CSV.write(joinpath(DATA, "sgg_decomp_final.csv"), DataFrame(sggrow))
CSV.write(joinpath(DATA, "cnc_share_daytype.csv"), DataFrame(cnc))
CSV.write(joinpath(DATA, "map_daytime_ratio.csv"), DataFrame(dayr))

# ── §8 validation (2022-09) ────────────────────────────────────────────────────
hr("8. 2022-09 validation")
sh = DataFrame(share)
tot3 = combine(groupby(sh, [:daytype, :seed]), :commute => sum => :C, :noncommute => sum => :NC,
               :nonmove => sum => :NM, :registered => sum => :Pg)
tot3.identity_dev = abs.(tot3.C .+ tot3.NC .+ tot3.NM .- tot3.Pg)
push!(ver, (daytype = "all", seed = "all", max_dest_dev = NaN, max_row_identity = maximum(tot3.identity_dev),
            clamp_lo = -1, clamp_hi = -1, commute_pct = NaN, noncommute_pct = NaN, nonmove_pct = NaN))
let ex = Dict("weekend" => summarise(C9, solve(C9, "weekend", CALSEED; φ = TH.φ)),
              "holiday" => summarise(C9, solve(C9, "holiday", CALSEED; φ = TH.φ, ψ = TH.ψ, χ = TH.χ)))
    got(x, c) = c === :c10 ? x.c10 : c === :commute ? x.c : c === :noncommute ? x.nc :
                c === :nonmove ? x.nm : x.q[String(c)[3:end]]
    d = 0.0
    for (dt, f) in (("weekend", c -> wev(c, TH.φ)), ("holiday", c -> hov(c, TH.φ, TH.ψ, TH.χ))), c in SUMCOL
        d = max(d, abs(f(c) - got(ex[dt], c)) / max(abs(got(ex[dt], c)), 1e-6))
    end
    @printf("  grid interpolation vs exact solution (θ̂): max rel diff %.3g\n", d)
    push!(ver, (daytype = "interp", seed = CALSEED, max_dest_dev = d, max_row_identity = NaN, clamp_lo = -1,
                clamp_hi = -1, commute_pct = NaN, noncommute_pct = NaN, nonmove_pct = NaN))
end
CSV.write(joinpath(DATA, "verification.csv"), DataFrame(ver))
for r in eachrow(tot3); @printf("  %-8s %-10s Σ(3 components) − ΣP = %.4g persons\n", r.daytype, r.seed, r.identity_dev); end

# seed comparison — same θ̂ with both initial guesses
seedcmp = NamedTuple[]
sg = DataFrame(sggrow)
for dt in DTS
    a = only(filter(r -> r.daytype == dt && r.seed == CALSEED, eachrow(tot3)))
    b = only(filter(r -> r.daytype == dt && r.seed == "mobile", eachrow(tot3)))
    for (mt, xa, xb) in (("commute_pct", a.C, b.C), ("noncommute_pct", a.NC, b.NC), ("nonmove_pct", a.NM, b.NM))
        push!(seedcmp, (daytype = dt, level = "national", scope = "all", metric = mt,
                        radiation = round(100xa/a.Pg, digits = 3), mobile = round(100xb/b.Pg, digits = 3),
                        diff_pp = round(100(xa - xb)/a.Pg, digits = 3)))
    end
    for (mt, col) in (("commute_pct", :commute_row), ("noncommute_pct", :noncommute_row), ("nonmove_pct", :nonmove))
        agg(seed) = sort(combine(groupby(sg[(sg.daytype .== dt) .& (sg.seed .== seed), :], :sgg_cd),
                                 col => sum => :v, :registered => sum => :p), :sgg_cd)
        A = agg(CALSEED); B = agg("mobile"); va = 100 .* A.v ./ A.p; vb = 100 .* B.v ./ B.p
        push!(seedcmp, (daytype = dt, level = "adm2", scope = "mae_pp", metric = mt, radiation = round(mean(va), digits = 3),
                        mobile = round(mean(vb), digits = 3), diff_pp = round(mean(abs.(va .- vb)), digits = 3)))
        push!(seedcmp, (daytype = dt, level = "adm2", scope = "cor", metric = mt, radiation = NaN, mobile = NaN,
                        diff_pp = round(cor(va, vb), digits = 5)))
    end
    for k in ("commute", "noncommute")
        A = HM[(dt, CALSEED, k)]; B = HM[(dt, "mobile", k)]
        push!(seedcmp, (daytype = dt, level = "od_cell", scope = "cor_log1p", metric = k,
                        radiation = round(sum(A)/1e6, digits = 4), mobile = round(sum(B)/1e6, digits = 4),
                        diff_pp = round(cor(log1p.(vec(A)), log1p.(vec(B))), digits = 5)))
        push!(seedcmp, (daytype = dt, level = "od_cell", scope = "share_abs_diff", metric = k,
                        radiation = round(sum(A)/1e6, digits = 4), mobile = round(sum(B)/1e6, digits = 4),
                        diff_pp = round(sum(abs.(A .- B)) / sum(A), digits = 5)))
    end
end
CSV.write(joinpath(DATA, "seed_comparison.csv"), DataFrame(seedcmp))
for r in seedcmp
    r.level == "national" && @printf("  seed diff %-8s %-15s radiation %6.2f · mobile %6.2f · diff %+.2f%%p\n",
                                     r.daytype, r.metric, r.radiation, r.mobile, r.diff_pp)
end

# Visitor-data held-out check — 2022-09 district data (ratios) + province data (JSD)
println("\n  Korea Tourism Data Lab visitor data — held-out")
const ECONF = joinpath(INP, "econ_active_rate_adm1_agegroup_2020.csv")
const A1 = sort(unique(adm1.(S)))
const TOSIDO = sido_code_map(ECONF)
const ADM1NM = load_adm1_names(ECONF)
nco(dt) = sum(HM[(dt, CALSEED, "noncommute")]) - inner(HM[(dt, CALSEED, "noncommute")])
nc_dest(dt, sd) = permutedims(HM[(dt, sd, "noncommute")])
const KTDB_DEST = permutedims(VALL)
let tocode = visitor_code_map(joinpath(INP, "kosis_registered_pop_sgg_sexage_202209.csv"), ECONF)
    TOUR, VND, VTOT, VDIAG = read_visitor(joinpath(INP, "visitor_202209.csv"), tocode, IDX, N)
    for dt in DTS; VND[dt] == C9.comp_days[dt].n || error("visitor data $dt day count ≠ calendar"); end
    vdaily = Dict(dt => VTOT[dt] / VND[dt] for dt in DTS)
    vrows = [(daytype = dt, n_days = VND[dt], total_visits = round(VTOT[dt]), daily_mean = round(vdaily[dt]),
              ratio_vs_weekday = round(vdaily[dt] / vdaily["weekday"], digits = 4),
              offdiag_zero_pct = round(100count(TOUR[dt][i, j] == 0.0 for i in 1:N, j in 1:N if i != j) / (N*(N-1)), digits = 1),
              model_nc_outer = mv(nco(dt))) for dt in DTS]
    CSV.write(joinpath(DATA, "visitor_daytype_totals.csv"), DataFrame(vrows))
    TOUR1, VND1, _, VDIAG1 = read_visitor_sido(joinpath(INP, "visitor_sido_202209.csv"), TOSIDO, A1)
    @assert VDIAG1 == 0.0
    a2 = Dict(dt => sum(TOUR[dt]) for dt in DTS); a1 = Dict(dt => sum(TOUR1[dt]) for dt in DTS)
    CSV.write(joinpath(DATA, "visitor_scale_ratio.csv"), DataFrame([
        (daytype = dt, adm2_daily = mv(a2[dt]), adm1_daily = mv(a1[dt]),
         adm2_ratio_vs_weekend = round(a2[dt] / a2["weekend"], digits = 4),
         adm1_ratio_vs_weekend = round(a1[dt] / a1["weekend"], digits = 4)) for dt in DTS]))
    @printf("    holiday/weekend: visitors district %.4f · province %.4f · expressway %.4f · ψ̂ %.4f · model NC outer %.4f\n",
            a2["holiday"]/a2["weekend"], a1["holiday"]/a1["weekend"], PSI_MAX, TH.ψ, nco("holiday")/nco("weekend"))
    jrows = NamedTuple[]; jcell1 = NamedTuple[]; hr1 = NamedTuple[]
    for dt in DTS
        w1 = vec(sum(TOUR1[dt], dims = 2))
        t1 = [("NC_" * CALSEED, to_adm1(nc_dest(dt, CALSEED), S, A1)), ("NC_mobile", to_adm1(nc_dest(dt, "mobile"), S, A1)),
              ("KTDB_V_seed", to_adm1(KTDB_DEST, S, A1))]
        for other in DTS; other == dt || push!(t1, ("tourism_" * other, TOUR1[other])); end
        for (lab, M1) in t1
            v, ok = jsd_rows(zero_diag_rownorm(TOUR1[dt]), zero_diag_rownorm(M1))
            s1 = jsd_summary(v, ok, w1)
            push!(jrows, (daytype = dt, target = lab, origins = "native", scale = "adm1", n_row = s1.n, mean = mv(s1.mean),
                          wmean = mv(s1.wmean), median = mv(s1.median), q25 = mv(s1.q25), q75 = mv(s1.q75), max = mv(s1.max)))
            lab == "NC_" * CALSEED || continue
            for i in eachindex(A1)
                ok[i] || continue
                push!(jcell1, (daytype = dt, adm1_code = A1[i], adm1_name = get(ADM1NM, A1[i], A1[i]),
                               jsd = round(v[i], digits = 5), visits = round(w1[i], digits = 1),
                               visit_share = round(w1[i] / sum(w1), digits = 5)))
            end
        end
        P1 = zero_diag_rownorm(TOUR1[dt]); Q1 = zero_diag_rownorm(to_adm1(nc_dest(dt, CALSEED), S, A1))
        for (src, M) in (("tourism", P1), ("model_NC", Q1)), i in eachindex(A1), j in eachindex(A1)
            push!(hr1, (daytype = dt, source = src, visit_code = A1[i], visit_name = get(ADM1NM, A1[i], A1[i]),
                        res_code = A1[j], res_name = get(ADM1NM, A1[j], A1[j]),
                        prob = isnan(M[i, j]) ? 0.0 : round(M[i, j], digits = 6)))
        end
    end
    CSV.write(joinpath(DATA, "visitor_jsd.csv"), DataFrame(jrows))
    CSV.write(joinpath(DATA, "visitor_jsd_rows_adm1.csv"), DataFrame(jcell1))
    CSV.write(joinpath(DATA, "visitor_heatmap_adm1.csv"), DataFrame(hr1))
    for r in jrows
        r.target in ("NC_" * CALSEED, "KTDB_V_seed") &&
            @printf("    JSD %-8s %-14s median %.4f · IQR %.4f\n", r.daytype, r.target, r.median, r.q75 - r.q25)
    end
end

# ── §9 figure summaries (2022-09) ─────────────────────────────────────────────
hr("9. 2022-09 figure summaries")
CSV.write(joinpath(DATA, "adm1_ticks.csv"), DataFrame([
    let ks = findall(j -> adm1(S[j]) == a, 1:N)
        (adm1_code = a, adm1_name = get(ADM1NM, a, a), i_min = minimum(ks), i_max = maximum(ks),
         i_mid = (minimum(ks) + maximum(ks)) / 2)
    end for a in A1]))
let rows = NamedTuple[]
    for dt in DTS, a3 in AG3S, k in ("commute", "noncommute"), j in 1:N, i in 1:N
        push!(rows, (daytype = dt, age3 = a3, kind = k, from_cd = S[j], to_cd = S[i], value = mv(HM3[(dt, a3, k)][j, i])))
    end
    CSV.write(joinpath(DATA, "heatmap_agegroup3_daytype.csv"), DataFrame(rows))
end
HM3 = nothing; GC.gc()
let pg = CSV.read(joinpath(INP, "sgg_polygons.csv"), DataFrame;
                  types = Dict(:sgg_cd=>String, :gid=>String, :ord=>Int, :x=>Float64, :y=>Float64))
    keep = falses(nrow(pg))
    for sub in groupby(pg, :gid)
        idx = parentindices(sub)[1]
        for (k, r) in enumerate(idx); (isodd(k) || k == length(idx)) && (keep[r] = true); end
    end
    pg = pg[keep, :]; pg.adm1 = adm1.(pg.sgg_cd)
    bigpart = Set{String}()
    for sub in groupby(pg, :gid)
        hypot(maximum(sub.x) - minimum(sub.x), maximum(sub.y) - minimum(sub.y)) >= 8_000 && push!(bigpart, first(sub.gid))
    end
    pg.cap_outline = in.(pg.adm1, Ref(["11","28","41"])) .& in.(pg.gid, Ref(bigpart))
    CSV.write(joinpath(DATA, "map_sgg_polygons.csv"), pg)
    cap = pg[in.(pg.adm1, Ref(["11","28","41"])), :]
    far = Set(g for g in unique(cap.gid) if startswith(g, "28") && maximum(cap[cap.gid .== g, :x]) < 885_000)
    cap = cap[.!in.(cap.gid, Ref(far)), :]
    nat = pg[.!in.(pg.sgg_cd, Ref(["47940","28720"])), :]
    box(d, m) = (minimum(d.x) - m*(maximum(d.x)-minimum(d.x)), maximum(d.x) + m*(maximum(d.x)-minimum(d.x)),
                 minimum(d.y) - m*(maximum(d.y)-minimum(d.y)), maximum(d.y) + m*(maximum(d.y)-minimum(d.y)))
    CSV.write(joinpath(DATA, "map_extent.csv"), DataFrame([(scope = s, xmin = b[1], xmax = b[2], ymin = b[3], ymax = b[4])
                                                           for (s, b) in (("national", box(nat, 0.02)), ("capital", box(cap, 0.03)))]))
end
let dr = DataFrame(dayr)
    dr = dr[dr.level .== "all", :]
    rt = Dict((r.daytype, r.sgg_cd) => r.ratio for r in eachrow(dr))
    rows = [(sgg_cd = S[i], sgg_nm = get(SGGNM, S[i], ""), adm1_code = adm1(S[i]), adm1_name = get(ADM1NM, adm1(S[i]), ""),
             registered = mv(sum(P[g][i] for g in AG)), ratio_weekday = rt[("weekday", S[i])],
             ratio_weekend = rt[("weekend", S[i])], ratio_holiday = rt[("holiday", S[i])],
             diff_weekend_weekday = round(rt[("weekend", S[i])] - rt[("weekday", S[i])], digits = 5),
             diff_holiday_weekday = round(rt[("holiday", S[i])] - rt[("weekday", S[i])], digits = 5)) for i in 1:N]
    sort!(rows, by = r -> r.diff_holiday_weekday)
    CSV.write(joinpath(DATA, "sgg_daytime_change.csv"), DataFrame(rows))
    iqr(v) = quantile(v, .75) - quantile(v, .25)
    @printf("  presence ratio IQR across districts — weekday %.3f · weekend %.3f · holiday %.3f\n",
            iqr([r.ratio_weekday for r in rows]), iqr([r.ratio_weekend for r in rows]), iqr([r.ratio_holiday for r in rows]))
end
let sgx = copy(sg), sc = NamedTuple[]
    sgx.age3 = [AG3[g] for g in sgx.age_group]
    for grp in groupby(sgx, [:daytype, :seed, :sgg_cd])
        for (lab, sub) in (("all", grp), ("20-64", grp[grp.age3 .== "20-64", :]))
            P0 = sum(sub.registered); P0 > 0 || continue
            push!(sc, (daytype = grp.daytype[1], seed = grp.seed[1], sgg_cd = grp.sgg_cd[1], scope = lab, registered = mv(P0),
                       commute_pct = round(100sum(sub.commute_row) / P0, digits = 3),
                       noncommute_pct = round(100sum(sub.noncommute_row) / P0, digits = 3),
                       nonmove_pct = round(100sum(sub.nonmove) / P0, digits = 3)))
        end
    end
    CSV.write(joinpath(DATA, "sgg_scatter.csv"), DataFrame(sc))
end
HM = nothing; sg = nothing; GC.gc()
@printf("\n  Part B done — elapsed %.1f min\n", (time() - T_START) / 60)

# ═════════════════════════════════════════════════════════════════════════════
# Part C — 72 months (radiation seed · θ̂)
# ═════════════════════════════════════════════════════════════════════════════
hr("10. Monthly OD — 2018-01 … 2023-12 (radiation · θ̂ = 2022-09 posterior median)")

"Write one OD set as gzip csv — columns: kind (commute/noncommute/nonmove) · from_cd (origin j) · to_cd (destination i) · age_group · value."
function write_od_gz(path, r)
    n = sum(count(>(0), x.sp.Chat) + count(>(0), x.sp.NChat) + count(>(0), x.sp.nm) for x in r)
    kind = Vector{String}(undef, n); fr = Vector{String}(undef, n); to = Vector{String}(undef, n)
    ag = Vector{String}(undef, n); val = Vector{Float64}(undef, n)
    k = 0
    for x in r
        for (kd, A) in (("commute", x.sp.Chat), ("noncommute", x.sp.NChat)), j in 1:N, i in 1:N
            v = A[j, i]; v > 0 || continue
            k += 1; kind[k] = kd; fr[k] = S[j]; to[k] = S[i]; ag[k] = x.g; val[k] = mv(v)
        end
        for j in 1:N
            v = x.sp.nm[j]; v > 0 || continue
            k += 1; kind[k] = "nonmove"; fr[k] = S[j]; to[k] = S[j]; ag[k] = x.g; val[k] = mv(v)
        end
    end
    CSV.write(path, DataFrame(kind = kind, from_cd = fr, to_cd = to, age_group = ag, value = val); compress = true)
    n
end

NAT = NamedTuple[]; VERM = NamedTuple[]
t_loop = time()
for (y, m) in YM_ALL
    tm = time()
    c = (y, m) == (Y0, M0) ? C9 : mkctx(y, m, nothing)
    ym = ymlabel(c)
    for dt in c.dts
        r = solve_full(c, dt, CALSEED, TH)
        nd = c.comp_days[dt]
        C = zeros(N, N); NC = zeros(N, N); NMv = zeros(N)
        for x in r
            C .+= x.sp.Chat; NC .+= x.sp.NChat; NMv .+= x.sp.nm
        end
        push!(NAT, (ym = ym, year = y, month = m, daytype = dt, n_days = nd.n, kind = nd.kind,
                    registered = mv(c.PTOT), commute = mv(sum(C)), commute_inner = mv(inner(C)),
                    noncommute = mv(sum(NC)), noncommute_inner = mv(inner(NC)), nonmove = mv(sum(NMv))))
        push!(VERM, (ym = ym, daytype = dt, max_row_identity = maximum(x.id for x in r),
                     max_dest_dev = maximum(x.dm for x in r), clamp_lo = sum(x.sp.n_lo for x in r),
                     clamp_hi = sum(x.sp.n_hi for x in r), nadd = c.nadd,
                     neg_offdiag = sum(zero_breakdown(x.sp.Chat, N).n_neg + zero_breakdown(x.sp.NChat, N).n_neg for x in r)))
        write_od_gz(joinpath(ODDIR, "od_radiation_$(ym)_$(dt).csv.gz"), r)
    end
    @printf("  %s  %-24s  %4.1f s  (cumulative %.1f min)\n", ym, join(c.dts, ","), time() - tm, (time() - t_loop) / 60)
    flush(stdout)
    c = nothing
    (m == 12) && GC.gc()
end
natdf = DataFrame(NAT)
natdf.commute_outer = natdf.commute .- natdf.commute_inner
natdf.noncommute_outer = natdf.noncommute .- natdf.noncommute_inner
for k in (:commute, :noncommute, :nonmove, :commute_outer, :noncommute_outer)
    natdf[!, Symbol(string(k) * "_pct")] = round.(100 .* natdf[!, k] ./ natdf.registered, digits = 4)
end
CSV.write(joinpath(DATA, "month_national.csv"), natdf)
CSV.write(joinpath(DATA, "month_verification.csv"), DataFrame(VERM))
@printf("  72 months done — %.1f min · %d OD files\n", (time() - t_loop) / 60, count(endswith(".csv.gz"), readdir(ODDIR)))
let v = DataFrame(VERM)
    @printf("  monthly checks: max row identity %.3g persons · max destination marginal %.3g persons · negative off-diagonal %d\n",
            maximum(v.max_row_identity), maximum(v.max_dest_dev), sum(v.neg_offdiag))
end

# ── §11 weekday three components (all ages · not age-standardised) ────────
hr("11. Monthly time series of the weekday three components (all ages, % of registered population)")
# Key social-distancing dates (figure annotation only, not used in the analysis). kind = tighten · relax.
# Dates and names were checked against Korean government policy briefings (korea.kr).
CSV.write(joinpath(DATA, "covid_events.csv"), DataFrame(
    date  = ["2020-02-23", "2020-08-30", "2020-12-08", "2021-07-12", "2021-11-01", "2022-04-18"],
    label = ["Crisis alert raised to Serious", "Capital region enhanced Level 2", "Third wave · capital region Level 2.5",
             "Capital region Level 4", "Phased return to normal", "Lifting of social distancing"],
    kind  = ["tighten", "tighten", "tighten", "tighten", "relax", "relax"]))
let wd = natdf[natdf.daytype .== "weekday", :]
    rows = [(ym = r.ym, year = r.year, month = r.month, component = k, pct = r[Symbol(k * "_pct")])
            for r in eachrow(wd) for k in ("commute", "noncommute", "nonmove")]
    wc = DataFrame(rows)
    CSV.write(joinpath(DATA, "weekday_components_crude.csv"), wc)
    # Summary — simple monthly means. Windows are defined by whether the **whole month** lies in the grey band (2020-02-23 – 2022-04-18):
    #   before = 2018-01 – 2020-01 · during = 2020-03 – 2022-03 · after = 2022-05 – 2023-12 (boundary months 2020-02 · 2022-04 excluded)
    win(y, m) = (y, m) <= (2020, 1) ? "before" : (2020, 3) <= (y, m) <= (2022, 3) ? "during" :
                (y, m) >= (2022, 5) ? "after" : "edge"
    srow = NamedTuple[]
    for k in ("commute", "noncommute", "nonmove")
        s = wc[wc.component .== k, :]
        w = win.(s.year, s.month)
        for (lab, sel) in (("year_2018", s.year .== 2018), ("year_2023", s.year .== 2023),
                           ("before", w .== "before"), ("during", w .== "during"), ("after", w .== "after"))
            push!(srow, (component = k, window = lab, n_months = count(sel), mean_pct = round(mean(s.pct[sel]), digits = 3),
                         min_pct = round(minimum(s.pct[sel]), digits = 3), max_pct = round(maximum(s.pct[sel]), digits = 3)))
        end
        i = argmin(s.pct)
        @printf("  %-11s 2018 %6.2f → 2023 %6.2f  | before %6.2f · during %6.2f · after %6.2f | min %s %.2f\n", k,
                (only(r.mean_pct for r in srow if r.component == k && r.window == x) for x in
                 ("year_2018", "year_2023", "before", "during", "after"))..., s.ym[i], s.pct[i])
    end
    CSV.write(joinpath(DATA, "weekday_components_summary.csv"), DataFrame(srow))
end

# ── §11b district code table for release ─────────────────────────────────────────────────
# from_cd · to_cd in the OD files are stable codes, fixed over 2018–2023, i.e. the pre-change codes:
# Michuhol-gu 28170, Gunwi 47720, Gangwon 42xxx, Jeonbuk 45xxx (35 districts differ from the December 2023 codes).
# The 2023-12 administrative code is listed alongside to link to the current code system.
hr("11b. District code table (district_codes.csv for release)")
let reg = CSV.read(joinpath(INP, "kosis_registered_pop_sgg_sexage_2018-2023.csv"), DataFrame;
                   select = [:year, :month, :sgg_cd, :sgg_cd_stable],
                   types = Dict(:sgg_cd => String, :sgg_cd_stable => String))
    cur = unique(reg[(reg.year .== 2023) .& (reg.month .== 12), [:sgg_cd_stable, :sgg_cd]])
    @assert nrow(cur) == N && length(unique(cur.sgg_cd_stable)) == N
    curd = Dict(r.sgg_cd_stable => r.sgg_cd for r in eachrow(cur))
    ADM1_EN = Dict("11"=>"Seoul","26"=>"Busan","27"=>"Daegu","28"=>"Incheon","29"=>"Gwangju","30"=>"Daejeon",
                   "31"=>"Ulsan","36"=>"Sejong","41"=>"Gyeonggi","42"=>"Gangwon","43"=>"Chungbuk",
                   "44"=>"Chungnam","45"=>"Jeonbuk","46"=>"Jeonnam","47"=>"Gyeongbuk","48"=>"Gyeongnam","50"=>"Jeju")
    dc = DataFrame([(district_code = S[i], district_name_ko = SGGNM[S[i]], province_code = adm1(S[i]),
                     province_name_ko = ADM1NM[adm1(S[i])], province_name_en = ADM1_EN[adm1(S[i])],
                     code_2023_12 = curd[S[i]], matrix_index = i) for i in 1:N])
    @assert all(!isempty, dc.district_name_ko) && length(unique(dc.province_code)) == 17
    CSV.write(joinpath(DATA, "district_codes.csv"), dc)
    @printf("  %d districts · %d provinces · %d districts with a different 2023-12 code\n", N, length(unique(dc.province_code)),
            count(dc.district_code .!= dc.code_2023_12))
end

# ═════════════════════════════════════════════════════════════════════════════
hr("Outputs")
for fn in sort(readdir(DATA))
    endswith(fn, ".csv") || continue
    p = joinpath(DATA, fn)
    @printf("  %-40s %9.1f MB  sha16 %s\n", fn, filesize(p) / 1e6, sha16(p))
end
@printf("  od/ : %d files · %.2f GB\n", length(readdir(ODDIR)), sum(filesize(joinpath(ODDIR, f)) for f in readdir(ODDIR)) / 1e9)
println("\n", "="^78)
@printf("Done — %.1f min · θ̂ = (φ %.4f, ψ %.4f, χ %.4f) applied to 72 months\n", (time() - T_START) / 60, TH.φ, TH.ψ, TH.χ)
println("="^78)

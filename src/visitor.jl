# src/visitor.jl — **held-out comparison** with the Korea Tourism Data Lab daily visitor data for 2022-09.
#
# ⚠️ **Used for validation, not as input.** This is itself a paid mobile-phone-based OD produced by the
#    agency, so it enters no seed, marginal or likelihood. Doing so would break the premise of this study
#    (see Methods): *"building population movement without paid mobile-phone OD"*. Used only in the comparison of run.jl §8.
#
# ── Data (input/data_descript.png) ──────────────────────────────────────────
#   날짜 (date) YYYYMMDD · day_type **weekday(1) · weekend(2) · holiday(3)** · 방문지_시군구 (visited district, 250) ·
#   거주지_시군구 (home district, 250) · 비율 (share of the visited district's visitors by home district, %) ·
#   방문자수 (visitors = total visitors to the district × share)
#
#   ★ The source data already uses the names `weekday`·`weekend`·`holiday` — our day-type names
#     match its terms exactly. Its date composition also matches our calendar (checked in run.jl §8).
#
# "Visitor" definition: commuters (work or school) **excluded** · stay of at least **30 min** ·
# a person visiting several districts is **counted multiple times**. The unit is **counts** of (person × visited district),
# so **levels cannot be compared** with our $\widehat{NC}$ (persons). We therefore use only two things:
#   (A) day-type **ratio** — holiday / weekend (second line of evidence for ψ)
#   (B) **spatial structure after row normalization** — JS divergence (run.jl §8)
# The diagonal (home = visited) is **exactly 0** in the data → only cross-district visits are counted.

const VIS_REN = Dict("강원특별자치도" => "강원도", "전북특별자치도" => "전라북도")
const VIS_DT  = Dict(1 => "weekday", 2 => "weekend", 3 => "holiday")

"""
    visitor_code_map(kosis_path, econ_path)

Visitor data `"province district"` names → our adm2 code. All 250 match 1:1.
  * province renamings (Gangwon · Jeonbuk Special Self-Governing Province) are reverted
  * non-autonomous districts are matched by their last token (`경기도 고양시 덕양구` → `덕양구`)
  * `세종특별자치시` has province = district; `군위군` follows the boundary **before** its 2023 transfer to Daegu, i.e. Gyeongsangbuk-do
"""
function visitor_code_map(kosis_path, econ_path)
    a1 = Dict{String,String}()
    for r in CSV.Rows(econ_path; types = String)
        a1[r.adm1_code] = r.adm1_name
    end
    code = Dict{Tuple{String,String},String}()
    for r in CSV.Rows(kosis_path; types = String)
        code[(a1[r.adm1_code], r.sgg_nm)] = r.sgg_cd_stable
    end
    function tocode(nm::AbstractString)
        nm == "세종특별자치시" && return code[("세종특별자치시", "세종특별자치시")]
        p = split(nm, ' ')
        last(p) == "군위군" && return code[("경상북도", "군위군")]
        code[(get(VIS_REN, p[1], p[1]), last(p))]
    end
    tocode
end

"""
    read_visitor(path, tocode, IDX, N)

Read the **daily-mean** visitor OD by day type as `A[i, j]` = (visited i, home j).
⚠️ Only this matrix has **row = destination** — the data's `비율` is defined per visited district, so we keep that axis
(all other pipeline matrices are `A[j,i]`, row = origin). It is named `TOUR` to avoid confusion.

Returns `(TOUR, ndays, total, ndiag)` — day type => 250×250 daily-mean matrix · number of days ·
total count (not daily mean) · diagonal sum (should be 0).
"""
function read_visitor(path, tocode, IDX, N)
    TOUR  = Dict(dt => zeros(N, N) for dt in values(VIS_DT))
    days  = Dict(dt => Set{String}() for dt in values(VIS_DT))
    total = Dict(dt => 0.0 for dt in values(VIS_DT))
    ndiag = 0.0
    cmap  = Dict{String,Int}()
    ix(nm) = get!(cmap, nm) do; IDX[tocode(nm)] end
    for r in CSV.Rows(path; types = Dict(:day_type => Int, :방문자수 => Float64))
        dt = VIS_DT[r.day_type]
        push!(days[dt], r.날짜)
        v = r.방문자수
        total[dt] += v
        i = ix(r.방문지_시군구); j = ix(r.거주지_시군구)
        i == j && (ndiag += v)
        TOUR[dt][i, j] += v
    end
    nd = Dict(dt => length(days[dt]) for dt in values(VIS_DT))
    for dt in values(VIS_DT); TOUR[dt] ./= nd[dt]; end
    (TOUR, nd, total, ndiag)
end

# ═════════════════════════════════════════════════════════════════════════════
# JS divergence — compares the home-district distribution of each row (= visited district)
# ═════════════════════════════════════════════════════════════════════════════
# Convention: diagonal set to 0 → row normalization → base-2 JSD per row.
#
# **The only scale is the 17 adm1 units.** At adm2 250, 19–27% of the visitor data's off-diagonal
# cells are exactly 0, so that sparsity floor is shared by all candidates, and the ratio between the raw
# KTDB shape and our output shrinks from 1.42–3.22× at adm1 to 1.06–1.21× at adm2.
# **A metric that cannot discriminate is not used.**
#
# **Why JSD** — bounded in $[0,1]$ (base 2), so comparable across scales; symmetric, so neither side
# has to be treated as the truth; and finite even when supports differ (KL would diverge).

"""
    dead_origins(TOUR) -> Vector{Int}

Find **home districts (columns) whose sum is exactly 0 in every day type**. In the 2022-09 data only **Sejong**
is caught — normal as a visited district, but **as a home district all 30 days × 250 destinations are 0**.
Sejong residents cannot have gone nowhere, so this is a **gap in the data**, not a real signal.

Left as is, only our model puts mass in that cell in every row, inflating JSD **systematically**
(that one cell adds about `0.5q` per row). The comparison therefore drops that column from both distributions
before normalizing — i.e. compares the **distributions conditional on home ≠ Sejong**. `run.jl` §8 reports both
versions to show the size of this effect.
"""
function dead_origins(TOUR)
    ks = collect(keys(TOUR)); n = size(TOUR[first(ks)], 2)
    [j for j in 1:n if all(sum(@view TOUR[k][:, j]) == 0 for k in ks)]
end

"Set the diagonal to 0 and normalize each row to sum to 1. Rows summing to 0 are left as `NaN`."
function zero_diag_rownorm(M)
    X = copy(M)
    for i in 1:size(X, 1); X[i, i] = 0.0; end
    for i in 1:size(X, 1)
        s = sum(@view X[i, :])
        X[i, :] .= s > 0 ? (@view(X[i, :]) ./ s) : NaN
    end
    X
end

kl_term(p, m) = (p > 0 && m > 0) ? p * log2(p / m) : 0.0
"base-2 Jensen–Shannon divergence. Both probability vectors must sum to 1."
function jsd(p, q)
    s = 0.0
    @inbounds for k in eachindex(p)
        m = 0.5 * (p[k] + q[k])
        s += 0.5 * kl_term(p[k], m) + 0.5 * kl_term(q[k], m)
    end
    s
end

"""
    jsd_rows(P, Q) -> (values, ok)

JSD per row. Rows whose sum was 0 on either side (`NaN`) are skipped and flagged in `ok`.
"""
function jsd_rows(P, Q)
    n = size(P, 1)
    v = fill(NaN, n); ok = falses(n)
    for i in 1:n
        (isnan(P[i, 1]) || isnan(Q[i, 1])) && continue
        v[i] = jsd(@view(P[i, :]), @view(Q[i, :])); ok[i] = true
    end
    (v, ok)
end
"250×250 (row = visited, column = home) → 17×17 adm1 aggregation."
function to_adm1(M, S, A1)
    m = Dict(a => k for (k, a) in enumerate(A1))
    n = length(A1); X = zeros(n, n)
    for i in eachindex(S), j in eachindex(S)
        X[m[adm1(S[i])], m[adm1(S[j])]] += M[i, j]
    end
    X
end

"Weighted/unweighted summary (weight = the visit mass of the row)."
function jsd_summary(v, ok, w)
    x = v[ok]; ww = w[ok]
    isempty(x) && return (; n = 0, mean = NaN, wmean = NaN, median = NaN,
                          q25 = NaN, q75 = NaN, max = NaN)
    s = sort(x); m = length(s)
    q(p) = s[clamp(ceil(Int, p * m), 1, m)]
    (; n = m, mean = sum(x) / m, wmean = sum(x .* ww) / sum(ww),
       median = q(0.5), q25 = q(0.25), q75 = q(0.75), max = s[end])
end

# ═════════════════════════════════════════════════════════════════════════════
# Province (adm1) level source data
# ═════════════════════════════════════════════════════════════════════════════
# ⚠️ **Summing adm2 data to provinces does not give adm1 data.**
#    A "visitor" is a **count** of *person × visited spatial unit*, so one person touring three districts in Gangwon
#    is 3 counts in the adm2 table but 1 in the adm1 table. Summing therefore **inflates inter-province visits**
#    — empirically by weekday 1.20 · weekend 1.28 · **holiday 1.31×**, and since the factor grows with day type
#    (multi-destination trips are common on holidays) it introduces **systematic bias into day-type comparisons**.
#    Our $\widehat{NC}$ is a **flow** matrix of persons, so summing it is exact — that is the asymmetry.
#    The adm1 comparison therefore uses **province-level source data aggregated separately by the agency**.

"Province name → our adm1 code (2 digits). Reverts the renamings in the data (Gangwon · Jeonbuk Special Self-Governing Province)."
function sido_code_map(econ_path)
    rev = Dict{String,String}()
    for r in CSV.Rows(econ_path; types = String)
        rev[r.adm1_name] = r.adm1_code
    end
    nm -> rev[get(VIS_REN, nm, nm)]
end

"""
    read_visitor_sido(path, tosido, A1)

Read one monthly province-level file as **daily-mean** visit matrices by day type, `A[i, j]` = (visited i, home j).
Same axis convention as the adm2 version (`read_visitor`).

Two differences from the adm2 version:
  * **no diagonal rows at all** (the adm2 version had rows with value 0) — `zero_diag_rownorm` removes them anyway.
  * **no Sejong-as-home gap** (3.5 million counts over 30 days). The all-zero Sejong column in the adm2 table was a gap in that table.
"""
function read_visitor_sido(path, tosido, A1)
    n = length(A1); idx = Dict(a => k for (k, a) in enumerate(A1))
    TOUR1 = Dict(dt => zeros(n, n) for dt in values(VIS_DT))
    days  = Dict(dt => Set{String}() for dt in values(VIS_DT))
    total = Dict(dt => 0.0 for dt in values(VIS_DT))
    ndiag = 0.0
    for r in CSV.Rows(path; types = Dict(:day_type => Int, :방문자수 => Float64))
        dt = VIS_DT[r.day_type]
        push!(days[dt], r.날짜); total[dt] += r.방문자수
        i = idx[tosido(r.방문지_시도)]; j = idx[tosido(r.거주지_시도)]
        i == j && (ndiag += r.방문자수)
        TOUR1[dt][i, j] += r.방문자수
    end
    nd = Dict(dt => length(days[dt]) for dt in values(VIS_DT))
    for dt in values(VIS_DT); nd[dt] > 0 && (TOUR1[dt] ./= nd[dt]); end
    (TOUR1, nd, total, ndiag)
end

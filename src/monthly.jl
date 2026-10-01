# src/monthly.jl — **monthly inputs** for 2018-01 … 2023-12 (four files in `input/`).
#
#   | What | File | Rule |
#   |---|---|---|
#   | Calendar | `daytype_2018-2023.csv` | runs of 1–2 consecutive non-working days = weekend, 3+ days = holiday |
#   | Registered population $P^a_j$ | `kosis_registered_pop_sgg_sexage_2018-2023.csv` | summed over sex · summed into age groups |
#   | De facto day-type mean $\bar d^{a,t}_i(h)$ | `defacto_daytype_2018-2023.csv` | summed over sex → per date → day-type mean, 2 decimals |
#   | Census 3-way split | `sgg_commute_agegroup_2018-2023.csv` | census 2020 shares × that month's registered population, 3 decimals |
#
#   run.jl §2 checks that, for 2022-09, these rules reproduce the 2022-09 inputs (`marginals_202209.csv` etc.).
#   Matrix convention as in `src/inputs.jl`: A[j,i], row = origin j, column = destination i.

const YM_ALL = [(y, m) for y in 2018:2023 for m in 1:12]
ymstr(y, m) = @sprintf("%04d%02d", y, m)
agebin(a) = a < 10 ? "00" : a < 15 ? "10" : a >= 75 ? "75" : string(5 * (a ÷ 5))

# ═════════════════════════════════════════════════════════════════════════════
# 1. Calendar
# ═════════════════════════════════════════════════════════════════════════════
const DTMAP = Dict("평일" => "weekday", "주말" => "weekend", "연휴" => "holiday")
const DOWNAME = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]     # dow: 0=Mon … 6=Sun

# Seollal and Chuseok days (lunar 1/1 · 8/15). A holiday run containing such a day (±1 day) is labelled with its name.
const SEOLLAL = ["20180216", "20190205", "20200125", "20210212", "20220201", "20230122"]
const CHUSEOK = ["20180924", "20190913", "20201001", "20210921", "20220910", "20230929"]

"""
    load_calendar(path)

Read the calendar and attach English day-type and weekday names and a **holiday run label** (`seollal` · `chuseok` · `other`).
A run is a block of consecutive non-working days; a run crossing a month boundary is split between the two months (e.g. Chuseok 2020-09-30 – 10-04).
"""
function load_calendar(path)
    cal = CSV.read(path, DataFrame; types = Dict(:date => String))
    cal.daytype_en = [DTMAP[d] for d in cal.daytype]
    cal.dow_name = [DOWNAME[d+1] for d in cal.dow]
    # run id — a new run starts wherever is_off breaks
    rid = zeros(Int, nrow(cal)); k = 0
    for t in 1:nrow(cal)
        cal.is_off[t] == 1 && (t == 1 || cal.is_off[t-1] == 0) && (k += 1)
        rid[t] = cal.is_off[t] == 1 ? k : 0
    end
    cal.run_id = rid
    near(d, L) = any(abs(Dates.value(Date(d, dateformat"yyyymmdd") - Date(x, dateformat"yyyymmdd"))) <= 1 for x in L)
    lab = Dict{Int,String}()
    for sub in groupby(cal[cal.run_id .> 0, :], :run_id)
        ds = sub.date
        lab[sub.run_id[1]] = any(near(d, SEOLLAL) for d in ds) ? "seollal" :
                             any(near(d, CHUSEOK) for d in ds) ? "chuseok" : "other"
    end
    cal.holiday_kind = [cal.daytype_en[t] == "holiday" ? lab[cal.run_id[t]] : "" for t in 1:nrow(cal)]
    cal
end

"Day-type composition of a month — used to weight the TUS day-of-week axis (weekday · Sat · Sun)."
function month_composition(cal, y, m)
    sub = cal[(cal.year .== y) .& (cal.month .== m), :]
    out = Dict{String,NamedTuple}()
    for dt in ("weekday", "weekend", "holiday")
        d = sub[sub.daytype_en .== dt, :]
        nrow(d) == 0 && continue
        n_sat = count(==("Sat"), d.dow_name); n_sun = count(==("Sun"), d.dow_name)
        kinds = sort(unique(d.holiday_kind))
        out[dt] = (n = nrow(d), n_sat = n_sat, n_sun = n_sun, n_wkh = nrow(d) - n_sat - n_sun,
                   kind = dt == "holiday" ? join(kinds, "+") : "")
    end
    out
end

# ═════════════════════════════════════════════════════════════════════════════
# 2. Registered population
# ═════════════════════════════════════════════════════════════════════════════
"""
    load_registered(path, S, IDX)

Returns `P[(y,m)][g]` (adm2 vectors), single-age sums `p03`·`p1214` for the child rules, national sums
`p7579`·`p80` for splitting the 75+ band, and adm2 names `sggnm`. All summed over sex, keyed by `sgg_cd_stable`.
"""
function load_registered(path, S, IDX)
    kp = CSV.read(path, DataFrame; select = [:year, :month, :sgg_cd_stable, :sgg_nm, :age, :pop],
                  types = Dict(:year => Int, :month => Int, :sgg_cd_stable => String,
                               :sgg_nm => String, :age => Int, :pop => Float64))
    N = length(S)
    P = Dict(ym => Dict(g => zeros(N) for g in AG) for ym in YM_ALL)
    p03 = Dict(ym => zeros(N) for ym in YM_ALL); p1214 = Dict(ym => zeros(N) for ym in YM_ALL)
    p7579 = Dict(ym => 0.0 for ym in YM_ALL); p80 = Dict(ym => 0.0 for ym in YM_ALL)
    sggnm = Dict{String,String}()
    for r in eachrow(kp)
        ym = (r.year, r.month); j = IDX[r.sgg_cd_stable]
        P[ym][agebin(r.age)][j] += r.pop
        r.age <= 3 && (p03[ym][j] += r.pop)
        12 <= r.age <= 14 && (p1214[ym][j] += r.pop)
        75 <= r.age <= 79 && (p7579[ym] += r.pop)
        r.age >= 80 && (p80[ym] += r.pop)
        ym == (2022, 9) && (sggnm[r.sgg_cd_stable] = r.sgg_nm)
    end
    (; P, p03, p1214, p7579, p80, sggnm)
end

# ═════════════════════════════════════════════════════════════════════════════
# 3. Day-type means of the de facto population
# ═════════════════════════════════════════════════════════════════════════════
"De facto day-type means → `RAW[(y,m)][(sgg, g)][(dt, tz)] = raw_mean`."
function load_defacto_daytype(path)
    df = CSV.read(path, DataFrame; types = Dict(:sgg_cd => String, :age_group => String,
                                               :timezn => String, :daytype => String))
    RAW = Dict{Tuple{Int,Int},Dict{Tuple{String,String},Dict{Tuple{String,String},Float64}}}()
    for r in eachrow(df)
        d1 = get!(RAW, (r.year, r.month)) do; Dict{Tuple{String,String},Dict{Tuple{String,String},Float64}}() end
        d2 = get!(d1, (r.sgg_cd, r.age_group)) do; Dict{Tuple{String,String},Float64}() end
        d2[(r.daytype, r.timezn)] = r.raw_mean
    end
    RAW
end

"""
    daytype_marginals(RAW_m, P, S, dts)

Eqs (1b)(1c) — take the correction factor from weekday tz03, apply it to tz15 of all three day types
(additive for age group `00`), and renormalize destinations to the origin total.
"""
function daytype_marginals(RAW_m, P, S, dts)
    N = length(S)
    VT = Dict(dt => Dict(g => zeros(N) for g in AG) for dt in dts)
    nadd = 0
    for g in AG, j in 1:N
        d = RAW_m[(S[j], g)]
        raw3wd = d[("weekday", "03")]
        add = (g == "00") || (raw3wd <= 0)
        add && (nadd += 1)
        c = add ? 0.0 : P[g][j] / raw3wd
        δ = add ? P[g][j] - raw3wd : 0.0
        for dt in dts
            r15 = d[(dt, "15")]
            VT[dt][g][j] = add ? r15 + δ : c * r15
        end
    end
    for dt in dts, g in AG
        VT[dt][g] .*= sum(P[g]) / sum(VT[dt][g])
    end
    (VT, nadd)
end

# ═════════════════════════════════════════════════════════════════════════════
# 4. Census 3-way split — census 2020 shares × that month's registered population
# ═════════════════════════════════════════════════════════════════════════════
const P_MID = 0.996   # 12–14 school participation rate (census)
const S_MID = 0.016   # 12–14 share commuting to school in another district (Personal Travel Survey)

"Monthly census commuting (work or school) file → `CMT[(y,m,sgg,g,category)] = pop`."
function load_census_comm(path)
    cm = CSV.read(path, DataFrame; types = Dict(:year => Int, :month => Int, :sgg_cd => String,
                                               :age_group => String, :category => String, :pop => Float64))
    Dict((r.year, r.month, r.sgg_cd, r.age_group, r.category) => r.pop for r in eachrow(cm))
end

"""
    month_components(CMT, y, m, P, p03, p1214, S)

c_in · nm · c_out of (2d-1). Ages 15+ from the census, under 15 by single-age rules
(0–3 non-move · 4–11 within district · 12–14 via p·s). Rounded to 3 decimals.
"""
function month_components(CMT, y, m, P, p03, p1214, S)
    N = length(S)
    c_in = Dict(g => zeros(N) for g in AG); nm = Dict(g => zeros(N) for g in AG)
    c_out = Dict(g => zeros(N) for g in AG)
    for g in AG, j in 1:N
        s = S[j]; R = P[g][j]
        if g == "00"
            V = p03[j]; X = 0.0; gw = R - V
        elseif g == "10"
            q = p1214[j]; X = q * P_MID * S_MID; V = q * (1 - P_MID); gw = R - X - V
        else
            C = get(CMT, (y, m, s, g, "통근통학"), 0.0); X = get(CMT, (y, m, s, g, "타시군구통근통학"), 0.0)
            V = get(CMT, (y, m, s, g, "비통근통학"), 0.0); gw = C - X
        end
        c_in[g][j] = round(gw, digits = 3); nm[g][j] = round(V, digits = 3); c_out[g][j] = round(X, digits = 3)
    end
    (; c_in, nm, c_out)
end

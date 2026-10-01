# src/inputs.jl — reads the fixed data in input/ into the arrays used by the model.
#   Every OD matrix is A[j,i] — **row = origin j (residence, night-time), column = destination i (activity, daytime)**.
#   The index order is the reverse of $O_{ij}$ in the paper's Methods, so this convention is fixed throughout the code.

const AG = ["00","10","15","20","25","30","35","40","45","50","55","60","65","70","75"]

const AGYEARS = Dict("00"=>0:9,  "10"=>10:14,"15"=>15:19,"20"=>20:24,"25"=>25:29,
                     "30"=>30:34,"35"=>35:39,"40"=>40:44,"45"=>45:49,"50"=>50:54,
                     "55"=>55:59,"60"=>60:64,"65"=>65:69,"70"=>70:74,"75"=>75:120)

const AG3  = Dict(g => (g in ("00","10","15") ? "0-19" : g in ("65","70","75") ? "65+" : "20-64") for g in AG)
const AG3S = ["0-19","20-64","65+"]

adm1(s::AbstractString) = s[1:2]        # first 2 digits of the adm2 region code = adm1 region code

"Dict of empty OD matrices, one per age group."
zeromats(N) = Dict(g => zeros(N, N) for g in AG)

# ── (1) marginals — registered population and de facto population ────────────────────
# input/marginals_202209.csv : sgg_cd, age_group, defacto3_corr, defacto15_corr, registered, raw3, raw15
#   raw3·raw15 = workday-mean de facto population (precomputed). defacto*_corr = values already corrected upstream.
#   Here the correction (see Methods) is **recomputed** and checked against the upstream values.
function load_marginals(path)
    df = CSV.read(path, DataFrame; types = Dict(:sgg_cd=>String, :age_group=>String))
    S   = sort(unique(df.sgg_cd)); N = length(S)
    idx = Dict(s => k for (k, s) in enumerate(S))
    P    = Dict(g => zeros(N) for g in AG)   # P^g_j  registered population
    d15  = Dict(g => zeros(N) for g in AG)   # \tilde d^g_j(15)  (computed here)
    ref15= Dict(g => zeros(N) for g in AG)   # upstream defacto15_corr (for comparison)
    nadd = 0
    for r in eachrow(df)
        j = idx[r.sgg_cd]; g = r.age_group
        P[g][j] = r.registered
        if g == "00" || r.raw3 <= 0            # additive correction
            d15[g][j] = r.raw15 + (r.registered - r.raw3); nadd += 1
        else                                    # ratio correction
            d15[g][j] = r.raw15 * (r.registered / r.raw3)
        end
        ref15[g][j] = r.defacto15_corr
    end
    (; S, N, idx, P, d15, ref15, n_additive = nadd)
end

# ── (2) Read OD matrices ────────────────────────────────────────────────────────────
"long csv → A[j,i] per age group. from/to are the code column names."
function load_od(path, idx, N; from::Symbol, to::Symbol, age::Symbol, val::Symbol, filt = nothing)
    df = CSV.read(path, DataFrame; types = Dict(from=>String, to=>String, age=>String, val=>Float64))
    filt === nothing || filter!(filt, df)
    A = zeromats(N)
    for r in eachrow(df)
        f = getproperty(r, from); t = getproperty(r, to)
        (haskey(idx, f) && haskey(idx, t)) || continue
        A[getproperty(r, age)][idx[f], idx[t]] = getproperty(r, val)
    end
    A
end

# ── (3) radiation shape — G_{ij} normalized per origin ───────────────────────
function load_radiation_shape(path, idx, N)
    df = CSV.read(path, DataFrame; types = Dict(:from_cd_stable=>String, :to_cd_stable=>String, :T_ij=>Float64))
    G = zeros(N, N)
    for r in eachrow(df)
        (haskey(idx, r.from_cd_stable) && haskey(idx, r.to_cd_stable)) || continue
        G[idx[r.from_cd_stable], idx[r.to_cd_stable]] = r.T_ij
    end
    for j in 1:N
        s = sum(@view G[j, :]); s > 0 && (G[j, :] ./= s)
    end
    G
end

# ── (4) census decomposition c_in · nm · c_out ───────────────────────────────────────────
function load_regular_components(path, idx, N)
    df = CSV.read(path, DataFrame; types = Dict(:sgg_cd=>String, :age_group=>String))
    c_in  = Dict(g => zeros(N) for g in AG)
    nm    = Dict(g => zeros(N) for g in AG)
    c_out = Dict(g => zeros(N) for g in AG)
    for r in eachrow(df)
        j = idx[r.sgg_cd]
        c_in[r.age_group][j]  = r.gwannae
        nm[r.age_group][j]    = r.bimove
        c_out[r.age_group][j] = r.tasgg
    end
    (; c_in, nm, c_out)
end

# ── (5) non-commute participation rate π^g ─────────────────────────────────────
# Extracts P_a_star (a constant independent of adm2 region) from input/nc_person_row_factor.csv per age group.
# input/person_factors.csv holds the same values rounded to 4 digits plus sample sizes, for documentation.
# The `kappa_shrunk` column of the same file is a coefficient not used here, so it is not read.
function load_pi(path)
    df = CSV.read(path, DataFrame; types = Dict(:from_cd_stable=>String, :age_group=>String))
    π = Dict{String,Float64}()
    for g in AG
        v = unique(round.(df[df.age_group .== g, :P_a_star], digits = 9))
        length(v) == 1 || error("π^$g varies across adm2 regions: $v")
        π[g] = v[1]
    end
    (; π)
end

# ── (6) activity rate f^g_j ────────────────────────────────────────────────────
# ρ(a, j) = "probability that a person of single age a in adm2 region j commutes (work or school)" (see Methods).
# Two variants are computed.
#
#   census_f    — **adopted**. Uses the census 2020 commuting (work or school) rate κ^g_j = (c_in+c_out)/P^g_j as is.
#                 It measures the same thing as the definition of ρ (work or school), at adm2 × age group resolution, with no extra assumptions.
#   obsolete_f  — **deprecated variant** used in earlier development. Census economic activity rate e_J for 20-69,
#                 e_J(70-74) (rejected upstream) for 70–74, and e_J(70-79) for 75–79.
#                 Not used in computation; kept **only to verify reproduction of the original build**.
#
# κ comes directly from the three components in `regular_components_202209.csv`. The source is census 2020
# "Population by commuting (work or school) type" T20/T00; for ages under 12 the upstream child rule is already applied.
# It is constant within an age group, so single-age weighting leaves it unchanged (f^g_j = κ^g_j).
const RHO_VARIANTS = ["census_f", "obsolete_f"]

"Age group containing single age a."
const AGOF = Dict{Int,String}(a => g for g in AG for a in AGYEARS[g])

"Reads the ρ inputs (single-age population · economic activity rate) from input/. Independent of variant."
function load_rate_inputs(kosis_path, econ_path)
    kp = CSV.read(kosis_path, DataFrame;
        select = [:sgg_cd_stable, :age, :pop],
        types  = Dict(:sgg_cd_stable=>String, :age=>Int, :pop=>Float64))
    pop1 = Dict{Tuple{String,Int},Float64}()
    for r in eachrow(kp)
        pop1[(r.sgg_cd_stable, r.age)] = get(pop1, (r.sgg_cd_stable, r.age), 0.0) + r.pop
    end
    ec = CSV.read(econ_path, DataFrame;
        types = Dict(:adm1_code=>String, :age_group=>String, :econ_active_rate=>Float64,
                     :rate_70_74=>Float64, :rate_70_79=>Float64))
    e     = Dict{Tuple{String,String},Float64}()
    e7074 = Dict{String,Float64}(); e7079 = Dict{String,Float64}()
    for r in eachrow(ec)
        e[(r.adm1_code, r.age_group)] = r.econ_active_rate
        if r.age_group == "70+"
            e7074[r.adm1_code] = r.rate_70_74; e7079[r.adm1_code] = r.rate_70_79
        end
    end
    (; pop1, e, e7074, e7079)
end

"Census commuting (work or school) rate κ^g_j = (c_in + c_out)/P^g_j."
function census_commute_rate(comp, P, N)
    Dict(g => [P[g][j] > 0 ? (comp.c_in[g][j] + comp.c_out[g][j]) / P[g][j] : 0.0 for j in 1:N]
         for g in AG)
end

"""
    activity_rate(variant, ri, κ, S, N)

Builds ρ and f. `ri` = output of `load_rate_inputs`, `κ` = output of `census_commute_rate`.
The returned `ρ(a, j)` takes the **adm2 region index j** (because the adopted variant is at adm2 resolution).
For `obsolete_f` only the adm1 region J containing j is used.
"""
function activity_rate(variant::String, ri, κ, S, N)
    variant in RHO_VARIANTS || error("unknown ρ variant: $variant")
    e, e7074, e7079, pop1 = ri.e, ri.e7074, ri.e7079, ri.pop1
    # Deprecated variant — 0–3 none · 4–18 enrolled · 19+ economic activity rate · 80+ none.
    # ⚠️ ① uses e_J(70-74) for 70–74, which was rejected upstream as an "impossible rebound", and
    #    ② misses school commuting of enrolled students (university) aged 19–29.
    ρ_obsolete(a, J) = a <= 3  ? 0.0 :
                       a <= 18 ? 1.0 :
                       a == 19 ? get(e, (J, "15-19"), 0.0) :
                       a <= 69 ? get(e, (J, "$(5*div(a,5))-$(5*div(a,5)+4)"), 0.0) :
                       a <= 74 ? get(e7074, J, 0.0) :
                       a <= 79 ? get(e7079, J, 0.0) : 0.0
    ρ = variant == "census_f" ? ((a, j) -> a <= 120 ? κ[AGOF[a]][j] : 0.0) :
                                ((a, j) -> ρ_obsolete(a, adm1(S[j])))
    f = Dict(g => zeros(N) for g in AG)
    for g in AG, j in 1:N
        s = S[j]; num = 0.0; den = 0.0
        for a in AGYEARS[g]
            p = get(pop1, (s, a), 0.0); den += p; num += p * ρ(a, j)
        end
        f[g][j] = den > 0 ? num / den : 0.0
    end
    (; f, ρ, variant)
end

# ── (7) adm1 region names (for heatmap axis ticks) ────────────────────────────────────
function load_adm1_names(econ_path)
    ec = CSV.read(econ_path, DataFrame; types = Dict(:adm1_code=>String, :adm1_name=>String))
    Dict(r.adm1_code => r.adm1_name for r in eachrow(ec))
end

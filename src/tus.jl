# src/tus.jl — extracts, directly from the PDF, the two items used here from the Time Use Survey (TUS)
#              (2024, National Data Office) report **Vol. 1-1, Time Use Volume (Tables 1-5)**.
#
#   (A) Table 2-2 "Participation rate by age group", ⑨ **travel block sub-codes** → builds λ·ω_base.
#   (B) the ⑨ **travel total row** of the same table             → builds target T2 (Δq).
#
# To keep all computation in Julia, PDF extraction is done in Julia too — we shell out to
# `pdftotext -layout` for the text and parse it here. The PDF is pinned in `contents/`.
#
# Table layout (verified):
#   Vol. 1-1 `Table 2-2` = PDF p.136–171. Each age band spans 3 pages; the ⑨ travel block is on the third.
#   Even pages = Korean version (code + label + 12 values), odd pages = English version (12 values + label). The two pages show different age bands.
#   12 values = day-of-week average (total·male·female) · weekday (t·m·f) · Saturday (t·m·f) · Sunday (t·m·f).
#   The ⑨ block has the same 22 rows in the same order in both versions:
#     91 910 92 921 922 923 924 93 930 94 940 95 951 952 953 954 96 960 97 970 98 980

const TUS_PDF = "제1-1권 생활시간량편(표1-5).pdf"

# Pages for sub-codes — 12 pages (8 single bands + aggregate bands 19+ · 60+ · 65+ · 15–64).
# Aggregate bands are not used for λ·ω, only for cross-checking against Table 1-2.
const TUS_PAGES = [140, 141, 146, 147, 152, 153, 158, 159, 164, 165, 170, 171]
# Pages for totals — **the 8 single age bands only.** 164·165·170·171 are the overlapping aggregate bands
# 19+ · 60+ · 65+ · 15–64; using them as targets would count the same people twice.
const TUS_TOTAL_PAGES = [140, 141, 146, 147, 152, 153, 158, 159]

const TRAVEL_ROWS = ["91","910","92","921","922","923","924","93","930","94","940",
                     "95","951","952","953","954","96","960","97","970","98","980"]
# Value tokens: 12.3 · 12.3* · 12.3** · - (missing)
const VALTOK = r"^(?:\d+\.\d+\*{0,2}|-)$"

parse_val(t) = t == "-" ? 0.0 : parse(Float64, replace(t, "*" => ""))

"Read the age-band label from the page header (`10~19세` · `80세 이상` etc.)."
function band_of(lines)
    for l in lines[1:min(6, length(lines))]
        m = match(r"(\d+~\d+세|\d+세 이상)", l)
        m === nothing || return m.captures[1]
    end
    nothing
end

pdflines(pdf, pg) = split(read(`pdftotext -f $pg -l $pg -layout $pdf -`, String), '\n')

"Extract the ⑨ travel block (22 rows × 12 values) from one page."
function travel_block(lines)
    hdr = findfirst(l -> occursin("TRAVEL", l) ||
                         (occursin("이동", l) && occursin("", l)) ||
                         occursin(r"^[^\d]*이동\s{3,}\d", l), lines)
    hdr === nothing && return nothing
    out = Vector{Vector{Float64}}()
    for l in lines[hdr+1:end]
        toks = [t for t in split(l) if occursin(VALTOK, t)]
        length(toks) < 12 && continue
        push!(out, parse_val.(toks[1:12]))
        length(out) == length(TRAVEL_ROWS) && break
    end
    length(out) == length(TRAVEL_ROWS) || return nothing
    Dict(TRAVEL_ROWS[k] => out[k] for k in eachindex(out))
end

"""
    read_tus_travel(pdfdir) -> (Dict((band, code) => (weekday, Sat, Sun)), bands)

Read the ⑨ travel **sub-codes** of Table 2-2 (total, in %). Input for building λ·ω_base.
"""
function read_tus_travel(pdfdir)
    pdf = joinpath(pdfdir, TUS_PDF)
    isfile(pdf) || error("TUS PDF not found: $pdf")
    res = Dict{Tuple{String,String},NTuple{3,Float64}}()
    bands = String[]
    for pg in TUS_PAGES
        lines = pdflines(pdf, pg)
        band = band_of(lines)
        band === nothing && continue
        blk = travel_block(lines)
        blk === nothing && error("Table 2-2 p.$pg: failed to parse the travel block")
        # Check parent = child rows (91/910 · 93/930 · 94/940 · 96/960 · 97/970 · 98/980)
        for (par, chi) in (("91","910"),("93","930"),("94","940"),
                           ("96","960"),("97","970"),("98","980"))
            maximum(abs.(blk[par] .- blk[chi])) < 1e-9 ||
                error("Table 2-2 p.$pg ($band): $par ≠ $chi — row alignment is off")
        end
        push!(bands, band)
        for (code, v) in blk
            res[(band, code)] = (v[4], v[7], v[10])   # weekday total · Saturday total · Sunday total
        end
    end
    res, unique(bands)
end

"""
    read_tus_travel_total(pdfdir) -> Dict(band => (weekday, Sat, Sun))

Read the ⑨ **travel total** of Table 2-2 (total, in %).
`1 − (⑨ total travel participation rate)` = **share of people who made no trip that day** = the counterpart of our N̂M;
its change relative to weekdays, Δq, is target T2. Raises an `error` unless all 8 bands are found —
failing silently would corrupt the target.
"""
function read_tus_travel_total(pdfdir)
    pdf = joinpath(pdfdir, TUS_PDF)
    isfile(pdf) || error("TUS PDF not found: $pdf")
    res = Dict{String,NTuple{3,Float64}}()
    for pg in TUS_TOTAL_PAGES
        lines = pdflines(pdf, pg)
        band = band_of(lines)
        band === nothing && error("Table 2-2 p.$pg: could not read the age-band header")
        # ⑨ block header line: Korean version has 12 values after `이동`, English version has `TRAVEL` after 12 values.
        k = findfirst(l -> (occursin("이동", l) || occursin("TRAVEL", l)) &&
                           length([t for t in split(l) if occursin(VALTOK, t)]) >= 12, lines)
        k === nothing && error("Table 2-2 p.$pg ($band): failed to parse the ⑨ travel total row")
        v = parse_val.([t for t in split(lines[k]) if occursin(VALTOK, t)][1:12])
        res[band] = (v[4], v[7], v[10])
    end
    length(res) == length(TUS_TOTAL_PAGES) ||
        error("only $(length(res)) bands found (expected 8)")
    res
end

# ── TUS bands ↔ the 15 age groups ───────────────────────────────────────────
# The 8 single bands used for targets. `00` (0-9) is outside the survey scope, so it is **not assigned to any band**
# — it is excluded from targets (for λ·ω it is proxied by `10~19세` with λ set to 0, run.jl §3).
const TUS_BANDS = ["10~19세","20~29세","30~39세","40~49세","50~59세","60~69세","70~79세","80세 이상"]
# For λ·ω (00 included). `75` (75+) is split into 70~79 / 80+ by registered population.
const BAND = Dict("00"=>["10~19세"], "10"=>["10~19세"], "15"=>["10~19세"],
                  "20"=>["20~29세"], "25"=>["20~29세"], "30"=>["30~39세"], "35"=>["30~39세"],
                  "40"=>["40~49세"], "45"=>["40~49세"], "50"=>["50~59세"], "55"=>["50~59세"],
                  "60"=>["60~69세"], "65"=>["60~69세"], "70"=>["70~79세"],
                  "75"=>["70~79세","80세 이상"])
# For target aggregation (00 excluded).
const G2BAND = Dict(g => BAND[g] for g in keys(BAND) if g != "00")

# ── Non-commute travel codes ─────────────────────────────────────────────────
#   wide   = 910(other work-related) + 940(shopping) + 95(medical·public) + 960(leisure) + 970(socializing) + 980(other)
#   narrow = 940(shopping) + 970(socializing) + 980(culture·leisure) — closest to KTDB trip purposes **(used)**
const NC_WIDE   = ["910","940","95","960","970","980"]
const NC_NARROW = ["940","970","980"]
"Union under independence (% → fraction). Correct by definition, since π^g is a **person** participation rate."
uni(ps) = 1 - prod(1 .- ps ./ 100)
"Simple sum (an upper bound, closer to trip counts)."
sm(ps) = sum(ps) / 100

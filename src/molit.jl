# src/molit.jl — two numbers from the MOLIT press release "2022 Chuseok Special Traffic Measures".
#   Source : 220906(석간)_2022년_추석_특별교통대책(교통정책총괄과).pdf (released 2022-09-05)
#
# ── Basis for ψ ────────────────────────────────────────────────────────────────
# Original text (translated): *"According to the Korea Expressway Corporation, the daily average number of vehicles
# using expressways during the special traffic measures period is expected to be about 5.42 million (up 13.4% from
# the previous year), about 20% higher than **an ordinary weekend (4.5 million)**."*
#
# **Why this pair.** It is a holiday value and a weekend value from the same agency and the same measurement
# (daily average expressway vehicles), so **a ratio can be formed.** Expressway traffic is a public administrative
# statistic, not a paid mobile-phone OD, so it respects this study's premise ("without paid mobile-phone OD").
#
# ⚠️ **ψ is an upper bound on the participation-rate multiplier.** 542/450 is a multiplier of **traffic volume**
#   (vehicle counts), whereas `π^g`, which ω multiplies, is a **person participation rate**. On holidays both the
#   number of people going out and trips per person increase, so the volume multiplier always exceeds the
#   participation multiplier. Hence the bracket [1.000, 1.204] with its geometric mean as the prior median (see Methods).
#
# ── Total travelers 30.17 million — **held-out** ──────────────────────────────────
# The total number of travelers in the same press release (30.17 million · 6.03 million per day) is **not used in the
# likelihood**. It is still computed and reported, but for validation only:
#   ① there is **no weekday counterpart** measured the same way, so no ratio can be formed.
#   ② its universe is **inter-regional long-distance travel** (90.6% by car), unlike our inter-district presence.
#   ③ it is a **pre-holiday survey forecast**, not a post-hoc measurement (Korea Transport Institute 2022 Chuseok holiday travel survey).
#   ④ using it as a target would lose a strong held-out validation — without changing the result.

const MOLIT = (hw_holiday  = 5.42e6,   # daily average expressway vehicles, Chuseok special traffic measures period
               hw_weekend  = 4.50e6,   # daily average expressway vehicles, ordinary weekend
               trips_total = 3.017e7,  # total travelers (5-day cumulative forecast) — held-out
               trips_daily = 6.03e6,   # daily average travelers — held-out
               n_days_policy = 5,
               daily = [5.74e6, 6.09e6, 7.58e6, 6.24e6, 4.52e6],   # 9/8 … 9/12
               dates = ["20220908","20220909","20220910","20220911","20220912"],
               car_share = 0.906)

"**Upper bound** ψ_max of the holiday amplification factor — daily average expressway vehicles holiday / weekend = 542/450 = 1.2044."
const PSI_MAX = MOLIT.hw_holiday / MOLIT.hw_weekend
"Lower bound — *\"people move on holidays exactly as on an ordinary weekend\"*."
const PSI_MIN = 1.0
"prior median — geometric mean of the bracket ends √(1.000 × 1.204) = 1.0975."
const PSI_MID = sqrt(PSI_MIN * PSI_MAX)

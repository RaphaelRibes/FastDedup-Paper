// ── Helpers ────────────────────────────────────────────────────────────────

#let _mean(a) = a.sum() / a.len()

// Sample standard deviation (Bessel-corrected, divides by n − 1)
#let _std(a) = {
  if a.len() <= 1 { 0.0 }
  else {
    let m = _mean(a)
    calc.sqrt(a.map(x => calc.pow(x - m, 2)).sum() / (a.len() - 1))
  }
}

#let _fmt(v, mem: false) = {
  if mem or v >= 100 { str(int(calc.round(v))) }
  else {
    let r = calc.round(v * 10) / 10
    let i = int(calc.floor(r))
    str(i) + "." + str(int(calc.round((r - i) * 10)))
  }
}

#let _fmt-reads(rs) = {
  let v = float(rs)
  let i = int(calc.floor(v))
  str(i) + "." + str(int(calc.round((v - i) * 10))) + "M"
}


// ── Aggregation (now also returns n) ──────────────────────────────────────

#let _agg(data, tool, rs) = {
  let rows = data.filter(r => r.tool == tool and r.reads == rs)
  let sv   = rows.map(r => r.s      ).filter(v => v != "OOM" and v != "").map(float)
  let rv   = rows.map(r => r.max_rss).filter(v => v != "OOM" and v != "").map(float)
  if sv.len() == 0 or rv.len() == 0 {
    (oom: true)
  } else {
    (oom: false,
     sm: _mean(sv), ss: _std(sv), sn: sv.len(),
     rm: _mean(rv), rstd: _std(rv), rn: rv.len())
  }
}


// ── One-sided t-critical value table ──────────────────────────────────────
// Each sub-dictionary maps integer df → critical value.
// For df > 30 the normal approximation is used (last entry = df 999).

#let _t-crit-005 = (
  // α = 0.05 one-sided  (used when k = 1 comparison, i.e. only 2 tools)
  "1": 6.314, "2": 2.920, "3": 2.353, "4": 2.132, "5": 2.015,
  "6": 1.943, "7": 1.895, "8": 1.860, "9": 1.833, "10": 1.812,
  "11": 1.796, "12": 1.782, "13": 1.771, "14": 1.761, "15": 1.753,
  "16": 1.746, "17": 1.740, "18": 1.734, "19": 1.729, "20": 1.725,
  "25": 1.708, "30": 1.697, "999": 1.645,
)

#let _t-crit-0025 = (
  // α = 0.025 one-sided  (used when k = 2)
  "1": 12.706, "2": 4.303, "3": 3.182, "4": 2.776, "5": 2.571,
  "6": 2.447, "7": 2.365, "8": 2.306, "9": 2.262, "10": 2.228,
  "11": 2.201, "12": 2.179, "13": 2.160, "14": 2.145, "15": 2.131,
  "16": 2.120, "17": 2.110, "18": 2.101, "19": 2.093, "20": 2.086,
  "25": 2.060, "30": 2.042, "999": 1.960,
)

#let _t-crit-001 = (
  // α = 0.01 one-sided  (used when k = 3–5)
  "1": 31.821, "2": 6.965, "3": 4.541, "4": 3.747, "5": 3.365,
  "6": 3.143, "7": 2.998, "8": 2.896, "9": 2.821, "10": 2.764,
  "11": 2.718, "12": 2.681, "13": 2.650, "14": 2.624, "15": 2.602,
  "16": 2.583, "17": 2.567, "18": 2.552, "19": 2.539, "20": 2.528,
  "25": 2.485, "30": 2.457, "999": 2.326,
)

#let _t-crit-0005 = (
  // α = 0.005 one-sided  (used when k = 6–10)
  "1": 63.657, "2": 9.925, "3": 5.841, "4": 4.604, "5": 4.032,
  "6": 3.707, "7": 3.499, "8": 3.355, "9": 3.250, "10": 3.169,
  "11": 3.106, "12": 3.055, "13": 3.012, "14": 2.977, "15": 2.947,
  "16": 2.921, "17": 2.898, "18": 2.878, "19": 2.861, "20": 2.845,
  "25": 2.787, "30": 2.750, "999": 2.576,
)


/// Pick the most conservative available α that is ≤ the Bonferroni-corrected
/// per-comparison α, then look up the critical value for the given df.
///
///   alpha-corr : 0.05 / k   (k = number of pairwise comparisons)
///   df         : Welch-Satterthwaite degrees of freedom (integer, ≥ 1)
///
#let _get-t-crit(alpha-corr, df) = {
  // Select table: largest tabled α that is still ≤ alpha-corr  (conservative)
  let tbl = if alpha-corr >= 0.05  { _t-crit-005  }
       else if alpha-corr >= 0.025 { _t-crit-0025 }
       else if alpha-corr >= 0.01  { _t-crit-001  }
       else                        { _t-crit-0005 }

  // Clamp df into the table's key set
  let df-key = if df >= 30 { "999" }
          else if df >= 25 { "25"  }
          else if df >= 20 { "20"  }
          else if df < 1   { "1"   }
          else             { str(df) }

  tbl.at(df-key)
}


// ── Welch's t-test ────────────────────────────────────────────────────────
///
/// One-sided test:  H₀ : μ_a ≤ μ_b   vs   H₁ : μ_a > μ_b
///
/// Returns `true` when we FAIL to reject H₀, i.e. tool A is NOT
/// significantly slower than tool B at the Bonferroni-corrected level.
///
/// Special cases handled:
///   • Identical stats          → not significant  (same tool or exact tie)
///   • Both SDs = 0             → compare means directly
///   • One group has n ≤ 1      → fall back to non-significant (too few data)
///
#let _welch-not-sig(mean-a, sd-a, n-a, mean-b, sd-b, n-b, alpha-corr) = {
  // Trivial: if A is faster or equal, it's certainly not worse
  if mean-a <= mean-b { return true }

  // Need at least 2 observations in each group for a variance estimate
  if n-a <= 1 or n-b <= 1 { return true }

  // Both variances zero ⇒ deterministic; means already compared above
  if sd-a == 0.0 and sd-b == 0.0 { return false }

  // Pooled standard error
  let v-a = calc.pow(sd-a, 2) / n-a
  let v-b = calc.pow(sd-b, 2) / n-b
  let se  = calc.sqrt(v-a + v-b)
  if se == 0.0 { return false }          // shouldn't happen after the check above

  // t statistic  (positive when A is slower than B)
  let t = (mean-a - mean-b) / se

  // Welch-Satterthwaite degrees of freedom  (floored to integer)
  let num = calc.pow(v-a + v-b, 2)
  let den = calc.pow(v-a, 2) / (n-a - 1) + calc.pow(v-b, 2) / (n-b - 1)
  let df  = if den == 0.0 { n-a + n-b - 2 }
            else { int(calc.floor(num / den)) }
  let df  = calc.max(df, 1)

  // Compare against critical value
  let t-crit = _get-t-crit(alpha-corr, df)
  t < t-crit
}


// ── Cell renderer ─────────────────────────────────────────────────────────
///
///  cd     : aggregated stats for this tool × read-count cell
///  cbest  : (time: stats-or-none, ram: stats-or-none) of the best tool
///  tool   : tool name string  (prinseqpp is excluded from comparison)
///  k      : number of Bonferroni comparisons
///
#let _render(cd, cbest, tool, k) = {
  if cd.oom { [OOM] }
  else {
    let t = _fmt(cd.sm)            + "±" + _fmt(cd.ss)
    let m = _fmt(cd.rm, mem: true) + "±" + _fmt(cd.rstd, mem: true)

    let alpha-corr = 0.05 / calc.max(k, 1)

    // ── Time comparison ──
    let is-t-best = (
      tool != "prinseqpp"
      and cbest.time != none
      and _welch-not-sig(
            cd.sm,  cd.ss,  cd.sn,
            cbest.time.sm, cbest.time.ss, cbest.time.sn,
            alpha-corr)
    )
    let tc = if is-t-best { text(fill: green, t) } else { t }

    // ── RAM comparison ──
    let is-m-best = (
      tool != "prinseqpp"
      and cbest.ram != none
      and _welch-not-sig(
            cd.rm,  cd.rstd, cd.rn,
            cbest.ram.rm, cbest.ram.rstd, cbest.ram.rn,
            alpha-corr)
    )
    let mc = if is-m-best { text(fill: green, m) } else { m }

    tc + linebreak() + mc
  }
}


// ── Public entry point ────────────────────────────────────────────────────

#let bench-table(
  csv-path,
  caption:    [Runtime _(s)_ and peak memory _(MB)_ of each tool on synthetic
               libraries of varying size. #text(fill: green)[Green] values are
               not significantly slower/larger than the best tool (one-sided
               Welch's _t_-test, family-wise #sym.alpha = 0.05, Bonferroni
               correction). OOM = Out Of Memory.],
  label:      none,
  tool-order: ("cdhitdup", "clumpify", "fastp", "fastuniq",
               "fdedup_h128", "fdedup_h64", "prinseqpp"),
) = {
  let data       = csv(csv-path, row-type: dictionary)
  let reads-strs = data.map(r => r.reads).dedup().sorted(key: float)
  let raw-tools  = data.map(r => r.tool).dedup()
  let tools      = {
    let ordered = tool-order.filter(t => raw-tools.contains(t))
    let rest    = raw-tools.filter(t => not tool-order.contains(t)).sorted()
    (..ordered, ..rest)
  }

  // Number of comparable tools (excluding prinseqpp)
  let comparable = tools.filter(t => t != "prinseqpp")
  // k = number of pairwise comparisons against the best (all others)
  let k = calc.max(comparable.len() - 1, 1)

  // Aggregate every tool × read-count cell
  let cells = tools.map(t => reads-strs.map(rs => _agg(data, t, rs)))

  // For each read-count column, find the best tool's full stats
  // (separately for time and RAM).  We keep the whole _agg record so
  // the t-test has access to mean, SD, *and* n.
  let col-bests = reads-strs.enumerate().map(pair => {
    let ci = pair.at(0)

    let valid = tools.zip(cells)
      .filter(t_c => t_c.at(0) != "prinseqpp")
      .map(t_c => t_c.at(1).at(ci))
      .filter(cd => not cd.oom)

    // Helper: find the record with the smallest value for field `f`
    let best(f) = if valid.len() == 0 { none }
      else {
        valid.fold(none, (acc, cd) =>
          if acc == none or cd.at(f) < acc.at(f) { cd } else { acc }
        )
      }

    (time: best("sm"), ram: best("rm"))
  })

  // ── Build table ──
  let header = (
    ([*Tool* \ _(s)_ \ _(MB)_],),
    reads-strs.map(rs => [*#_fmt-reads(rs)*]),
  ).flatten()

  let body = tools.enumerate().map(pair => {
    let ti   = pair.at(0)
    let tool = pair.at(1)
    (
      ([#raw(tool)],),
      reads-strs.enumerate().map(ip =>
        _render(cells.at(ti).at(ip.at(0)), col-bests.at(ip.at(0)), tool, k)
      ),
    ).flatten()
  }).flatten()

  [
    #figure(
      align(center)[
        #table(columns: reads-strs.len() + 1, align: center + horizon,
               ..header, ..body)
      ],
      caption: caption,
    )
    #if label != none { label }
  ]
}


// ==========================================================================
//  Duplication-rate table  (unchanged logic, only _std / _mean bug-fixed)
// ==========================================================================

// Redefine helpers in this scope with the same Bessel-corrected _std
#let _mean(a) = if a.len() == 0 { 0 } else { a.sum() / a.len() }
#let _std(a)  = {
  let m = _mean(a)
  if a.len() <= 1 { 0 }
  else { calc.sqrt(a.map(x => calc.pow(x - m, 2)).sum() / (a.len() - 1)) }
}

#let _fmt(v) = {
  let r = calc.round(v * 100) / 100
  let i = int(calc.floor(r))
  let decimals = int(calc.round((r - i) * 100))
  str(i) + "." + (if decimals < 10 { "0" } else { "" }) + str(decimals)
}

#let _fmt-reads(rs) = {
  let v = float(rs)
  if v >= 1 { str(int(v)) + "M" }
  else { str(v) + "M" }
}

#let _agg(data, tool, rs) = {
  let rows  = data.filter(r => r.tool == tool and r.reads == rs)
  let rates = rows.map(r => r.tool_rate).filter(v => v != "").map(float)

  if rates.len() == 0 { (empty: true) }
  else { (empty: false, mean: _mean(rates), std: _std(rates)) }
}


#let duplication-rate-table(
  csv-path,
  caption: [Mean duplication rate (%) across different library sizes.],
  label: none,
  tool-order: ("fdedup_h128", "fdedup_h64", "fastp", "clumpify",
               "cdhitdup", "prinseqpp", "seqkit", "fastuniq"),
) = {
  let data       = csv(csv-path, row-type: dictionary)
  let reads-strs = data.map(r => r.reads).dedup().sorted(key: float)
  let raw-tools  = data.map(r => r.tool).dedup()

  let tools = {
    let ordered = tool-order.filter(t => raw-tools.contains(t))
    let rest    = raw-tools.filter(t => not tool-order.contains(t)).sorted()
    (..ordered, ..rest).filter(t => t != "fastuniq")
  }

  let header = (
    ([*Tool\ Duplication (%)*],),
    reads-strs.map(rs => [*#_fmt-reads(rs)*]),
  ).flatten()

  let body = tools.map(tool => {
    let row-cells = reads-strs.map(rs => {
      let res = _agg(data, tool, rs)
      if res.empty { [---] }
      else { _fmt(res.mean) }
    })
    (([#raw(tool)],), row-cells).flatten()
  }).flatten()

  let fig = figure(
    table(
      columns: reads-strs.len() + 1,
      align: center + horizon,
      stroke: 0.5pt + gray,
      inset: 6pt,
      ..header,
      ..body,
    ),
    caption: caption,
  )

  [#fig #label]
}
// ==========================================================================
//  comparison-tables.typ  –  Welch's t-test edition
//
//  Changes from original:
//    1. Single `build-comparison-table` replaces both `build-comparaison-table`
//       and `build-ram-comparaison-table` via the `metric` parameter.
//    2. The arbitrary 5 % threshold is replaced by a one-sided Welch's t-test
//       with Bonferroni correction (family-wise α = 0.05).  A tool's color
//       appears in a cell iff it is NOT significantly slower/larger than the
//       best tool in that cell.
//    3. Raw replicate values are kept so we can compute sample SD and run
//       the test — the old sum/count accumulator is gone.
//    4. OOM / ERR cells are rendered with a distinct hatched pattern instead
//       of being left blank (which was indistinguishable from "not tested").
//    5. The legend is always shown (was commented-out for the RAM table).
//    6. `_std` uses Bessel's correction (n − 1).
//
//  Backwards-compatible wrappers at the bottom preserve the old call-sites.
// ==========================================================================


// ── Statistics helpers ────────────────────────────────────────────────────

#let _mean(a) = a.sum() / a.len()

#let _std(a) = {
  if a.len() <= 1 { 0.0 }
  else {
    let m = _mean(a)
    calc.sqrt(a.map(x => calc.pow(x - m, 2)).sum() / (a.len() - 1))
  }
}


// ── One-sided t-critical value tables ─────────────────────────────────────
// (identical to the ones in bench-table.typ so both files are self-contained)

#let _t-crit-005 = (
  "1": 6.314, "2": 2.920, "3": 2.353, "4": 2.132, "5": 2.015,
  "6": 1.943, "7": 1.895, "8": 1.860, "9": 1.833, "10": 1.812,
  "11": 1.796, "12": 1.782, "13": 1.771, "14": 1.761, "15": 1.753,
  "16": 1.746, "17": 1.740, "18": 1.734, "19": 1.729, "20": 1.725,
  "25": 1.708, "30": 1.697, "999": 1.645,
)
#let _t-crit-0025 = (
  "1": 12.706, "2": 4.303, "3": 3.182, "4": 2.776, "5": 2.571,
  "6": 2.447, "7": 2.365, "8": 2.306, "9": 2.262, "10": 2.228,
  "11": 2.201, "12": 2.179, "13": 2.160, "14": 2.145, "15": 2.131,
  "16": 2.120, "17": 2.110, "18": 2.101, "19": 2.093, "20": 2.086,
  "25": 2.060, "30": 2.042, "999": 1.960,
)
#let _t-crit-001 = (
  "1": 31.821, "2": 6.965, "3": 4.541, "4": 3.747, "5": 3.365,
  "6": 3.143, "7": 2.998, "8": 2.896, "9": 2.821, "10": 2.764,
  "11": 2.718, "12": 2.681, "13": 2.650, "14": 2.624, "15": 2.602,
  "16": 2.583, "17": 2.567, "18": 2.552, "19": 2.539, "20": 2.528,
  "25": 2.485, "30": 2.457, "999": 2.326,
)
#let _t-crit-0005 = (
  "1": 63.657, "2": 9.925, "3": 5.841, "4": 4.604, "5": 4.032,
  "6": 3.707, "7": 3.499, "8": 3.355, "9": 3.250, "10": 3.169,
  "11": 3.106, "12": 3.055, "13": 3.012, "14": 2.977, "15": 2.947,
  "16": 2.921, "17": 2.898, "18": 2.878, "19": 2.861, "20": 2.845,
  "25": 2.787, "30": 2.750, "999": 2.576,
)

#let _get-t-crit(alpha-corr, df) = {
  let tbl = if alpha-corr >= 0.05  { _t-crit-005  }
       else if alpha-corr >= 0.025 { _t-crit-0025 }
       else if alpha-corr >= 0.01  { _t-crit-001  }
       else                        { _t-crit-0005 }
  let df-key = if df >= 30 { "999" }
          else if df >= 25 { "25"  }
          else if df >= 20 { "20"  }
          else if df <  1  { "1"   }
          else             { str(df) }
  tbl.at(df-key)
}


// ── Welch's one-sided t-test ──────────────────────────────────────────────
/// Returns `true` when we FAIL to reject H₀ : μ_a ≤ μ_b,
/// i.e. tool A is NOT significantly worse than tool B.

#let _welch-not-sig(vals-a, vals-b, alpha-corr) = {
  let n-a = vals-a.len()
  let n-b = vals-b.len()
  if n-a <= 1 or n-b <= 1 { return true }

  let mean-a = _mean(vals-a)
  let mean-b = _mean(vals-b)
  // If A is at least as good, trivially not worse
  if mean-a <= mean-b { return true }

  let sd-a = _std(vals-a)
  let sd-b = _std(vals-b)
  if sd-a == 0.0 and sd-b == 0.0 { return false }   // means differ, no variance

  let v-a = calc.pow(sd-a, 2) / n-a
  let v-b = calc.pow(sd-b, 2) / n-b
  let se  = calc.sqrt(v-a + v-b)
  if se == 0.0 { return false }

  let t = (mean-a - mean-b) / se

  let num = calc.pow(v-a + v-b, 2)
  let den = calc.pow(v-a, 2) / (n-a - 1) + calc.pow(v-b, 2) / (n-b - 1)
  let df  = if den == 0.0 { n-a + n-b - 2 }
            else { int(calc.floor(num / den)) }
  let df  = calc.max(df, 1)

  t < _get-t-crit(alpha-corr, df)
}


// ── Unified comparison table builder ──────────────────────────────────────
///
///  files         : array of CSV file paths
///  metric        : "time" or "ram"
///  exclude-tools : tools to exclude entirely  (e.g. ("prinseqpp",) for RAM)
///  show-legend   : whether to render the color legend
///
#let build-comparison-table(
  files,
  metric:        "time",
  exclude-tools: (),
  show-legend:   true,
) = {

  // ── 1. Parse CSVs, keep raw replicate values per group ──────────────
  let raw-data = ()
  for f in files {
    let lines = csv(f)
    for row in lines.slice(1) {
      let min-cols = if metric == "ram" { 6 } else { 5 }
      if row.len() >= min-cols {
        raw-data.push((
          reads:  row.at(0),
          mode:   row.at(1),
          in_out: row.at(2),
          tool:   row.at(3),
          time:   row.at(4),
          value:  if metric == "ram" { row.at(5) } else { row.at(4) },
        ))
      }
    }
  }

  // grouped:  key → (values: array-of-float)
  // oom-keys: set of keys where at least one OOM / ERR was seen
  let grouped  = (:)
  let oom-keys = (:)

  for entry in raw-data {
    if exclude-tools.contains(entry.tool) { continue }
    let norm-reads = str(float(entry.reads))
    let key = norm-reads + "|" + entry.in_out + "|" + entry.mode + "|" + entry.tool
    // Track OOM / ERR
    if (entry.time == "OOM" or entry.time == "ERR"
       or entry.value == "OOM" or entry.value == "ERR") {
      oom-keys.insert(key, true)
      continue
    }
    if key not in grouped { grouped.insert(key, ()) }
    grouped.insert(key, (..grouped.at(key), float(entry.value)))
  }

  // ── 2. Discover unique read sizes and tools ─────────────────────────
  let unique-reads-str = ()
  let unique-tools     = ()
  for entry in raw-data {
    if exclude-tools.contains(entry.tool) { continue }
    let norm-reads = str(float(entry.reads))
    if norm-reads not in unique-reads-str { unique-reads-str.push(norm-reads) }
    if entry.tool not in unique-tools     { unique-tools.push(entry.tool)     }
  }
  let unique-reads = unique-reads-str.map(r => float(r)).sorted()
  unique-reads-str = unique-reads.map(r => str(r))

  let modes   = ("paired", "single")
  let in_outs = ("cc", "dd", "cd")

  // Number of Bonferroni comparisons  (all tools minus the best one)
  let k-bonf = calc.max(unique-tools.len() - 1, 1)
  let alpha-corr = 0.05 / k-bonf

  // ── 3. Color palette ────────────────────────────────────────────────
  let oom-color = rgb("E0E0E0") // distinct grey for OOM cells

  // Statically map tools to colors so they don't shift when tools are excluded
  let base-colors = (
    "fastp":     rgb("9DBCFF"),
    "clumpify":  rgb("9DFFB0"),
    "prinseqpp": rgb("FFB89D"),
    "cdhitdup":  rgb("9DFFF5"),
    "fastuniq":  rgb("FFE49D"),
    "seqkit":    rgb("D0B0FF"),
  )

  let fallback-palette = (rgb("C9FF9D"), rgb("FFB89D"), rgb("FFE49D"))

  let tool-colors = (:)
  let ci = 0

  for t in unique-tools {
    if t.starts-with("fdedup") {
      tool-colors.insert(t, rgb("FF9D9D"))
    }
    else if t in base-colors {
      tool-colors.insert(t, base-colors.at(t))
    }
    else {
      // Fallback if an unexpected tool appears in the CSV
      tool-colors.insert(t, fallback-palette.at(calc.rem(ci, fallback-palette.len())))
      ci += 1
    }
  }

  // ── 4. Legend ───────────────────────────────────────────────────────
  let legend-entries = ()
  let seen-colors    = ()
  for t in unique-tools {
    let c = tool-colors.at(t)
    if c not in seen-colors {
      seen-colors.push(c)
      let label = if t.starts-with("fdedup") { "fdedup" } else { t }
      legend-entries.push((color: c, label: label))
    }
  }

  let legend = align(center)[
    #legend-entries.map(entry => {
      box(width: 1em, height: 1em, fill: entry.color, radius: 2pt)
      h(0.5em)
      entry.label
    }).join(h(2em))
    #v(1em)
  ]

  // ── 5. Table header ─────────────────────────────────────────────────
  let total-cols = 1 + unique-reads-str.len() * 3
  let header = (table.cell(rowspan: 2, [*Mode*]),)
  for r in unique-reads-str {
    header.push(table.cell(colspan: 3, [*#r M*]))
  }
  let subheader = ()
  for r in unique-reads-str {
    subheader.push([*cc*])
    subheader.push([*dd*])
    subheader.push([*cd*])
  }

  // ── 6. Build table body ─────────────────────────────────────────────
  let table-rows = ()
  for m in modes {
    let print-m = if m == "paired" { "PE" } else { "SE" }
    table-rows.push([*#print-m*])

    for r in unique-reads-str {
      for io in in_outs {
        // Collect every tool that has data for this cell
        let cell-tools = ()         // ( (tool, values), … )
        let any-oom    = false
        for t in unique-tools {
          let key = r + "|" + io + "|" + m + "|" + t
          if key in grouped and grouped.at(key).len() > 0 {
            cell-tools.push((tool: t, vals: grouped.at(key)))
          }
          if key in oom-keys { any-oom = true }
        }

        if cell-tools.len() == 0 {
          // No tool succeeded → mark OOM if any tried, else empty
          if any-oom {
            table-rows.push(table.cell(fill: oom-color, []))
          } else {
            table-rows.push([])
          }
        } else {
          // ── Find the best tool (lowest mean) ──
          let best-idx = 0
          let best-mean = _mean(cell-tools.at(0).vals)
          for idx in range(1, cell-tools.len()) {
            let m-cur = _mean(cell-tools.at(idx).vals)
            if m-cur < best-mean {
              best-mean = m-cur
              best-idx  = idx
            }
          }
          let best-vals = cell-tools.at(best-idx).vals

          // ── Test every tool against the best ──
          let best-colors = ()
          for entry in cell-tools {
            let is-not-sig = _welch-not-sig(entry.vals, best-vals, alpha-corr)
            if is-not-sig {
              let c = tool-colors.at(entry.tool)
              if c not in best-colors { best-colors.push(c) }
            }
          }

          // Render the cell
          if best-colors.len() == 0 {
            // Shouldn't happen (the best tool always passes), but safety net
            table-rows.push([])
          } else if best-colors.len() == 1 {
            table-rows.push(table.cell(fill: best-colors.at(0), []))
          } else {
            table-rows.push(
              table.cell(fill: gradient.linear(dir: ltr, ..best-colors), [])
            )
          }
        }
      }
    }
  }

  // ── 7. Assemble ─────────────────────────────────────────────────────
  [
    #if show-legend { legend }
    #table(
      columns: total-cols,
      align: center + horizon,
      ..header,
      ..subheader,
      ..table-rows,
    )
  ]
}


// ── Backwards-compatible wrappers ─────────────────────────────────────────
// Drop-in replacements so existing call-sites keep compiling unchanged.

#let build-comparaison-table(files) = {
  build-comparison-table(files, metric: "time")
}

#let build-ram-comparaison-table(files) = {
  build-comparison-table(
    files,
    metric:        "ram",
    exclude-tools: ("prinseqpp",),
    show-legend:   true,
  )
}

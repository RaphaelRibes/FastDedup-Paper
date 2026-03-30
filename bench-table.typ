// ╔══════════════════════════════════════════════════════════════════════╗
// ║  bench-table.typ  –  reusable benchmark table module                 ║
// ║                                                                      ║
// ║  Usage in your main file:                                            ║
// ║    #import "bench-table.typ": bench-table                            ║
// ║    #bench-table("dd.csv", caption: [...], label: <PE_dd_table>)      ║
// ║                                                                      ║
// ║  Required CSV columns: reads, tool, s, max_rss                       ║
// ║  (any extra columns are ignored)                                     ║
// ╚══════════════════════════════════════════════════════════════════════╝


// ── Internal helpers (not exported) ─────────────────────────────────────

#let _mean(a)    = a.sum() / a.len()
#let _std(a)     = { let m = _mean(a); calc.sqrt(a.map(x => calc.pow(x - m, 2)).sum() / a.len()) }

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

#let _agg(data, tool, rs) = {
  let rows  = data.filter(r => r.tool == tool and r.reads == rs)
  let sv    = rows.map(r => r.s      ).filter(v => v != "OOM" and v != "").map(float)
  let rv    = rows.map(r => r.max_rss).filter(v => v != "OOM" and v != "").map(float)
  if sv.len() == 0 or rv.len() == 0 { (oom: true) }
  else { (oom: false, sm: _mean(sv), ss: _std(sv), rm: _mean(rv), rstd: _std(rv)) }
}

#let _render(cd, cmin) = {
  if cd.oom { [OOM] }
  else {
    let t = _fmt(cd.sm)            + "±" + _fmt(cd.ss)
    let m = _fmt(cd.rm, mem: true) + "±" + _fmt(cd.rstd, mem: true)
    let tc = if cmin.sm != none and _fmt(cd.sm) == _fmt(cmin.sm) { text(fill: green, t) } else { t }
    let mc = if cmin.rm != none and _fmt(cd.rm, mem: true) == _fmt(cmin.rm, mem: true) { text(fill: green, m) } else { m }
    tc + linebreak() + mc
  }
}


// ── Public API ───────────────────────────────────────────────────────────

/// Render a benchmark figure from a CSV file.
///
/// Parameters:
///   csv-path   (str)      – path to the CSV file
///   caption    (content)  – figure caption
///   label      (label)    – optional figure label, e.g. <PE_dd_table>
///   tool-order (array)    – preferred row order; unknown tools are appended
#let bench-table(
  csv-path,
  caption:    [Runtime and memory usage benchmarking of existing tools and
               `FastDedup` on synthetic datasets of varying sizes.
               OOM indicates Out Of Memory errors.],
  label:      none,
  tool-order: ("cdhitdup", "clumpify", "fastp", "fastuniq",
               "fdedup_h128", "fdedup_h64", "prinseqpp"),
) = {
  let data        = csv(csv-path, row-type: dictionary)
  let reads-strs  = data.map(r => r.reads).dedup().sorted(key: float)
  let raw-tools   = data.map(r => r.tool).dedup()
  let tools       = {
    let ordered = tool-order.filter(t => raw-tools.contains(t))
    let rest    = raw-tools.filter(t => not tool-order.contains(t)).sorted()
    (..ordered, ..rest)
  }
  
  let cells    = tools.map(t => reads-strs.map(rs => _agg(data, t, rs)))
  let col-mins = reads-strs.enumerate().map(pair => {
    let ci  = pair.at(0)
    let col = cells.map(tc => tc.at(ci)).filter(cd => not cd.oom)
    let mn(f) = if col.len() == 0 { none }
      else { col.map(cd => cd.at(f)).fold(none, (acc, v) => if acc == none or v < acc { v } else { acc }) }
    (sm: mn("sm"), rm: mn("rm"))
  })

  let header = (([*Tool* \ _(s) \ (MB)_],), reads-strs.map(rs => [*#_fmt-reads(rs)*])).flatten()
  let body    = tools.enumerate().map(pair => {
    let ti = pair.at(0); let tool = pair.at(1)
    (([#raw(tool)],), reads-strs.enumerate().map(ip =>
      _render(cells.at(ti).at(ip.at(0)), col-mins.at(ip.at(0)))
    )).flatten()
  }).flatten()

  [
    #figure(
      align(center)[
        #table(columns: reads-strs.len() + 1, align: center + horizon, ..header, ..body)
      ],
      caption: caption,
    )
    #if label != none { label }
  ]
}

// ── Internal helpers ────────────────────────────────────────────────────

#let _mean(a) = if a.len() == 0 { 0 } else { a.sum() / a.len() }
#let _std(a)  = { 
  let m = _mean(a)
  if a.len() <= 1 { 0 } 
  else { calc.sqrt(a.map(x => calc.pow(x - m, 2)).sum() / a.len()) }
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
  let rows = data.filter(r => r.tool == tool and r.reads == rs)
  let rates = rows.map(r => r.tool_rate).filter(v => v != "").map(float)
  
  if rates.len() == 0 { (empty: true) }
  else { (empty: false, mean: _mean(rates), std: _std(rates)) }
}

// ── Public API ──────────────────────────────────────────────────────────

#let duplication-rate-table(
  csv-path,
  caption: [Mean duplication rate (%) across different library sizes.],
  label: none,
  tool-order: ("fdedup_h128", "fdedup_h64", "fastp", "clumpify", "cdhitdup", "prinseqpp", "seqkit", "fastuniq"),
) = {
  let data = csv(csv-path, row-type: dictionary)
  let reads-strs = data.map(r => r.reads).dedup().sorted(key: float)
  let raw-tools  = data.map(r => r.tool).dedup()
  
  let tools = {
    let ordered = tool-order.filter(t => raw-tools.contains(t))
    let rest    = raw-tools.filter(t => not tool-order.contains(t)).sorted()
    (..ordered, ..rest).filter(t => t != "fastuniq")
  }

  let header = (([*Tool\ Duplication (%)*],), reads-strs.map(rs => [*#_fmt-reads(rs)*])).flatten()
  
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
      ..body
    ),
    caption: caption,
  )

  // Use a content block to group the figure and the label together
  [#fig #label] 
}
#let build-comparaison-table(files) = {
  let raw-data = ()
  for f in files {
    let lines = csv(f)
    for row in lines.slice(1) {
      if row.len() >= 5 {
        raw-data.push((
          reads: row.at(0), mode: row.at(1), in_out: row.at(2),
          tool: row.at(3), time: row.at(4)
        ))
      }
    }
  }

  let grouped = (:)
  for entry in raw-data {
    if entry.time == "OOM" or entry.time == "ERR" { continue }
    let norm-reads = str(float(entry.reads))
    let key = norm-reads + "|" + entry.in_out + "|" + entry.mode + "|" + entry.tool
    if key not in grouped { grouped.insert(key, (sum: 0.0, count: 0)) }
    let stats = grouped.at(key)
    stats.sum += float(entry.time)
    stats.count += 1
    grouped.insert(key, stats)
  }

  let averages = (:)
  for (k, v) in grouped { averages.insert(k, v.sum / v.count) }

  let unique-reads-str = ()
  let unique-tools = ()
  for entry in raw-data {
    let norm-reads = str(float(entry.reads))
    if norm-reads not in unique-reads-str { unique-reads-str.push(norm-reads) }
    if entry.tool not in unique-tools { unique-tools.push(entry.tool) }
  }

  let unique-reads = unique-reads-str.map(r => float(r)).sorted()
  unique-reads-str = unique-reads.map(r => str(r))
  
  let modes = ("paired", "single")
  let in_outs = ("cc", "dd", "cd")
  
  let palette = (
    rgb("9DBCFF"), rgb("9DFFB0"), rgb("FFE49D"), 
    rgb("9DFFF5"), rgb("FFB89D"), rgb("C9FF9D")
  )
  let fdedup-color = rgb("FF9D9D")
  let seqkit-color = rgb("D0B0FF")
  
  let tool-colors = (:)
  let i = 0
  for t in unique-tools {
    if t.starts-with("fdedup") { tool-colors.insert(t, fdedup-color) } 
    else if t == "seqkit" { tool-colors.insert(t, seqkit-color) } 
    else { tool-colors.insert(t, palette.at(calc.rem(i, palette.len()))); i += 1 }
  }

  let legend-entries = ()
  let seen-colors = ()
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

  let total-cols = 1 + unique-reads-str.len() * 3
  let header = (table.cell(rowspan: 2, [*Mode*]),)
  for r in unique-reads-str { header.push(table.cell(colspan: 3, [*#r M*])) }
  
  let subheader = ()
  for r in unique-reads-str {
    subheader.push([*cc*]); subheader.push([*dd*]); subheader.push([*cd*])
  }

  let table-rows = ()
  for m in modes {
    let print-m = if m == "paired" { "PE" } else { "SE" }
    table-rows.push([*#print-m*])
    for r in unique-reads-str {
      for io in in_outs {
        // Étape 1 : Trouver le temps absolu le plus rapide
        let min-time = 99999999.0
        for t in unique-tools {
          let k = r + "|" + io + "|" + m + "|" + t
          if k in averages and averages.at(k) < min-time {
            min-time = averages.at(k)
          }
        }
        
        // Étape 2 : Récupérer les couleurs de tous les outils à 10% de ce min
        if min-time != 99999999.0 {
          let threshold = min-time * 1.05
          let best-colors = ()
          
          for t in unique-tools {
            let k = r + "|" + io + "|" + m + "|" + t
            if k in averages and averages.at(k) <= threshold {
              let c = tool-colors.at(t)
              if c not in best-colors {
                best-colors.push(c)
              }
            }
          }
          
          // Dessin de la case
          if best-colors.len() == 1 {
            table-rows.push(table.cell(fill: best-colors.at(0), []))
          } else {
            table-rows.push(table.cell(fill: gradient.linear(dir: ltr, ..best-colors), []))
          }
        } else {
          table-rows.push([])
        }
      }
    }
  }

  [
    #legend
    #table(
      columns: total-cols,
      align: center + horizon,
      ..header,
      ..subheader,
      ..table-rows
    )
  ]
}

// =========================================================================

#let build-ram-comparaison-table(files) = {
  let raw-data = ()
  for f in files {
    let lines = csv(f)
    for row in lines.slice(1) {
      if row.len() >= 6 {
        raw-data.push((
          reads: row.at(0), mode: row.at(1), in_out: row.at(2),
          tool: row.at(3), time: row.at(4), ram: row.at(5)
        ))
      }
    }
  }

  let grouped = (:)
  for entry in raw-data {
    if entry.time == "OOM" or entry.time == "ERR" { continue }
    if entry.tool == "prinseqpp" { continue } // Exclusion de prinseq++
    let norm-reads = str(float(entry.reads))
    let key = norm-reads + "|" + entry.in_out + "|" + entry.mode + "|" + entry.tool
    if key not in grouped { grouped.insert(key, (sum: 0.0, count: 0)) }
    let stats = grouped.at(key)
    stats.sum += float(entry.ram)
    stats.count += 1
    grouped.insert(key, stats)
  }

  let averages = (:)
  for (k, v) in grouped { averages.insert(k, v.sum / v.count) }

  let unique-reads-str = ()
  let unique-tools = ()
  for entry in raw-data {
    if entry.tool == "prinseqpp" { continue }
    let norm-reads = str(float(entry.reads))
    if norm-reads not in unique-reads-str { unique-reads-str.push(norm-reads) }
    if entry.tool not in unique-tools { unique-tools.push(entry.tool) }
  }

  let unique-reads = unique-reads-str.map(r => float(r)).sorted()
  unique-reads-str = unique-reads.map(r => str(r))
  
  let modes = ("paired", "single")
  let in_outs = ("cc", "dd", "cd")
  
  let palette = (
    rgb("9DBCFF"), rgb("9DFFB0"), rgb("FFE49D"), 
    rgb("9DFFF5"), rgb("FFB89D"), rgb("C9FF9D")
  )
  let fdedup-color = rgb("FF9D9D")
  let seqkit-color = rgb("D0B0FF")
  
  let tool-colors = (:)
  let i = 0
  for t in unique-tools {
    if t.starts-with("fdedup") { tool-colors.insert(t, fdedup-color) } 
    else if t == "seqkit" { tool-colors.insert(t, seqkit-color) } 
    else { tool-colors.insert(t, palette.at(calc.rem(i, palette.len()))); i += 1 }
  }

  let total-cols = 1 + unique-reads-str.len() * 3
  let header = (table.cell(rowspan: 2, [*Mode*]),)
  for r in unique-reads-str { header.push(table.cell(colspan: 3, [*#r M*])) }
  
  let subheader = ()
  for r in unique-reads-str {
    subheader.push([*cc*]); subheader.push([*dd*]); subheader.push([*cd*])
  }

  let table-rows = ()
  for m in modes {
    let print-m = if m == "paired" { "PE" } else { "SE" }
    table-rows.push([*#print-m*])
    for r in unique-reads-str {
      for io in in_outs {
        // Étape 1 : Trouver la RAM absolue la plus basse
        let min-ram = 9999999999.0 
        for t in unique-tools {
          let k = r + "|" + io + "|" + m + "|" + t
          if k in averages and averages.at(k) < min-ram {
            min-ram = averages.at(k)
          }
        }
        
        // Étape 2 : Récupérer les couleurs de tous les outils à 10% de ce min
        if min-ram != 9999999999.0 {
          let threshold = min-ram * 1.05
          let best-colors = ()
          
          for t in unique-tools {
            let k = r + "|" + io + "|" + m + "|" + t
            if k in averages and averages.at(k) <= threshold {
              let c = tool-colors.at(t)
              if c not in best-colors {
                best-colors.push(c)
              }
            }
          }
          
          // Dessin de la case
          if best-colors.len() == 1 {
            table-rows.push(table.cell(fill: best-colors.at(0), []))
          } else {
            table-rows.push(table.cell(fill: gradient.linear(dir: ltr, ..best-colors), []))
          }
        } else {
          table-rows.push([])
        }
      }
    }
  }

  [
    #table(
      columns: total-cols,
      align: center + horizon,
      ..header,
      ..subheader,
      ..table-rows
    )
    /* Légende pour le tableau RAM
    #align(center)[
      #v(1em)
      #unique-tools.map(t => {
        let c = tool-colors.at(t)
        box(width: 1em, height: 1em, fill: c, radius: 2pt)
        h(0.5em)
        if t.starts-with("fdedup") { "fdedup" } else { t }
      }).join(h(2em))
    ] */
  ]
}
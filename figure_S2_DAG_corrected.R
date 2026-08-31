#!/usr/bin/env Rscript

# Corrected conceptual DAG for Figure S2.
# This script generates artwork only; it does not load data or fit any model.

required_packages <- c("DiagrammeR", "DiagrammeRsvg", "rsvg")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

output_svg <- "Figure_S2_DAG_corrected.svg"
output_png <- "Figure_S2_DAG_corrected.png"
output_pdf <- "Figure_S2_DAG_corrected.pdf"

# The nine causal nodes and sixteen causal edges are listed explicitly below.
# The title and bottom note are graph annotations, not causal nodes.
causal_nodes <- c(
  "period", "origin", "social", "infection", "prenatal",
  "diagnosis", "congenital", "ascertainment", "recorded"
)

causal_edges <- data.frame(
  from = c(
    "period", "period", "period", "period", "period", "period",
    "origin", "social", "social", "social", "infection", "prenatal",
    "infection", "diagnosis", "congenital", "ascertainment"
  ),
  to = c(
    "origin", "social", "infection", "prenatal", "diagnosis", "ascertainment",
    "social", "infection", "prenatal", "ascertainment", "diagnosis", "diagnosis",
    "congenital", "congenital", "ascertainment", "recorded"
  ),
  stringsAsFactors = FALSE
)

stopifnot(
  length(causal_nodes) == 9L,
  length(unique(causal_nodes)) == 9L,
  nrow(causal_edges) == 16L,
  all(causal_edges$from %in% causal_nodes),
  all(causal_edges$to %in% causal_nodes),
  !any(causal_edges$from == "origin" & causal_edges$to == "infection"),
  !any(causal_edges$from == "prenatal" & causal_edges$to == "infection"),
  !any(causal_edges$from == "origin" & causal_edges$to == "recorded"),
  any(causal_edges$from == "prenatal" & causal_edges$to == "diagnosis")
)

dot <- paste0(
  "digraph corrected_dag {\n",
  "  graph [layout=dot, rankdir=LR, bgcolor=white, pad=0.25, nodesep=0.45, ",
  "ranksep=0.75, splines=polyline, outputorder=edgesfirst, fontname=Arial, ",
  "labelloc=t, labeljust=l, fontsize=24, ",
  "label=\"Conceptual DAG of the Maternal-Origin Contrast\"];\n",
  "  node [shape=box, style=\"rounded,filled\", color=\"#333333\", ",
  "penwidth=1.2, fontname=Arial, fontsize=12, margin=\"0.16,0.10\"];\n",
  "  edge [color=\"#555555\", penwidth=1.25, arrowsize=0.72];\n\n",

  "  period [label=\"Calendar period /\\nmigration / surveillance\", fillcolor=\"#E5E5E5\"];\n",
  "  origin [label=\"Maternal origin\", fillcolor=\"#D9EAF7\"];\n",
  "  social [label=\"Unmeasured social determinants\\nand healthcare access barriers\", fillcolor=\"#EFE3F7\"];\n",
  "  infection [label=\"Maternal syphilis infection\", fillcolor=\"#EFE3F7\"];\n",
  "  prenatal [label=\"Prenatal care utilization\", fillcolor=\"#FFF2CC\"];\n",
  "  diagnosis [label=\"Syphilis diagnosis and\\ntreatment during pregnancy\", fillcolor=\"#FFF2CC\"];\n",
  "  congenital [label=\"Congenital syphilis\", fillcolor=\"#F7D9D5\"];\n",
  "  ascertainment [label=\"Neonatal hospitalization /\\ncase ascertainment\", fillcolor=\"#E5E5E5\"];\n",
  "  recorded [label=\"Hospital-recorded\\ncongenital syphilis\", fillcolor=\"#F7D9D5\"];\n\n",

  # Rank groups organize a left-to-right flow. They do not add causal edges.
  "  { rank=same; period; origin; }\n",
  "  { rank=same; infection; prenatal; }\n",
  "  { rank=same; diagnosis; congenital; }\n\n",

  "  period -> origin;\n",
  "  period -> social;\n",
  "  period -> infection;\n",
  "  period -> prenatal;\n",
  "  period -> diagnosis;\n",
  "  period -> ascertainment;\n",
  "  origin -> social;\n",
  "  social -> infection;\n",
  "  social -> prenatal;\n",
  "  social -> ascertainment;\n",
  "  infection -> diagnosis;\n",
  "  prenatal -> diagnosis;\n",
  "  infection -> congenital;\n",
  "  diagnosis -> congenital;\n",
  "  congenital -> ascertainment;\n",
  "  ascertainment -> recorded;\n",
  "}\n"
)

render_fallback_svg <- function(filename) {
  # Vector-only R renderer used when the Debian DiagrammeR build lacks Viz.js.
  # Coordinates are normalized and correspond to the same validated DAG above.
  grDevices::svg(filename, width = 18, height = 6.8, bg = "white", pointsize = 12)
  grid::grid.newpage()

  xy <- list(
    period = c(.42, .86), origin = c(.065, .50), social = c(.215, .50),
    infection = c(.385, .64), prenatal = c(.385, .38), diagnosis = c(.54, .50),
    congenital = c(.68, .50), ascertainment = c(.81, .50), recorded = c(.945, .50)
  )
  wh <- list(
    period = c(.20, .095), origin = c(.10, .085), social = c(.18, .12),
    infection = c(.14, .085), prenatal = c(.14, .085), diagnosis = c(.16, .12),
    congenital = c(.10, .085), ascertainment = c(.14, .12), recorded = c(.11, .105)
  )
  labels <- list(
    period = "Calendar period /\nmigration / surveillance",
    origin = "Maternal origin",
    social = "Unmeasured social determinants\nand healthcare access barriers",
    infection = "Maternal syphilis infection",
    prenatal = "Prenatal care utilization",
    diagnosis = "Syphilis diagnosis and\ntreatment during pregnancy",
    congenital = "Congenital syphilis",
    ascertainment = "Neonatal hospitalization /\ncase ascertainment",
    recorded = "Hospital-recorded\ncongenital syphilis"
  )
  fills <- c(period = "#E5E5E5", origin = "#D9EAF7", social = "#EFE3F7",
             infection = "#EFE3F7", prenatal = "#FFF2CC", diagnosis = "#FFF2CC",
             congenital = "#F7D9D5", ascertainment = "#E5E5E5", recorded = "#F7D9D5")

  grid::grid.text("Conceptual DAG of the Maternal-Origin Contrast",
                  x = .03, y = .975, just = c("left", "top"),
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 24, fontface = "bold"))

  edge_gp <- grid::gpar(col = "#555555", lwd = 1.8, lineend = "round", linejoin = "round")
  draw_path <- function(x, y) {
    grid::grid.lines(x, y, default.units = "npc", gp = edge_gp,
                     arrow = grid::arrow(type = "closed", length = grid::unit(0.12, "inches")))
  }
  left <- function(n) xy[[n]][1] - wh[[n]][1] / 2
  right <- function(n) xy[[n]][1] + wh[[n]][1] / 2
  top <- function(n) xy[[n]][2] + wh[[n]][2] / 2
  bottom <- function(n) xy[[n]][2] - wh[[n]][2] / 2

  # Period arrows use separate upper lanes. The prenatal route enters from the
  # left so it cannot be mistaken for an infection-to-prenatal-care arrow.
  draw_path(c(.345, .345, xy$origin[1], xy$origin[1]),
            c(bottom("period"), .775, .775, top("origin")))
  draw_path(c(.375, .375, xy$social[1], xy$social[1]),
            c(bottom("period"), .755, .755, top("social")))
  draw_path(c(.405, .405, xy$infection[1], xy$infection[1]),
            c(bottom("period"), .735, .735, top("infection")))
  draw_path(c(.435, .435, .305, .305, left("prenatal")),
            c(bottom("period"), .715, .715, xy$prenatal[2], xy$prenatal[2]))
  draw_path(c(.465, .465, xy$diagnosis[1], xy$diagnosis[1]),
            c(bottom("period"), .695, .695, top("diagnosis")))
  draw_path(c(.495, .495, xy$ascertainment[1], xy$ascertainment[1]),
            c(bottom("period"), .675, .675, top("ascertainment")))

  draw_path(c(right("origin"), left("social")), c(xy$origin[2], xy$social[2]))
  draw_path(c(right("social"), left("infection")), c(.53, .64))
  draw_path(c(right("social"), left("prenatal")), c(.47, .38))
  draw_path(c(xy$social[1], xy$social[1], xy$ascertainment[1], xy$ascertainment[1]),
            c(bottom("social"), .19, .19, bottom("ascertainment")))
  draw_path(c(right("infection"), left("diagnosis")), c(.62, .53))
  draw_path(c(right("prenatal"), left("diagnosis")), c(.40, .47))
  draw_path(c(right("infection"), .47, .47, xy$congenital[1], xy$congenital[1]),
            c(.64, .64, .72, .72, top("congenital")))
  draw_path(c(right("diagnosis"), left("congenital")), c(.50, .50))
  draw_path(c(right("congenital"), left("ascertainment")), c(.50, .50))
  draw_path(c(right("ascertainment"), left("recorded")), c(.50, .50))

  for (n in causal_nodes) {
    grid::grid.roundrect(x = xy[[n]][1], y = xy[[n]][2],
                         width = wh[[n]][1], height = wh[[n]][2], r = grid::unit(.08, "snpc"),
                         gp = grid::gpar(fill = fills[n], col = "#333333", lwd = 1.4))
    grid::grid.text(labels[[n]], x = xy[[n]][1], y = xy[[n]][2],
                    gp = grid::gpar(fontfamily = "Arial", fontsize = 11.5, col = "#222222"))
  }
  grDevices::dev.off()
  invisible(filename)
}

dag <- DiagrammeR::grViz(dot, engine = "dot")
viz_js <- system.file("htmlwidgets/lib/viz/viz.js", package = "DiagrammeR")
if (nzchar(viz_js) && file.exists(viz_js)) {
  svg_text <- tryCatch(
    DiagrammeRsvg::export_svg(dag),
    error = function(e) {
      message(
        "DiagrammeRsvg could not render the SVG (", conditionMessage(e),
        "). Falling back to the built-in vector R renderer."
      )
      svg_file <- tempfile(fileext = ".svg")
      on.exit(unlink(svg_file), add = TRUE)
      render_fallback_svg(svg_file)
      paste(readLines(svg_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    }
  )
} else {
  message("The Debian DiagrammeR package lacks bundled Viz.js; using the built-in vector R renderer.")
    svg_file <- tempfile(fileext = ".svg")
    on.exit(unlink(svg_file), add = TRUE)
    render_fallback_svg(svg_file)
  svg_text <- paste(readLines(svg_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

# Add the required explanatory note below the Graphviz canvas. This is SVG
# annotation text, not a node and not part of the causal graph.
note_line_1 <- paste0(
  "Prenatal care utilization was available only among recorded congenital syphilis cases; therefore, its mediating role cannot be evaluated"
)
note_line_2 <- paste0(
  "at the population level. Key individual-level determinants were unavailable, and the maternal-origin contrast should not be interpreted causally."
)

viewbox_match <- regexec(
  'viewBox="([0-9.]+) ([0-9.]+) ([0-9.]+) ([0-9.]+)"', svg_text
)
viewbox_values <- regmatches(svg_text, viewbox_match)[[1]]
if (length(viewbox_values) != 5L) stop("Could not read the exported SVG viewBox.")

vb_x <- as.numeric(viewbox_values[2])
vb_y <- as.numeric(viewbox_values[3])
vb_w <- as.numeric(viewbox_values[4])
vb_h <- as.numeric(viewbox_values[5])
note_space <- 58
new_vb_h <- vb_h + note_space

svg_text <- sub(
  'viewBox="[0-9.]+ [0-9.]+ [0-9.]+ [0-9.]+"',
  sprintf('viewBox="%s %s %s %s"', vb_x, vb_y, vb_w, new_vb_h),
  svg_text
)

svg_note <- sprintf(
  paste0(
    '<g id="figure-note" font-family="Arial, sans-serif" font-size="11" fill="#333333">',
    '<text x="%0.2f" y="%0.2f">%s</text>',
    '<text x="%0.2f" y="%0.2f">%s</text>',
    '</g>'
  ),
  vb_x + 8, vb_y + vb_h + 25, note_line_1,
  vb_x + 8, vb_y + vb_h + 42, note_line_2
)
svg_text <- sub("</svg>", paste0(svg_note, "\n</svg>"), svg_text, fixed = TRUE)

writeLines(svg_text, output_svg, useBytes = TRUE)

# 6000 px across provides at least 300 dpi at widths up to 20 inches.
rsvg::rsvg_png(output_svg, file = output_png, width = 6000)
rsvg::rsvg_pdf(output_svg, file = output_pdf)

stopifnot(
  file.exists(output_svg), file.info(output_svg)$size > 0,
  file.exists(output_png), file.info(output_png)$size > 0,
  file.exists(output_pdf), file.info(output_pdf)$size > 0
)

caption <- paste0(
  "Figure S2. Conceptual DAG of the analytical framework underlying the maternal-origin contrast and the restriction of the secondary analysis to recorded cases.\n\n",
  "Note: Prenatal care utilization was available only among recorded congenital syphilis cases; therefore, its mediating role cannot be evaluated at the population level. Key individual-level determinants were unavailable, and the maternal-origin contrast should not be interpreted causally. DAG, directed acyclic graph."
)

cat(caption, "\n")
cat("\nGenerated files:\n", normalizePath(output_svg), "\n",
    normalizePath(output_png), "\n", normalizePath(output_pdf), "\n", sep = "")

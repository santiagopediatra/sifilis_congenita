#!/usr/bin/env Rscript

# Alternative predominantly vertical layout for the corrected Figure S2 DAG.
# Artwork only: no data are loaded and no statistical model is evaluated.

if (!requireNamespace("rsvg", quietly = TRUE)) {
  install.packages("rsvg", repos = "https://cloud.r-project.org")
}

nodes <- c("period", "origin", "social", "infection", "prenatal",
           "diagnosis", "congenital", "ascertainment", "recorded")
edges <- data.frame(
  from = c("period", "period", "period", "period", "period", "period",
           "origin", "social", "social", "social", "infection", "prenatal",
           "infection", "diagnosis", "congenital", "ascertainment"),
  to = c("origin", "social", "infection", "prenatal", "diagnosis", "ascertainment",
         "social", "infection", "prenatal", "ascertainment", "diagnosis", "diagnosis",
         "congenital", "congenital", "ascertainment", "recorded")
)
stopifnot(length(nodes) == 9L, nrow(edges) == 16L,
          !any(edges$from == "origin" & edges$to == "infection"),
          !any(edges$from == "prenatal" & edges$to == "infection"),
          !any(edges$from == "origin" & edges$to == "recorded"))

svg_file <- "Figure_S2_DAG_corrected_vertical.svg"
png_file <- "Figure_S2_DAG_corrected_vertical.png"
pdf_file <- "Figure_S2_DAG_corrected_vertical.pdf"

grDevices::svg(svg_file, width = 10, height = 13, bg = "white", pointsize = 12)
grid::grid.newpage()

xy <- list(
  period = c(.50, .91), origin = c(.25, .79), social = c(.69, .79),
  infection = c(.30, .63), prenatal = c(.70, .63), diagnosis = c(.50, .49),
  congenital = c(.50, .36), ascertainment = c(.50, .23), recorded = c(.50, .10)
)
wh <- list(
  period = c(.36, .07), origin = c(.24, .06), social = c(.38, .085),
  infection = c(.29, .06), prenatal = c(.29, .06), diagnosis = c(.37, .085),
  congenital = c(.26, .06), ascertainment = c(.36, .085), recorded = c(.30, .075)
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

left <- function(n) xy[[n]][1] - wh[[n]][1] / 2
right <- function(n) xy[[n]][1] + wh[[n]][1] / 2
top <- function(n) xy[[n]][2] + wh[[n]][2] / 2
bottom <- function(n) xy[[n]][2] - wh[[n]][2] / 2

grid::grid.text("Conceptual DAG of the Maternal-Origin Contrast",
                x = .04, y = .985, just = c("left", "top"),
                gp = grid::gpar(fontfamily = "Arial", fontsize = 20, fontface = "bold"))

edge_gp <- grid::gpar(col = "#555555", lwd = 1.7, lineend = "round", linejoin = "round")
path <- function(x, y) {
  grid::grid.lines(x, y, default.units = "npc", gp = edge_gp,
                   arrow = grid::arrow(type = "closed", length = grid::unit(.11, "inches")))
}

# Calendar-period arrows are routed around the outside of the central flow.
path(c(.43, xy$origin[1]), c(bottom("period"), top("origin")))
path(c(.57, xy$social[1]), c(bottom("period"), top("social")))
path(c(.39, .10, .10, left("infection")),
     c(bottom("period"), .855, xy$infection[2], xy$infection[2]))
path(c(.61, .90, .90, right("prenatal")),
     c(bottom("period"), .84, xy$prenatal[2], xy$prenatal[2]))
path(c(.64, .93, .93, right("diagnosis")),
     c(bottom("period"), .825, xy$diagnosis[2], xy$diagnosis[2]))
path(c(.67, .96, .96, right("ascertainment")),
     c(bottom("period"), .81, xy$ascertainment[2], xy$ascertainment[2]))

# Main conceptual flow.
path(c(right("origin"), left("social")), c(.79, .79))
path(c(.62, xy$infection[1]), c(bottom("social"), top("infection")))
path(c(.75, xy$prenatal[1]), c(bottom("social"), top("prenatal")))
path(c(xy$social[1], xy$social[1], .055, .055, left("ascertainment")),
     c(bottom("social"), .72, .72, xy$ascertainment[2], xy$ascertainment[2]))
path(c(xy$infection[1], .42), c(bottom("infection"), top("diagnosis")))
path(c(xy$prenatal[1], .58), c(bottom("prenatal"), top("diagnosis")))
path(c(left("infection"), .12, .12, left("congenital")),
     c(xy$infection[2], xy$infection[2], xy$congenital[2], xy$congenital[2]))
path(c(.50, .50), c(bottom("diagnosis"), top("congenital")))
path(c(.50, .50), c(bottom("congenital"), top("ascertainment")))
path(c(.50, .50), c(bottom("ascertainment"), top("recorded")))

# Nodes are drawn after edges so lines cannot pass visibly through boxes.
for (n in nodes) {
  grid::grid.roundrect(x = xy[[n]][1], y = xy[[n]][2],
                       width = wh[[n]][1], height = wh[[n]][2],
                       r = grid::unit(.06, "snpc"),
                       gp = grid::gpar(fill = fills[n], col = "#333333", lwd = 1.4))
  grid::grid.text(labels[[n]], x = xy[[n]][1], y = xy[[n]][2],
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 10.5, col = "#222222"))
}

note_1 <- paste0("Prenatal care utilization was available only among recorded congenital syphilis cases; therefore, its mediating role cannot be")
note_2 <- paste0("evaluated at the population level. Key individual-level determinants were unavailable, and the maternal-origin contrast should not be interpreted causally.")
grid::grid.text(note_1, x = .04, y = .040, just = c("left", "center"),
                gp = grid::gpar(fontfamily = "Arial", fontsize = 8.5, col = "#333333"))
grid::grid.text(note_2, x = .04, y = .022, just = c("left", "center"),
                gp = grid::gpar(fontfamily = "Arial", fontsize = 8.5, col = "#333333"))

grDevices::dev.off()
rsvg::rsvg_png(svg_file, png_file, width = 3600)
rsvg::rsvg_pdf(svg_file, pdf_file)

stopifnot(all(file.exists(c(svg_file, png_file, pdf_file))),
          all(file.info(c(svg_file, png_file, pdf_file))$size > 0))
cat(normalizePath(png_file), "\n")

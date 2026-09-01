# English STROBE flow diagram for the audited analytic sample.
# The source data are opened read-only; the Spanish outputs are not modified.

source_data <- "SIFILIS_hgoia2009-2024.csv"
denominator_file <- "output/base_anual_auditada.csv"
png_file <- "graficos/STROBE_patient_selection_flowchart_English.png"
pdf_file <- "graficos/STROBE_patient_selection_flowchart_English.pdf"
caption_file <- "graficos/STROBE_patient_selection_flowchart_English_caption.txt"

if (!file.exists(source_data)) stop("Source data not found: ", source_data)
if (!file.exists(denominator_file)) stop("Audited denominator file not found: ", denominator_file)

data <- read.csv(source_data, check.names = FALSE, stringsAsFactors = FALSE)
denominators <- read.csv(denominator_file, check.names = FALSE)
if (!"Extrangera" %in% names(data)) stop("Variable not found: Extrangera")
if (!"nacimientos_total" %in% names(denominators)) stop("Variable not found: nacimientos_total")

origin <- tolower(trimws(as.character(data[["Extrangera"]])))
n_initial <- nrow(data)
n_final <- nrow(data)
n_excluded <- n_initial - n_final
n_national <- sum(origin == "no", na.rm = TRUE)
n_foreign <- sum(origin == "si", na.rm = TRUE)
n_unclassified <- sum(is.na(origin) | !origin %in% c("no", "si"))
n_live_births <- sum(denominators[["nacimientos_total"]], na.rm = TRUE)
pct_national <- round(100 * n_national / n_final, 1)
pct_foreign <- round(100 * n_foreign / n_final, 1)

stopifnot(
  n_live_births == 111742L,
  n_initial == 255L,
  n_excluded == 0L,
  n_final == 255L,
  n_national == 219L,
  n_foreign == 36L,
  pct_national == 85.9,
  pct_foreign == 14.1,
  n_unclassified == 0L,
  n_final == n_national + n_foreign
)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 is required")
dir.create("graficos", recursive = TRUE, showWarnings = FALSE)

caption <- paste0(
  "Figure 1. STROBE flow diagram of the selection of recorded congenital syphilis cases included in the analysis. ",
  "Between January 2009 and June 2024, 111,742 institutional live births and 255 congenital syphilis cases were ",
  "recorded in the Perinatal Information System (SIP). No cases were excluded, and all 255 cases were included in ",
  "the final analytic sample. By maternal origin, 219 (85.9%) cases were born to Ecuadorian mothers and 36 (14.1%) ",
  "to foreign-born mothers."
)
writeLines(caption, caption_file, useBytes = TRUE)

library(ggplot2)

boxes <- data.frame(
  xmin = 1, xmax = 9,
  ymin = c(7.2, 4.5, 1.0), ymax = c(9.2, 6.2, 3.5),
  label = c(
    paste0("Institutional data, January 2009–June 2024\n",
           "Institutional census: live births, n = 111,742\n",
           "Perinatal Information System (SIP): recorded congenital syphilis cases\n",
           "assessed for eligibility, n = ", n_initial),
    paste0("Excluded: n = ", n_excluded,
           "\nNo eligible recorded congenital syphilis cases were excluded",
           "\nfrom the primary analysis"),
    paste0("Primary analytic sample\nn = ", n_final,
           "\nEcuadorian maternal origin: n = ", n_national, " (85.9%)",
           "\nForeign maternal origin: n = ", n_foreign, " (14.1%)")
  )
)

flowchart <- ggplot() +
  geom_rect(data = boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", colour = "black", linewidth = 0.65) +
  geom_text(data = boxes,
            aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
            family = "sans", size = 4.0, lineheight = 1.12) +
  annotate("segment", x = 5, xend = 5, y = 7.2, yend = 6.4,
           colour = "black", linewidth = 0.6,
           arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed")) +
  annotate("segment", x = 5, xend = 5, y = 4.5, yend = 3.7,
           colour = "black", linewidth = 0.6,
           arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed")) +
  coord_cartesian(xlim = c(0.4, 9.6), ylim = c(0.4, 9.8), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(20, 20, 20, 20))

ggsave(png_file, flowchart, width = 7, height = 8, units = "in", dpi = 300, bg = "white")
ggsave(pdf_file, flowchart, width = 7, height = 8, units = "in",
       device = grDevices::cairo_pdf, bg = "white")

stopifnot(all(file.exists(c(png_file, pdf_file, caption_file))),
          all(file.info(c(png_file, pdf_file, caption_file))$size > 0))
cat(normalizePath(png_file), "\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(magick)
  library(officer)
})

dir.create("graficos", recursive = TRUE, showWarnings = FALSE)

# Figure 1: use the saved Poisson predictions; no model is loaded or evaluated.
pred <- read_csv("output/articulo/poisson/predicciones_figura1_poisson.csv",
                 show_col_types = FALSE) %>%
  mutate(origen = factor(origen, levels = c("nacional", "extranjera")))

year_labels <- as.character(2009:2024)
year_labels[2009:2024 == 2024] <- "2024\n(Jan–Jun)"
blue <- "#2166AC"; orange <- "#D6604D"; gray <- "#4D4D4D"

fig1 <- ggplot(pred, aes(anio, tasa_mediana, color = origen, fill = origen)) +
  geom_ribbon(aes(ymin = ICr95_inf, ymax = ICr95_sup),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(aes(y = tasa_observada), size = 2) +
  scale_color_manual(values = c(nacional = blue, extranjera = orange),
                     labels = c(nacional = "Ecuadorian", extranjera = "Foreign")) +
  scale_fill_manual(values = c(nacional = blue, extranjera = orange),
                    labels = c(nacional = "Ecuadorian", extranjera = "Foreign")) +
  scale_x_continuous(breaks = 2009:2024, labels = year_labels) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Hospital-Recorded Congenital Syphilis Cases by Maternal Origin",
    subtitle = "Observed Rates and Posterior Predictions From the Bayesian Poisson Model With a Temporal Function",
    x = "Year",
    y = "Hospital-recorded cases per 1,000 live births",
    color = "Maternal origin", fill = "Maternal origin"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text = element_text(color = gray))

ggsave("graficos/Figure1_English.png", fig1, width = 9, height = 6,
       dpi = 600, bg = "white")

# Figure S1: retain the original raster and replace only its Spanish text.
src_s1 <- image_read("graficos/figura_s1_ppc_poisson.png")
canvas <- image_draw(src_s1)
par(xpd = NA)
rect(0, 0, 4800, 345, col = "white", border = NA)
rect(0, 345, 4800, 465, col = "white", border = NA)
rect(0, 2185, 4800, 2355, col = "white", border = NA)
rect(0, 1110, 155, 3020, col = "white", border = NA)
rect(1500, 3970, 3350, 4200, col = "white", border = NA)
text(325, 95, "Posterior Predictive Checks", adj = c(0, 0.5),
     family = "Arial", cex = 9.0)
text(325, 235, "Orange line: observed statistic", adj = c(0, 0.5),
     family = "Arial", cex = 6.5)
text(2550, 410, "2021 Cases, Ecuadorian Origin", family = "Arial", cex = 5.5)
text(2550, 2280, "Variance of the 32 Counts", family = "Arial", cex = 5.5)
text(80, 2100, "Frequency", family = "Arial", cex = 6.0, srt = 90)
text(2400, 4100, "Replicated Statistic", family = "Arial", cex = 6.5)
dev.off()
image_write(canvas, "graficos/FigureS1_English.png", format = "png")

stopifnot(file.exists("graficos/FigureS2_DAG_vertical.png"))

# Word document: figures and their required captions/notes only.
doc <- read_docx()
sec <- prop_section(page_size = page_size(orient = "landscape"),
                    page_margins = page_mar(top = 0.45, bottom = 0.45,
                                            left = 0.6, right = 0.6))
doc <- body_set_default_section(doc, sec)

add_caption <- function(doc, title, text, note = NULL) {
  doc <- body_add_fpar(doc, fpar(ftext(title, fp_text(bold = TRUE)),
                                 ftext(text, fp_text())))
  if (!is.null(note)) {
    doc <- body_add_fpar(doc, fpar(ftext("Note: ", fp_text(italic = TRUE)),
                                   ftext(note, fp_text())))
  }
  doc
}

doc <- body_add_img(doc, "graficos/Figure1_English.png", width = 8.7, height = 5.8)
doc <- add_caption(doc,
  "Figure 1. ",
  "Observed rates and posterior predictions of hospital-recorded congenital syphilis cases by maternal origin, 2009–June 2024.\nPoints represent observed rates per 1,000 institutional live births. Lines represent posterior expected rates estimated using the Bayesian Poisson regression model with a spline temporal function and the logarithm of live births as an offset; shaded bands represent 95% credible intervals. The maternal-origin contrast was modeled as constant over time. Data for 2024 include January through June only."
)
doc <- body_add_break(doc)
doc <- body_add_img(doc, "graficos/FigureS1_English.png", width = 7.2, height = 6.3)
doc <- add_caption(doc, "Figure S1. ",
  "Posterior predictive checks of the primary Bayesian model.",
  "The orange line represents the observed value; the distribution represents values obtained from data replicated by the model using posterior_predict()."
)
doc <- body_add_break(doc)
doc <- body_add_img(doc, "graficos/FigureS2_DAG_vertical.png", width = 5.3, height = 6.9)
doc <- add_caption(doc, "Figure S2. ",
  "Conceptual DAG for interpretation of the maternal-origin contrast.",
  "Prenatal care utilization was available only among recorded congenital syphilis cases; therefore, its mediating role cannot be evaluated at the population level. Key individual-level determinants were unavailable, and the maternal-origin contrast should not be interpreted causally. DAG, directed acyclic graph."
)

print(doc, target = "graficos/Figures_English_Final.docx")
cat("Created graficos/Figures_English_Final.docx\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(igraph)
  library(ggraph)
  library(grid)
  library(magick)
  library(officer)
})

dir.create("figures/english", recursive = TRUE, showWarnings = FALSE)

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

ggsave("figures/english/Figure1_English.png", fig1, width = 9, height = 6,
       dpi = 600, bg = "white")

# Figure S1: retain the original raster and replace only its Spanish text.
src_s1 <- image_read("figures/articulo/poisson/figura_s1_ppc_poisson.png")
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
image_write(canvas, "figures/english/FigureS1_English.png", format = "png")

# Figure S2: fixed nodes and fixed directed edges from the original source.
nodes <- tribble(
  ~node, ~label, ~x, ~y, ~type,
  "origen", "Maternal\norigin", 0, 2, "exposure",
  "social", "Social determinants\nand access barriers", 2, 3, "unmeasured",
  "infeccion", "Maternal syphilis:\nexposure and treatment", 4, 2.8, "unmeasured",
  "cpn", "Number of prenatal\nvisits", 3, 1.2, "cases_only",
  "anio", "Calendar year /\nmigration / surveillance", 0.5, 0, "measured",
  "resultado", "Hospital-recorded\ncongenital syphilis", 6, 2, "outcome"
)
edges <- tribble(
  ~from, ~to,
  "origen", "social", "social", "infeccion", "social", "cpn",
  "cpn", "infeccion", "infeccion", "resultado", "origen", "infeccion",
  "anio", "origen", "anio", "social", "anio", "infeccion", "anio", "resultado"
)
graph <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
dag <- ggraph(graph, layout = "manual", x = nodes$x, y = nodes$y) +
  geom_edge_link(arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                 end_cap = circle(9, "mm"), start_cap = circle(7, "mm"),
                 linewidth = 0.6, color = "grey35") +
  geom_node_label(aes(x = x, y = y, label = label, fill = type),
                  size = 3.2, linewidth = 0.3,
                  label.padding = unit(0.25, "lines")) +
  scale_fill_manual(values = c(exposure = "#D9EAF7", outcome = "#F7D9D5",
                               measured = "#E5E5E5", cases_only = "#FFF2CC",
                               unmeasured = "#EFE3F7"), guide = "none") +
  coord_equal(xlim = c(-0.7, 6.8), ylim = c(-0.6, 3.7), clip = "off") +
  labs(
    title = "Conceptual DAG of the Maternal-Origin Contrast",
    caption = paste0(
      "The number of prenatal visits is observed only among cases and does not allow assessment of mediation ",
      "at the population level. Key determinants were not measured."
    )
  ) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 8))

ggsave("figures/english/FigureS2_English.png", dag, width = 9, height = 5.5,
       dpi = 600, bg = "white")

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

doc <- body_add_img(doc, "figures/english/Figure1_English.png", width = 8.7, height = 5.8)
doc <- add_caption(doc,
  "Figure 1. ",
  "Observed rates and posterior predictions of hospital-recorded congenital syphilis cases by maternal origin, 2009–June 2024.\nPoints represent observed rates per 1,000 institutional live births. Lines represent posterior expected rates estimated using the Bayesian Poisson regression model with a spline temporal function and the logarithm of live births as an offset; shaded bands represent 95% credible intervals. The maternal-origin contrast was modeled as constant over time. Data for 2024 include January through June only."
)
doc <- body_add_break(doc)
doc <- body_add_img(doc, "figures/english/FigureS1_English.png", width = 7.2, height = 6.3)
doc <- add_caption(doc, "Figure S1. ",
  "Posterior predictive checks of the primary Bayesian model.",
  "The orange line represents the observed value; the distribution represents values obtained from data replicated by the model using posterior_predict()."
)
doc <- body_add_break(doc)
doc <- body_add_img(doc, "figures/english/FigureS2_English.png", width = 8.8, height = 5.38)
doc <- add_caption(doc, "Figure S2. ",
  "Conceptual DAG of the analytical restriction of the secondary analysis among cases.",
  "The number of prenatal visits was observed only among recorded congenital syphilis cases. Consequently, conditioning the analysis on the presence of the outcome allows estimation of the association between maternal origin and the number of prenatal visits within the set of cases, but does not allow assessment of mediation at the population level. DAG, directed acyclic graph."
)

print(doc, target = "SYPHILIS_FIGURES_ENGLISH_FINAL.docx")
cat("Created SYPHILIS_FIGURES_ENGLISH_FINAL.docx\n")

## plotting themes --------------------------------

theme_cowplot2 <- function(...) {
  theme_cowplot(font_size = 16, font_family = "sans", ...) %+replace%
    theme(strip.background = element_blank(),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.border = element_blank())
}
theme_set(theme_cowplot2())

remove_xaxis <- theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.line.x = element_blank())

remove_yaxis <- theme(axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.line.y = element_blank())

remove_guides <- guides(color = F, fill = F, shape = F, alpha = F, size = F)


## ggsave wrapper suppressing dingbats symbols 
## for adobe illustrator compatibility
ggsave_pdf <- function(filename, plot = last_plot(), device = NULL, path = NULL, 
                       scale = 1, width = NA, height = NA, units = "in",#units = c("in", "cm", "mm"), 
                       dpi = 300, limitsize = TRUE, ...) {
  ggsave(filename = filename, plot = plot, #device = cairo_pdf, 
         path = path, 
         scale = scale, width = width, height = height, units = units, 
         dpi = dpi, limitsize = limitsize, ...)
}

ggsave_png <- function(filename, plot = last_plot(), device = NULL, path = NULL, 
                       scale = 1, width = NA, height = NA, units = "in",#units = c("in", "cm", "mm"), 
                       dpi = 300, limitsize = TRUE, type = "cairo", ...) {
  ggsave(filename = filename, plot = plot, device = device, path = path, 
         scale = scale, width = width, height = height, units = units, 
         dpi = dpi, limitsize = limitsize, type = type, ...)
}


## umap helpers --------------------------------------

arrow <- arrow(angle = 20, type = "closed", length = unit(0.1, "npc"))
umap_coord_anno <- ggplot(tibble(group = c("UMAP1", "UMAP2"),
                                 x = c(0, 0), xend = c(1, 0),
                                 y = c(0, 0), yend = c(0, 1),
                                 lx = c(0.5, -0.15), ly = c(-0.15, 0.5),
                                 angle = c(0, 90))) +
  geom_segment(aes(x, y, xend = xend, yend = yend, group = group),
               arrow = arrow, size = 1, lineend = "round") +
  geom_text(aes(lx, ly, label = group, angle = angle), size = 4) +
  theme_void() +
  coord_fixed(xlim = c(-0.3, 1), ylim = c(-0.3, 1))

add_umap_coord <- function(gg_obj) {
  p <- ggdraw() + 
    draw_plot(gg_obj, x = 0, y = 0, width = 1, height = 1) +
    draw_plot(umap_coord_anno, x = -0.015, y = -0.02, width = 0.4, height = 0.4)
  return(p)
}


## cohort marker genes ----------------------------

# markers_v7 <- yaml::read_yaml("resources/annotation/hgsc_v7_major.yaml")
# 
# helper_markers <- function(x) dplyr::select(unnest(enframe(x, "subtype", "gene"), cols = gene), gene, subtype)
# markers_v7_super <- lapply(yaml::read_yaml("resources/annotation/hgsc_v7_super.yaml"), helper_markers)

## load color code --------------------------------

clrs <- yaml::read_yaml("resources/annotation/colors.yaml") %>%
  lapply(function(x) map_depth(x, vec_depth(x)-2, unlist))

# clrs$patient_id_short <- clrs$patient_id
# names(clrs$patient_id_short) <- str_remove_all(names(clrs$patient_id), "SPECTRUM-OV-")

# shps <- yaml::read_yaml("resources/annotation/shapes.yaml") %>% 
#   lapply(function(x) map_depth(x, vec_depth(x)-2, unlist))

## load database ----------------------------------

# db <- readr::read_rds("resources/db/tme/SPECTRUM.rds")

## cell type sort fraction -------------------------

cell_type_super_lookup <- c(
  B.cell = "Immune",
  Plasma.cell = "Immune",
  T.cell = "Immune", 
  Myeloid.cell = "Immune", 
  Mast.cell = "Immune", 
  Dendritic.cell = "Immune", 
  Endothelial.cell = "Stromal",
  Fibroblast = "Stromal", 
  Ovarian.cancer.cell = "Stromal", 
  Ov.cancer.cell = "Stromal", 
  Other = "Other"
)

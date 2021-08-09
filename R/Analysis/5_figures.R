#### ------------ This script generates all the figures in the paper and in appendix S2 ------------ ####


#### ------------ Tibble of intestinal traits (for figures 1 and 2) ------------ ####

# Use posterior predictions
# Predict values of intestinal traits at a common length (SL=15cm)
# Get posterior predictions:
newdata <- tibble(phylo = int_moor_traits_f$phylo, 
                  species = int_moor_traits_f$species,
                  spec_mean_sl_log = int_moor_traits_f$spec_mean_sl_log, 
                  troph = int_moor_traits_f$troph, 
                  elon_log = int_moor_traits_f$elon_log,
                  stomach = int_moor_traits_f$stomach,
                  durophagy = int_moor_traits_f$durophagy,
                  within_spec_sl_log = log(150) - int_moor_traits_f$spec_mean_sl_log)

# Intestinal length
intercept_leng <- as_tibble(fitted(m3_leng, newdata = newdata))
# Intestinal diameter
intercept_diam <- as_tibble(fitted(m3_diam, newdata = newdata))
# Intestinal surface
intercept_surf <- as_tibble(fitted(m3_surf, newdata = newdata))

# Create tibble by combining the traits with the taxonomic information and trophic level
traits <- tibble(phylo = int_moor_traits_f$phylo, 
                 family = int_moor_traits_f$family, 
                 genus = int_moor_traits_f$genus, 
                 species = int_moor_traits_f$species, 
                 leng = intercept_leng$Estimate, 
                 diam = intercept_diam$Estimate, 
                 surf = intercept_surf$Estimate, 
                 troph = int_moor_traits_f$troph, 
                 trophic_guild = int_moor_traits_f$trophic_guild_predicted) %>% 
  unique() %>% 
  arrange(family, genus, species)

# Separate the parrotfishes (tribe Scarini) from other Labridae
traits <- traits %>%
  mutate(family = if_else(genus == "Scarus" | genus == "Chlorurus" | genus == "Calotomus", 
                          true = paste0(family, "\nScarini"), 
			  false = as.character(family)))

# Example of two species with contrasting intestinal morphology for the manuscript
# Get mean intestinal length and diameter (in cm) for Acanthurus guttatus and Cephalopholis argus
example_traits <- traits %>%
  filter(species == "Acanthurus guttatus" | species == "Cephalopholis argus") %>%
  select(species:diam) %>%
  mutate(leng = round(exp(leng)/10, 2), 
         diam = round(exp(diam)/10, 2))
example_traits
# species              leng  diam
# Acanthurus guttatus  95.55 0.72
# Cephalopholis argus  14.26 0.25


#### ------------ Figure 1 - phylogenetic tree ------------ ####

# Use "tree_int", which include only species with valid phylogenetic information (n=139)

# Select tip labels
tip_label <- data.frame(phylo = tree_int$tip.label)
tip_label <- tip_label %>%
  rename(phylo = tree_int.tip.label)

# Join traits with tip labels to order the data frame
traits_phylo <- left_join(tip_label, traits)
rownames(traits_phylo) <- traits_phylo$phylo # need for gheatmap()

# Plot circular tree
tree_plot <- ggtree(tree_int, layout = "circular")

# Join traits data with tree view
tree_plot <- tree_plot %<+% traits_phylo

# Process the updated tree view data
# Need this step to be able to place segments and family labels around the tree view
# Alternatively use ggtree::geom_cladelabel(), but need to repeat for each family
tree_dt <- tree_plot$data
head(tree_dt)

# Select only the tip labels and order by coord y
tree_dt <- tree_dt %>%
  filter(isTip == TRUE) %>%
  arrange(y)

# Make table with y coords for each family in consecutive order
# This helps for drawing and labelling segments
coord_families <- tree_dt %>%
  group_by(family) %>%
  mutate(y1 = y[1],
         y2 = y[n()],
         angle = mean(angle),
         n = n()) %>%
  ungroup() %>%
  select(family, y1, y2, angle, n) %>%
  unique()

coord_families

# Compute the middle y - will be used for placing the family label
# The mean angle was computed above
coord_families <- coord_families %>%
  mutate(y_mid = rowMeans(select(., y1:y2)))

# Adjust y coordinates so that a segment is also drawn for one-species families where y1=y2
# If not, no segments are drawn for such cases
# To force a segment add and subtract a small amount
# Do this for all families to keep the same distance among all segments
coord_families <- coord_families %>%
  mutate(y1_adj = y1 - 0.3,
         y2_adj = y2 + 0.3)

# Adjust label's angle for cases between 90 and 270 degrees
coord_families <- coord_families %>%
  mutate(angle_adj = if_else(between(angle, 90, 180), 
                             true = angle + 180,
                             false = if_else(angle > 180 & angle <= 270,
                               		     true = angle - 180,
                                             false = angle)))

# Make horizontal adjustment from 0 to 1 for labels with angles between 90 and 270 degrees
coord_families <- coord_families %>%
  mutate(hjust_adj = if_else(between(angle, 90, 270), true = 1, false = 0))

coord_families

# Define variable to control x coordinate of segments and labels
my_x <- max(tree_dt$x) + 32

# Add tip points for trophic level
tree_plot1 <- tree_plot +
  geom_tippoint(aes(x = x + 1.6, fill = troph), shape = 21, size = 1.5) + 
  scale_fill_gradient(low = "white", high = "black", "Trophic<br>level",
                      guide = guide_colorbar(ticks = FALSE, order = 1, title.position = "top", title.hjust = 0.5), 
                      breaks = c(2, 3, 4))

# Add heatmaps of traits
# Add heatmap 1 - intestinal length
tree_plot2 <- tree_plot1 + new_scale_fill()
tree_plot2 <- gheatmap(tree_plot2, traits_phylo[, "leng", drop = F], offset = 5, width = .05, colnames = 0) + 
  scale_fill_gradientn("*ln* intestinal<br>length (mm)", colors = bluecols(100), 
                       guide = guide_colorbar(ticks = FALSE, order = 2, title.position = "top", title.hjust = 0.5))

# Add heatmap 2 - intestinal diameter
tree_plot3 <- tree_plot2 + new_scale_fill()
tree_plot3 <- gheatmap(tree_plot3, traits_phylo[, "diam", drop = F], offset = 12, width = .05, colnames = 0) + 
  scale_fill_gradientn("*ln* intestinal<br>diameter (mm)", colors = redcols(100), 
                       guide = guide_colorbar(ticks = FALSE, order = 3, title.position = "top", title.hjust = 0.5))

# Add heatmap 3 - intestinal surface
tree_plot4 <- tree_plot3 + new_scale_fill()
tree_plot4 <- gheatmap(tree_plot4, traits_phylo[, "surf", drop = F], offset = 19, width = .05, colnames = 0) + 
  scale_fill_gradientn("*ln* intestinal<br>surface (mm<sup>2</sup>)", colors = purplecols(100), 
                       guide = guide_colorbar(ticks = FALSE, order = 4, title.position = "top", title.hjust = 0.5), 
                       breaks = c(5, 7, 9))

# Add segments and family labels
tree_plot5 <- tree_plot4 +
  # Add segments for each family
  geom_segment(data = coord_families, 
               aes(x = my_x, y = y1_adj, xend = my_x, yend = y2_adj), 
               color = "black", lineend = "butt", size = 1) +
  # Add family labels
  geom_text(data = coord_families,
            aes(x = my_x, y = y_mid, angle = angle_adj, hjust = hjust_adj, label = family, family = "serif"),
            size = 2.8, nudge_x = 5, color = "black")

# Adjust legend
tree_plot6 <- tree_plot5 +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_markdown(size = 9, face = "bold", family = "serif"),
        legend.text = element_text(size = 8, family = "serif"),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.spacing = unit(-2, "cm"), # reduce space between legend and plot
        plot.margin = unit(c(t = -2, r = -2, b = 0.1, l = -2), "cm"))

# Create tibble with silhouettes for most speciose families (n > 2)
# Use silhouettes from the package "fishualize" (Schiettekatte et al. 2019)
images <- tibble(family = c("Acanthuridae", "Balistidae", "Chaetodontidae", "Cirrhitidae", "Holocentridae", "Labridae", "Labridae\nScarini", "Lutjanidae", "Monacanthidae", "Mullidae", "Pomacanthidae", "Pomacentridae", "Serranidae", "Tetraodontidae"), 
                 image = c("https://github.com/simonjbrandl/fishape/raw/master/shapes/Acanthuridae_Zebrasoma.scopas.png",
                           "https://github.com/simonjbrandl/fishape/raw/master/shapes/Balistidae_Balistapus.undulatus.png",
                           "https://github.com/simonjbrandl/fishape/raw/master/shapes/Chaetodontidae_Chaetodon.vagabundus.png",
                           "https://github.com/simonjbrandl/fishape/raw/master/shapes/Cirrhitidae_Paracirrhites.forsteri.png",
                           "https://github.com/simonjbrandl/fishape/raw/master/shapes/Holocentridae_Sargocentron.spiniferum.png",
                           "https://github.com/simonjbrandl/fishape/raw/master/shapes/Labridae_Epibulus.insidiator.png",
                           "https://github.com/simonjbrandl/fishape/raw/master/shapes/Labridae_Chlorurus.microrhinus.png",
                           "https://github.com/simonjbrandl/fishape/raw/master/shapes/Lutjanidae_Lutjanus.gibbus.png",
                           "https://github.com/simonjbrandl/fishape/raw/master/shapes/Monacanthidae_Aluterus.scriptus.png",
                           "https://github.com/simonjbrandl/fishape/raw/master/shapes/Mullidae_Parupeneus.multifasciatus.png",
                           "https://github.com/simonjbrandl/fishape/raw/master/shapes/Pomacanthidae_Pomacanthus.imperator.png",
                           "https://github.com/simonjbrandl/fishape/raw/master/shapes/Pomacentridae_Abudefduf.sexfasciatus.png",
                           "https://github.com/simonjbrandl/fishape/raw/master/shapes/Serranidae_Cephalopholis.argus.png",
                           "https://github.com/simonjbrandl/fishape/raw/master/shapes/Tetraodontidae_Arothron.meleagris.png"),
                 image_size = c(0.035, 0.045, 0.035, 0.04, 0.05, 0.045, 0.065, 0.05, 0.05, 0.055, 0.04, 0.035, 0.06, 0.055),
                 x = c(my_x+70, my_x+65, my_x+80, my_x+60, my_x+80, my_x+55, my_x+60, my_x+65, my_x+85, my_x+60, my_x+75, my_x+75, my_x+60, my_x+85))

# Add mean y coord
images <- left_join(images, coord_families[, c("family", "y_mid")])

# Add silhouettes to the plot
tree_plot7 <- tree_plot6 + 
  ggimage::geom_image(data = images, aes(x = x, y = y_mid, image = image), size = images$image_size)

# Save Fig. 1
ggsave("Output/Figures/Fig. 1.png", plot = tree_plot7, width = 16, height = 16, units = "cm", dpi = 600, type = "cairo")



#### ------------ Figure 2 - morphospace ------------ ####

## Fig. 2a
traits <- traits %>%
  group_by(family) %>%
  mutate(n = n()) %>%
  ungroup()

# Calculate average intestinal length and diameter (i.e. model intercept) at SL=15cm
# These will be placed as horizontal and vertical lines in the morphospace
# Use the most common categories for stomach and durophagy
newdata_mean <- tibble(spec_mean_sl_log = mean(int_moor_traits_f$spec_mean_sl_log), 
                       troph = mean(int_moor_traits_f$troph),
                       elon_log = mean(int_moor_traits_f$elon_log),
                       stomach = "Present",
                       durophagy = "Not durophagous",
                       within_spec_sl_log = log(150) - mean(int_moor_traits_f$spec_mean_sl_log))
# Intestinal length
intercept_leng_mean <- fitted(m3_leng, newdata = newdata_mean, re_formula = NA)
# Intestinal diameter
intercept_diam_mean <- fitted(m3_diam, newdata = newdata_mean, re_formula = NA)

# Colors for families with at least 3 species based on fish palette "Scarus_quoyi" from package "fishualize" (Schiettekatte et al. 2019)
morphospace_family_cols <- fishualize::fish(option = "Scarus_quoyi", n = 14, direction = -1)

# Create tibble with x and y coordinates for text, segment and image for each family to annotate the morphospace
# No segment for Balistidae, Chaetodontidae, Labridae and Tetraodontidae
morphospace_family_coord <- tibble(family = images$family,
                                   x_text = c(2.2, 1.8, 0.8, 0.5, 0.25, 1.55, 2.3, 0.9, 1.97, 1.2, 1.65, 0.65, 0.65, 2.2),
                                   y_text = c(7.2, 5.45, 7.1, 3.9, 5.2, 4.1, 5, 3.9, 4.5, 3.7, 7.2, 6.1, 5, 6.45),
                                   x_segm = c(2.03, NA, NA, 0.5, 0.43, NA, 2.3, 0.9, 1.97, 1.2, 1.65, 0.84, 0.79, NA),
                                   y_segm = c(7.15, NA, NA, 4, 5.18, NA, 5.25, 4, 4.6, 3.8, 7.1, 6.08, 4.99, NA),
                                   xend_segm = c(1.99, NA, NA, 0.98, 1.15, NA, 2.26, 1.01, 1.92, 1.2, 1.61, 1.21, 0.88, NA),
                                   yend_segm = c(6.92, NA, NA, 4.96, 5.2, NA, 5.93, 4.75, 5.85, 4.6, 6.77, 5.9, 4.98, NA),
                                   image = images$image,
                                   image_size = c(0.05, 0.06, 0.05, 0.06, 0.07, 0.065, 0.08, 0.07, 0.08, 0.08, 0.06, 0.05, 0.08, 0.075))

morphospace_family_coord <- morphospace_family_coord %>%
  mutate(x_image = x_text,
         y_image = y_text + c(0.3, - 0.3, 0.25, -0.3, -0.3, -0.3, -0.5, -0.3, -0.3, -0.3, 0.25, 0.25, -0.3, 0.25))

# Plot morphospace with families
morphospace_family <- ggplot() +
  geom_vline(xintercept = intercept_diam_mean[1,1], linetype = 2, color = "grey50") + # add vertical line for mean diameter (model intercept)
  geom_hline(yintercept = intercept_leng_mean[1,1], linetype = 2, color = "grey50") + # add horizontal line for mean length (model intercept)
  geom_point(data = subset(traits, n < 3), aes(x = diam, y = leng, size = surf), shape = 16, color = "grey") + # add gray dots for families with less than 3 species (n = 18)
  ggpubr::stat_chull(data = subset(traits, n > 2), aes(x = diam, y = leng, fill = family), alpha = 0.2, geom = "polygon") + # add polygons for families with at least 3 species (n = 13)
  geom_point(data = subset(traits, n > 2), aes(x = diam, y = leng, color = family, size = surf), shape = 16, alpha = 0.9) + # add dots for families with at least 3 species (n = 13)
  scale_size(range = c(0, 3), limits = c(3.5, 10), breaks = c(4, 6, 8, 10)) +
  labs(x = "*ln* intestinal diameter (mm)", y = "*ln* intestinal length (mm)", size = "*ln* intestinal<br>surface (mm<sup>2</sup>)") + # use package "ggtext"
  scale_fill_manual(values = morphospace_family_cols) + scale_color_manual(values = morphospace_family_cols) +
  guides(color = FALSE, fill = FALSE) +
  xlim(0.11, 2.36) + ylim(3, 7.65) +
  theme(legend.position = c(0.1, 0.8), # top-left
        legend.title = element_markdown(size = 10, hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.background = element_blank(),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown())

# Add family names
morphospace_family <- morphospace_family +
  geom_text(data = morphospace_family_coord, 
            aes(x = x_text, y = y_text, label = family, family = "serif"), 
            size = 4, color = "black")
  
# Add arrows
morphospace_family <- morphospace_family +
  geom_segment(data = morphospace_family_coord, 
               aes(x = x_segm, y = y_segm, xend = xend_segm, yend = yend_segm),
               color = "black", lineend = "butt", size = 0.4)

# Add silhouettes
morphospace_family <- morphospace_family +
  ggimage::geom_image(data = morphospace_family_coord, 
                      aes(x = x_image, y = y_image, image = image), 
                      size = morphospace_family_coord$image_size, color = morphospace_family_cols, asp = 1.7)

# Align legend
morphospace_family <- cowplot::ggdraw(align_legend(morphospace_family))


## Fig. 2b
# Add names and colors to trophic guilds
trophic_guild_categories <- tibble(trophic_guild = c(1:8),
                                   trophic_guild_name = c("Sessile invertivores", "HMD", "Corallivores", "Piscivores", "Microinvertivores", "Macroinvertivores", "Crustacivores", "Planktivores"),
                                   color = c("hotpink3", "chartreuse2", "slateblue4", "orangered", "lightskyblue", "burlywood", "purple1", "yellow"))
# Add trophic guild names to traits
traits <- left_join(traits, trophic_guild_categories[, 1:2])

# Plot morphospace with trophic guilds
morphospace_trophic_guild <- ggplot() +
  geom_vline(xintercept = intercept_diam_mean[1,1], linetype = 2, color = "grey50") + # add vertical line for mean diameter (model intercept)
  geom_hline(yintercept = intercept_leng_mean[1,1], linetype = 2, color = "grey50") + # add horizontal line for mean length (model intercept)
  geom_point(data = subset(traits, is.na(trophic_guild_name)), aes(x = diam, y = leng, size = surf), shape = 16, color = "grey") + # add gray dots for species without a predicted trophic guild
  ggpubr::stat_chull(data = subset(traits, !is.na(trophic_guild_name)), aes(x = diam, y = leng, fill = trophic_guild_name), alpha = 0.2, geom = "polygon") + # add polygons for trophic guilds
  geom_point(data = subset(traits, !is.na(trophic_guild_name)), aes(x = diam, y = leng, color = trophic_guild_name, size = surf), shape = 16, alpha = 0.9) + # add dots for species with predicted trophic guild
  scale_size(range = c(0, 3), limits = c(3.5, 10)) +
  labs(x = "*ln* intestinal diameter (mm)", y = "*ln* intestinal length (mm)", color = "Trophic guild") + # use package "ggtext"
  scale_color_manual(values = trophic_guild_categories$color, 
                     breaks = trophic_guild_categories$trophic_guild_name) +
  scale_fill_manual(values = trophic_guild_categories$color, 
                    breaks = trophic_guild_categories$trophic_guild_name) +
  guides(fill = FALSE, size = FALSE) +
  xlim(0.11, 2.36) + ylim(3, 7.65) +
  theme(legend.position = c(0.12, 0.74), # top-left
        legend.title = element_text(size = 10, hjust = 0.5),
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill = alpha('white', 0.6)),
        legend.key.size = unit(1, "line"),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        axis.text = element_text())

# Align legend
morphospace_trophic_guild <- cowplot::ggdraw(align_legend(morphospace_trophic_guild))


# Arrange the plots
fig_2 <- ggpubr::ggarrange(morphospace_family, morphospace_trophic_guild, labels = c("a", "b"),
                          font.label = list(size = 13), nrow = 2, ncol = 1, align = c("v"))
 
# Save Fig. 2
ggsave("Output/Figures/Fig. 2.png", plot = fig_2, width = 18, height = 22, units = "cm", dpi = 600, type = "cairo")



#### ------------ Figure 3 - relationships between intestinal traits and elongation and trophic level ------------ ####

# Fig. 3a - intestinal length vs elongation
fig_3a <- plot_fig3_elon(m3_leng, int_moor_traits_f, 1000, xlab = "", ylab = "Intestinal length (mm)", 
                         lim = c(0, 800), col1 = bluecols(5)[5], col2 = bluecols(5)[3]) +
  theme(plot.margin = margin(t = 0.1, r = 0, b = 0, l = 0.1))

# Fig. 3b - intestinal length vs trophic level
fig_3b <- plot_fig3_troph(m3_leng, int_moor_traits_f, 1000, xlab = "", ylab = "", 
                          lim = c(0, 800), col1 = bluecols(5)[5], col2 = bluecols(5)[3]) +
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0, l = 0))

# Fig. 3c - intestinal diameter vs elongation
fig_3c <- plot_fig3_elon(m3_diam, int_moor_traits_f, 1000, xlab = "", ylab = "Intestinal diameter (mm)", 
                         lim = c(0, 7), col1 = redcols(5)[5], col2 = redcols(5)[3]) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0.1))

# Fig. 3d - intestinal diameter vs trophic level
fig_3d <- plot_fig3_troph(m3_diam,int_moor_traits_f, 1000,  xlab = "", ylab = "", 
                          lim = c(0, 7), col1 = redcols(5)[5], col2 = redcols(5)[3]) +
  theme(plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0))

# Fig. 3e - intestinal surface vs elongation
fig_3e <- plot_fig3_elon(m3_surf, int_moor_traits_f, 1000, xlab = "Elongation", ylab = "Intestinal surface (mm<sup>2</sup>)", 
                         lim = c(0, 15000), col1 = purplecols(5)[5], col2 = purplecols(5)[3]) +
  theme(axis.title.y = element_markdown(),
        plot.margin = margin(t = 0, r = 0, b = 0.1, l = 0.1))

# Fig. 3f - intestinal surface vs trophic level
fig_3f <- plot_fig3_troph(m3_surf, int_moor_traits_f, 1000,  xlab = "Trophic level", ylab = "", 
                          lim = c(0, 15000), col1 = purplecols(5)[5], col2 = purplecols(5)[3]) +
  theme(plot.margin = margin(t = 0, r = 0.1, b = 0.1, l = 0))

# Arrange the plots
fig_3 <- ggpubr::ggarrange(fig_3a, fig_3b, fig_3c, fig_3d, fig_3e, fig_3f, 
                          labels = c("a", "b", "c", "d", "e", "f"),
                          font.label = list(size = 12),
                          label.x = 0.9, nrow = 3, ncol = 2, align = c("hv"))

# Save Fig. 3
ggsave("Output/Figures/Fig. 3.png", plot = fig_3, width = 16, height = 18, units = "cm", dpi = 600, type = "cairo")



#### ------------ Figure 4 - stomach and durophagy ------------ ####

# INTESTINAL LENGTH
fig_4a <- plot_fig4(int_moor_traits_f, m3_leng, "*ln* intestinal length (mm)")

# INTESTINAL DIAMETER
fig_4b <- plot_fig4(int_moor_traits_f, m3_diam, "*ln* intestinal diameter (mm)")

# INTESTINAL SURFACE
fig_4c <- plot_fig4(int_moor_traits_f, m3_surf, "*ln* intestinal surface (mm<sup>2</sup>)")

fig_4 <- ggarrange(fig_4a, 
                   fig_4b + rremove("y.text") + rremove("y.ticks"),
                   fig_4c + rremove("y.text") + rremove("y.ticks"),
                   nrow = 1, ncol = 3, common.legend = TRUE, widths = c(1, 0.8, 0.8), align = "h",
                   labels = c("a", "b", "c"), font.label = list(size = 12), label.x = c(0.88, 0.85, 0.85), label.y = 0.97)

# Save plot
ggsave("Output/Figures/Fig. 4.png", plot = fig_4 , width = 17, height = 8, units = "cm", dpi = 600, type = "cairo")



#### ------------ Figure 5 - species-specific slopes ------------ ####

# Use the tibbles obtained with the selected slopes
# Need to add three columns to each tibble:
# - a column TRAIT to create 3 separate panel with facet_grid
# - a column ISOMETRY with the respective value of isometry (1 for length and diameter, 2 for surface)
# - a column ALLOMETRY with the type of relationship (isometric or allometric) to color CIs

# INTESTINAL LENGTH
slope_leng_selection <- slope_leng_selection %>%
  mutate(trait = "Intestinal length",
         isometry = 1,
         allometry = if_else(.upper < 1, true = "Neg", # negative allometry
                             false = if_else(.lower > 1, true = "Pos", # positive allometry
                                            false = "Iso") # isometry
         )
  )

# INTESTINAL DIAMETER
slope_diam_selection <- slope_diam_selection %>%
  mutate(trait = "Intestinal diameter",
         isometry = 1,
         allometry = if_else(.upper < 1, true = "Neg", # negative allometry
                             false = if_else(.lower > 1, true = "Pos", # positive allometry
                                            false = "Iso") # isometry
         )
  )

# INTESTINAL SURFACE
slope_surf_selection <- slope_surf_selection %>%
  mutate(trait = "Intestinal surface",
         isometry = 2,
         allometry = if_else(.upper < 2, true = "Neg", # negative allometry
                             false = if_else(.lower > 2, true = "Pos", # positive allometry
                                             false = "Iso") # isometry
         )
  )

# Plot
fig_5 <- rbind(slope_leng_selection, slope_diam_selection, slope_surf_selection) %>%
  
  # Need to reorder the traits (default is in alphabetical order)
  mutate(trait = factor(trait, levels = c("Intestinal length", "Intestinal diameter", "Intestinal surface"))) %>%
  ggplot() +
  geom_vline(aes(xintercept = isometry), linetype = 2) +
  tidybayes::geom_pointinterval(aes(x = spec_slope, y = forcats::fct_rev(species), xmin = .lower, xmax = .upper, color = allometry),
                                interval_size_range = c(0.4, 0.9), fatten_point = 1.5) +
  facet_grid(.~trait, scales = "free") +
  scale_color_manual(breaks = c("Iso", "Neg", "Pos"), values = c("grey50", "firebrick", "dodgerblue4"), guide = FALSE) +
  labs(x = "Slope", y = "") + 
  theme(panel.grid.major = element_line(color = "grey90"),
        axis.text.y = element_text(face = "italic", size = 8),
        axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10),
        strip.text.x = element_text(face = "bold", size = 10),
        strip.background = element_blank())

# Save plot
ggsave("Output/Figures/Fig. 5.png", plot = fig_5, width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")



#### ------------ Appendix S1: figure 1 - delta15N vs trophic level ------------ ####

# Plot
TL_vs_d15n <- m1_d15n %>% 
  tidybayes::spread_draws(b_Intercept, b_scaletroph, n = 1000, seed = 29) %>% # extract 1000 random draws
  mutate(troph = list(seq(min(d15n_moor$troph), max(d15n_moor$troph), 0.01))) %>% # the observed value range of trophic level
  tidyr::unnest(troph) %>%
  mutate(pred = b_Intercept + (b_scaletroph/sd(d15n_moor$troph))*(troph - mean(troph))) %>% # predict using unscaled slope and centred trophic level
  group_by(troph) %>%
  mutate(pred_m = mean(pred, na.rm = TRUE)) %>%
  ggplot(aes(x = troph)) +
  geom_line(aes(y = pred, group = .draw), color = "grey80", alpha = 0.2) + # add 1000 predictions to show uncertainty around the mean
  geom_point(data = d15n_moor, aes(x = troph, y = delta15N), color = "grey30", size = 0.6, alpha = 0.8, shape = 1) + # add raw data
  geom_line(aes(y = pred_m), color = "black") + # add mean prediction
  ylab(expression(delta^15*"N")) + xlab("Trophic level")

# Save figure
ggsave("Output/Figures/TL_vs_d15n.png", plot = TL_vs_d15n, width = 12, height = 8, units = "cm", dpi = 600, type = "cairo")



#### ------------ Appendix S2: Figure 1 - posterior predictive checks ------------ ####

# Intestinal length
set.seed(29)
pp_check_leng <- pp_check(m3_leng, type = "dens_overlay", nsamples = 100) +
  xlab("*ln* intestinal length (mm)") +
  scale_color_manual(values = c(bluecols(5)[c(5,3)])) +
  theme(axis.title.x = element_markdown())

# Intestinal diameter
set.seed(29)
pp_check_diam <- pp_check(m3_diam, type = "dens_overlay", nsamples = 100) + 
  xlab("*ln* intestinal diameter (mm)") +
  scale_color_manual(values = c(redcols(5)[c(5,3)])) +
  theme(axis.title.x = element_markdown())

# Intestinal surface
set.seed(29)
pp_check_surf <- pp_check(m3_surf, type = "dens_overlay", nsamples = 100) +
  xlab("*ln* intestinal surface (mm<sup>2</sup>)") +
  scale_color_manual(values = c(purplecols(5)[c(5,3)])) +
  theme(axis.title.x = element_markdown())

# Arrange the plots
pp_check_plot <- ggpubr::ggarrange(pp_check_leng, pp_check_diam, pp_check_surf, 
                                   legend = "none", nrow = 1, ncol = 3, align = "h",
                                   labels = c("a", "b", "c"), font.label = list(size = 11), label.x = 0.1, label.y = 0.97)
# Add title
pp_check_plot <- ggpubr::annotate_figure(pp_check_plot, top = ggpubr::text_grob("Posterior predictive draws (thin lines) vs observed data (thick line)", color = "black", size = 13, family = "serif"))

# Save plot
ggsave("Output/Figures/pp_check.png", plot = pp_check_plot, width = 17, height = 8, units = "cm", dpi = 600, type = "cairo")



#### ------------ Appendix S2: Figure 2 - sensitivity plot trophic level ------------ ####

# Plot relationship between intestinal traits and trophic level for models of sensitivity analysis

# INTESTINAL LENGTH
# Model with 5 individuals
sens_leng_tl_5 <- plot_fig3_troph(m_leng_5, int_moor_traits_f_5, 1000, xlab = "", ylab = "", 
                                  lim = c(0, 800), col1 = bluecols(5)[5], col2 = bluecols(5)[3]) +
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0, l = 0))

# Model with 8 individuals
sens_leng_tl_8 <- plot_fig3_troph(m_leng_8, int_moor_traits_f_8, 1000, xlab = "", ylab = "", 
                                  lim = c(0, 800), col1 = bluecols(5)[5], col2 = bluecols(5)[3]) +
  theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0, l = 0))

# INTESTINAL DIAMETER
# Model with 5 individuals
sens_diam_tl_5 <- plot_fig3_troph(m_diam_5,int_moor_traits_f_5, 1000,  xlab = "", ylab = "", 
                                  lim = c(0, 7), col1 = redcols(5)[5], col2 = redcols(5)[3]) +
  theme(plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0))

# Model with 8 individuals
sens_diam_tl_8 <- plot_fig3_troph(m_diam_8,int_moor_traits_f_8, 1000,  xlab = "", ylab = "", 
                                  lim = c(0, 7), col1 = redcols(5)[5], col2 = redcols(5)[3]) +
  theme(plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0))

# INTESTINAL SURFACE
# Model with 5 individuals
sens_surf_tl_5 <- plot_fig3_troph(m_surf_5, int_moor_traits_f_5, 1000,  xlab = "Trophic level", ylab = "", 
                                  lim = c(0, 15000), col1 = purplecols(5)[5], col2 = purplecols(5)[3]) +
  theme(plot.margin = margin(t = 0, r = 0.1, b = 0.1, l = 0))

# Model with 8 individuals
sens_surf_tl_8 <- plot_fig3_troph(m_surf_8, int_moor_traits_f_8, 1000,  xlab = "Trophic level", ylab = "", 
                                  lim = c(0, 15000), col1 = purplecols(5)[5], col2 = purplecols(5)[3]) +
  theme(plot.margin = margin(t = 0, r = 0.1, b = 0.1, l = 0))

# Arrange the plots
sensitivity_plot_TL <- ggpubr::ggarrange(fig_3b + ylab("Intestinal length (mm)"), sens_leng_tl_5, sens_leng_tl_8, 
                                         fig_3d + ylab("Intestinal diameter (mm)"), sens_diam_tl_5, sens_diam_tl_8,
                                         fig_3f + ylab("Intestinal surface (mm<sup>2</sup>)") + theme(axis.title.y = element_markdown()), sens_surf_tl_5, sens_surf_tl_8,
                                         labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i"),
                                         label.x = 0.9, font.label = list(size = 12),
                                         nrow = 3, ncol = 3, align = c("hv"))

# Save plot
ggsave("Output/Figures/sensitivity_plot_TL.png", plot = sensitivity_plot_TL, width = 17, height = 20, units = "cm", dpi = 600, type = "cairo")



#### ------------ Appendix S2: Figure 3 - sensitivity plot stomach durophagy ------------ ####

# Plot stomach and durophagy effects for models of sensitivity analysis

# INTESTINAL LENGTH
# Model with 5 individuals
sens_leng_st_du_5 <- plot_fig4(int_moor_traits_f_5, m_leng_5, "") +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0.1))

# Model with 8 individuals
sens_leng_st_du_8 <- plot_fig4(int_moor_traits_f_8, m_leng_8, "*ln* intestinal length (mm)") +
  theme(plot.margin = margin(t = 0, r = 0, b = 0.1, l = 0.1))

# INTESTINAL DIAMETER
# Model with 5 individuals
sens_diam_st_du_5 <- plot_fig4(int_moor_traits_f_5, m_diam_5, "") +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

# Model with 8 individuals
sens_diam_st_du_8 <- plot_fig4(int_moor_traits_f_8, m_diam_8, "*ln* intestinal diameter (mm)") +
  theme(plot.margin = margin(t = 0, r = 0, b = 0.1, l = 0))

# INTESTINAL SURFACE
# Model with 5 individuals
sens_surf_st_du_5 <- plot_fig4(int_moor_traits_f_5, m_surf_5, "") +
  theme(plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0))

# Model with 8 individuals
sens_surf_st_du_8 <- plot_fig4(int_moor_traits_f_8, m_surf_8, "*ln* intestinal surface (mm<sup>2</sup>)") +
  theme(plot.margin = margin(t = 0, r = 0.1, b = 0.1, l = 0))

# Arrange the plots
sensitivity_plot_st_du <- ggpubr::ggarrange(fig_4a + rremove("xlab") + theme(plot.margin = margin(t = 0.1, r = 0, b = 0, l = 0.1)),  
                                            fig_4b + rremove("xlab") + rremove("y.text") + rremove("y.ticks") + theme(plot.margin = margin(t = 0.1, r = 0, b = 0, l = 0)),
                                            fig_4c + rremove("xlab") + rremove("y.text") + rremove("y.ticks") + theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0, l = 0)),
                                            sens_leng_st_du_5, 
                                            sens_diam_st_du_5 + rremove("y.text") + rremove("y.ticks"),
                                            sens_surf_st_du_5 + rremove("y.text") + rremove("y.ticks"), 
                                            sens_leng_st_du_8, 
                                            sens_diam_st_du_8 + rremove("y.text") + rremove("y.ticks"), 
                                            sens_surf_st_du_8 + rremove("y.text") + rremove("y.ticks"),
                                            labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i"),
                                            font.label = list(size = 12), label.x = 0.9,
                                            common.legend = TRUE, legend = "top",
                                            nrow = 3, ncol = 3, align = c("hv"))

# Save plot
ggsave("Output/Figures/sensitivity_plot_st_du.png", plot = sensitivity_plot_st_du, width = 17, height = 15, units = "cm", dpi = 600, type = "cairo")

 

#### ------------ Appendix S2: Figure 4 - phylogenetic effect ------------ ####

# Extract phylogenetic effect for each species

# INTESTINAL LENGTH
phylo_leng <- m3_leng %>%
  tidybayes::spread_draws(r_phylo[species,]) %>% 
  group_by(species) %>% 
  tidybayes::median_qi(.width = c(.50, .95)) %>% 
  ungroup() %>%
  mutate(species = gsub("_", " ", species),
         trait = "Intestinal length") %>%
  left_join(traits[, c("species", "family")])

# INTESTINAL DIAMETER
phylo_diam <- m3_diam %>%
  tidybayes::spread_draws(r_phylo[species,]) %>% 
  group_by(species) %>% 
  tidybayes::median_qi(.width = c(.50, .95)) %>% 
  ungroup() %>%
  mutate(species = gsub("_", " ", species),
         trait = "Intestinal diameter") %>%
  left_join(traits[, c("species", "family")])

# INTESTINAL SURFACE
phylo_surf <- m3_surf %>%
  tidybayes::spread_draws(r_phylo[species,]) %>% 
  group_by(species) %>% 
  tidybayes::median_qi(.width = c(.50, .95)) %>% 
  ungroup() %>%
  mutate(species = gsub("_", " ", species),
         trait = "Intestinal surface") %>%
  left_join(traits[, c("species", "family")])


# Plot 
phylo_effect_plot <- bind_rows(phylo_leng, phylo_diam, phylo_surf) %>%
  
  # Rename Acanthurus nigroris as A. nigros
  mutate(species = recode(species, "Acanthurus nigroris" = "Acanthurus nigros")) %>%
  
  # Need to order the traits (default is in alphabetical order)
  mutate(trait = factor(trait, levels = c("Intestinal length", "Intestinal diameter", "Intestinal surface"))) %>%
  ggplot() +
  tidybayes::geom_pointinterval(aes(x = r_phylo, y = forcats::fct_rev(species), xmin = .lower, xmax = .upper),
                                interval_size_range = c(0.3, 0.7), fatten_point = 1.2) +
  facet_grid(family~trait, scales = "free", space = "free_y") +
  labs(x = "Deviation from the global intercept on natural-log scale", y = "") + 
  geom_vline(xintercept = 0) +
  theme(panel.grid.major = element_line(color = "grey90"),
        axis.text.y = element_text(face = "italic", color = "black", size = 6.5),
        axis.text.x = element_text(color = "black", size = 8),
        axis.title.x = element_text(size = 11),
        strip.text.x = element_text(face = "bold", color = "black", size = 12),
        strip.background.x = element_blank(),
        strip.text.y = element_text(angle = 0, color = "black", size = 7),
        strip.background.y = element_rect(fill = "white"),
        panel.spacing.y = unit(0.2, "lines"))

# Save plot
ggsave("Output/Figures/phylo effect.png", plot = phylo_effect_plot, width = 20, height = 30, units = "cm", dpi = 600, type = "cairo")

#### ------------ END ------------ ####
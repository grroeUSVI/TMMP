
export_tables_figures <- function(output_dir) {
  

## Save tree measurements table
tree_measurements_table <- create_tree_measurement_table(tree_measurements = tree_measurements, tree_heights = tree_heights) %>% 
  gtsave(filename = paste(output_dir,"tree_measurements_table.pdf"))


## Site size frequency figures


## Save densiometer figure
densiometer_figure <- densiometer_chart(densiometer_data)
ggsave(filename =  paste(output_dir,"densiometer_figure.jpg"), plot = densiometer_figure, width = 7.5, height = 9.8, dpi = 300, units = "in", device = "jpg")

## seedling density figures
seedling_figure_nobreak <- seedling_density(regen = regen_data, densio = densiometer_data)
seedling_figure_break <- seedling_density(regen = regen_data, densio = densiometer_data, breaks = 15)
ggsave(filename =  paste(output_dir,"seedling_figure_nobreak.jpg"), plot = seedling_figure_nobreak, width = 9.8, height = 7.5, dpi = 300, units = "in", device = "jpg")
ggsave(filename =  paste(output_dir,"seedling_figure_break.jpg"), plot = seedling_figure_break, width = 9.8, height = 7.5, dpi = 300, units = "in", device = "jpg")

## seedling relative abundance table
seedling_density_table <- seedling_rel_abundance_table(regen_data = regen_data)%>%
  gtsave(filename =  paste(output_dir,"seedling_density_table.pdf"))

## Sapling density figures
rhma_density <- sapling_density(regen = regen_data, sapling = sapling_data, species = "RHMA", densio = densiometer_data)
avge_density <- sapling_density(regen = regen_data, sapling = sapling_data, species = "AVGE", densio = densiometer_data)
lara_density <- sapling_density(regen = regen_data, sapling = sapling_data, species = "LARA", densio = densiometer_data)
ggsave(filename =  paste(output_dir,"rhma_density.jpg"), plot = rhma_density, width = 9.8, height = 7.5, dpi = 300, units = "in", device = "jpg")
ggsave(filename =  paste(output_dir,"avge_density.jpg"), plot = avge_density, width = 9.8, height = 7.5, dpi = 300, units = "in", device = "jpg")
ggsave(filename =  paste(output_dir,"lara_density.jpg"), plot = lara_density, width = 9.8, height = 7.5, dpi = 300, units = "in", device = "jpg")

## Water quality table
water_quality_table <- water_qual_table(YSI = YSI_data) %>%
  gtsave(filename =  paste(output_dir,"water_quality_table.pdf"))

for (site in unique(tree_measurements$Site)) {
  
  p <- site_LF(tree_measurements = tree_measurements, site = site, bin_size = 3)
  ggsave(filename =  paste(output_dir, site,"_size_frequency.jpg"), plot = p, width = 9.8, height = 7.5, dpi = 300, units = "in", device = "jpg")
}

}
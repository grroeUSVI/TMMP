library(tidyverse)
library(gt)
library(readxl)


## Function parameters refer to the corresponding worksheet in the imported excel file


################################################################################################
## Table and Figure creation functions
################################################################################################


# ------------------------------------------------------------------------------
# Function: site_LF
# Purpose: Plot the size frequency distribution (SFD) of tree diameters (DBH)
#          by site and year, separated by mortality status (alive/dead).
#
# Parameters:
#   - tree_measurements: Data frame of individual tree records, including DBH,
#                        status, and year of measurement.
#   - site: String. Name of the site to filter the data by.
#   - bin_size: Numeric. Width of the DBH size classes (bins).
#
# Returns:
#   - A ggplot object showing stacked bar plots of DBH distributions per year,
#     faceted by year, with colors representing mortality status.
# ------------------------------------------------------------------------------

site_LF <- function(tree_measurements, site, bin_size) {
  
  a <- tree_measurements %>% filter(Site == site & DBH_cm >= 0 & Mortality != "NA") %>%  
    mutate(Mortality = if_else(Mortality == "Dying", "Alive", Mortality)) %>% 
    group_by(SY, DBH_cm, Mortality) %>% 
    summarise(n = n(), .groups = "drop") %>%
    group_by(SY) %>% 
    mutate(freq = n / sum(n)) %>% 
    ungroup()
  
  
  bins <- seq(0, max(a$DBH_cm + bin_size), by = bin_size)
  
  
  b <- a %>% 
    group_by(SY, Mortality) %>% 
    nest() %>% 
    mutate(LF = map(data, ~ .x %>%
                      data.frame() %>% 
                      mutate(Bin = cut(DBH_cm, breaks = bins, right = FALSE)) %>%
                      group_by(Bin) %>%
                      summarise(total_freq = sum(freq, na.rm = TRUE)))) %>% 
    unnest(LF) %>% 
    ungroup()
  
  
  br <- unique(b$Bin)
  la <- labeler(bin_num = length(unique(b$Bin)) , bin_size = bin_size)
  
  ggplot(b, aes(x=Bin, y=total_freq, fill=Mortality)) +
    geom_bar(stat="identity", position = "stack", width = .9, color="black", linewidth=.5) +
    scale_x_discrete(breaks = br,
                     labels = la[1:(length(unique(b$Bin)))],
                     limits = factor(br)) +
    scale_fill_manual(values = c("Alive" = "darkolivegreen3", "Dead" = "burlywood4")) +
    xlab(label = "DBH Size Class (cm)") + 
    ylab(label = "Relative Frequency") + 
    ggtitle(site) +
    theme_Publication() +
    facet_wrap(~SY)
  
}

# ------------------------------------------------------------------------------
# Function: all_sites_LF
# Purpose: Plot the size frequency distribution (SFD) of tree diameters (DBH)
#          across all sites for a given year, separated by mortality status.
#
# Parameters:
#   - tree_measurements: Data frame of individual tree records, including DBH,
#                        status, site, and year of measurement.
#   - year: Numeric or string. The sampling year (SY) to filter by.
#   - bin_size: Numeric. Width of the DBH size classes (bins).
#
# Returns:
#   - A ggplot object showing stacked bar plots of DBH distributions by site
#     for the given year. Each facet represents a site and shows mortality
#     status in separate color fills.
# ------------------------------------------------------------------------------

all_sites_LF <- function(tree_measurements, year, bin_size) {
  
  a <- tree_measurements %>% 
    filter(SY == year & DBH_cm >= 0 & Mortality != "NA") %>%  
    mutate(Mortality = if_else(Mortality == "Dying", "Alive", Mortality)) %>% 
    group_by(Site, DBH_cm, Mortality) %>% 
    summarise(n = n(), .groups = "drop") %>%
    group_by(Site) %>% 
    mutate(freq = n / sum(n)) %>% 
    ungroup()
  
  bins <- seq(0, max(a$DBH_cm + bin_size), by = bin_size)
  
  b <- a %>% 
    group_by(Site, Mortality) %>% 
    nest() %>% 
    mutate(LF = map(data, ~ .x %>%
                      data.frame() %>% 
                      mutate(Bin = cut(DBH_cm, breaks = bins, right = FALSE)) %>%
                      group_by(Bin) %>%
                      summarise(total_freq = sum(freq, na.rm = TRUE)))) %>% 
    unnest(LF) %>% 
    ungroup()
  
  br <- unique(b$Bin)
  la <- labeler(bin_num = length(unique(b$Bin)) , bin_size = bin_size)
  
  ggplot(b, aes(x = Bin, y = total_freq, fill = Mortality)) +
    geom_bar(stat = "identity", position = "stack", width = 0.9, color = "black", linewidth = 0.5) +
    scale_x_discrete(breaks = br,
                     labels = la[1:length(unique(b$Bin))],
                     limits = factor(br)) +
    scale_fill_manual(values = c("Alive" = "darkolivegreen3", "Dead" = "burlywood4")) +
    xlab("DBH Size Class (cm)") +
    ylab("Relative Frequency") +
    ggtitle(paste("Size Frequency Distribution -", year)) +
    theme_Publication() +
    facet_wrap(~ Site) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}



# ------------------------------------------------------------------------------
# Function: create_tree_measurement_table
# Purpose: Combine multiple forest structure metrics (e.g., basal area, height,
#          stem density, species contribution) into a single summary table
#          formatted for reporting.
#
# Parameters:
#   - tree_measurements: Data frame of tree records, including DBH and species.
#   - tree_heights: Data frame of tree height measurements by plot and year.
#
# Returns:
#   - A gt table object summarizing:
#       - Basal area (m²/ha)
#       - Mean tree height (m)
#       - Stem density (stems/ha)
#       - Top 3 species by basal area (% contribution)
#     Grouped by site and year with "mean ± SE" format for numeric values.
# ------------------------------------------------------------------------------

create_tree_measurement_table <- function(tree_measurements, tree_heights) {
  
  a <- mean_basal_area(tree_measurements)
  b <- spp_contibution(tree_measurements, mortality = "alive")
  c <- spp_contibution(tree_measurements, mortality = "dead")
  d <- mean_tree_height(tree_heights)
  e <- mean_stem_density(tree_measurements)
  
  x <- reduce(.x = list(a,b,c,d,e), .f = full_join)
  
  format_x <- x %>% 
    mutate(across(where(is.numeric), round, 1)) %>% 
    mutate(mean_stem_density    = round(mean_stem_density, digits = 0),
           mean_stem_density_SE = round(mean_stem_density_SE, digits = 0)) %>% 
    mutate(
      Height = paste(mean_height, "±", mean_height_SE),
      Basal_Area = paste(mean_basal_area, "±", mean_basal_area_SE),
      Density = paste(mean_stem_density, "±", mean_stem_density_SE)
    ) %>% 
    select(SY, Site, Height, Basal_Area, Density, RHMA_alive, AVGE_alive, LARA_alive, RHMA_dead, AVGE_dead, LARA_dead)
  
# Table using gt
 table <- format_x %>% 
    # gt(groupname_col = c("Site")) %>%
    gt(rowname_col = "SY", groupname_col = "Site") %>%
    # cols_add(empty = NA_character_, .before = "SY") %>%
    # sub_missing(columns = empty, missing_text = "     ") %>%
    fmt_markdown() %>% 
    cols_label(
      # empty = md('&emsp;&emsp;&emsp;'),
      SY = md("Year <br/> &emsp; "),
      RHMA_alive = md("*R. mangle* <br/> (%)"),
      AVGE_alive = md("*A. germinians* <br/> (%)"),
      LARA_alive = md("*L. racemosa* <br/> (%)"),
      RHMA_dead = md("*R. mangle* <br/> (%)"),
      AVGE_dead = md("*A. germinians* <br/> (%)"),
      LARA_dead = md("*L. racemosa* <br/> (%)"),
      Height = md("Height <br/> (m)"),
      Basal_Area = md("Basal area <br/> (m²/ha)"),
      Density = md("Density <br/> (stems/ha)")
    ) %>% 
    tab_spanner(label = "Relative Distribution of Species (Alive)", columns = RHMA_alive:LARA_alive) %>% 
    tab_spanner(label = "Relative Distribution of Species (Dead)", columns = RHMA_dead:LARA_dead) %>%
    cols_align(align = "center", everything()) %>%
    # cols_width(empty ~ px(30),
    #            Density ~ px(120),
    #            everything() ~ px(100)) %>%
   cols_width(everything() ~ px(100)) %>%
   tab_options(
     table.font.size = px(14),
     row_group.as_column = FALSE,
     data_row.padding = px(5)
   )
   # opt_interactive(
   # use_search = TRUE,
   # use_filters = FALSE,
   # use_resizers = TRUE,
   # use_highlight = TRUE,
   # use_compact_mode = FALSE,
   # use_text_wrapping = FALSE,
   # use_page_size_select = TRUE
   # )
 
 return(table)
  
}

# ------------------------------------------------------------------------------
# Function: densiometer_data_table
# Purpose: Summarize canopy structure components (vegetation, wood, open) as
#          percent cover using data from convex densiometer readings.
#
# Parameters:
#   - densiometer_data: Data frame of densiometer readings including plot, year,
#                       and type of cover (vegetation, wood, open).
#
# Returns:
#   - A gt table object showing mean ± SE of percent cover for each cover type,
#     grouped by site and year.
# ------------------------------------------------------------------------------

densiometer_data_table <- function(densiometer_data) {
  a <- densiometer_data %>% select(Site, Plot, SY, starts_with("Calc"))
  b <- df_long <- a %>%
    mutate(across(starts_with("Calc_"), ~ suppressWarnings(as.numeric(.)))) %>% 
    pivot_longer(cols = starts_with("Calc_"),
                 names_to = c("Direction", "Type"),
                 names_pattern = "Calc_([NSEW])_([A-Za-z]+)")
  
  c <- b %>% group_by(SY, Site, Plot, Type) %>% 
    summarise(perc_plot = mean(value, na.rm = T), .groups = "drop")
  d <- c %>% group_by(SY, Site, Type) %>% 
    summarise(perc = mean(perc_plot, na.rm = T), perc_SE = standard_error(perc_plot, na.rm = T), .groups = "drop")
  
  
  d2 <- d %>%
    pivot_wider(
      id_cols = c(SY, Site),
      names_from = Type,
      values_from = c(perc, perc_SE),
      values_fill = 0,
      names_glue = "{Type}_{.value}"
    ) %>% 
    select(SY, Site, Veg_perc, Veg_perc_SE, Wood_perc, Wood_perc_SE, Open_perc, Open_perc_SE)
  
  format_x <- d2 %>% 
    mutate(across(where(is.numeric), round, 1)) %>% 
    mutate(
      Alive = paste(Veg_perc, "±", Veg_perc_SE),
      Dead = paste(Wood_perc, "±", Wood_perc_SE),
      Open = paste(Open_perc, "±", Open_perc_SE)
    ) %>% 
    select(SY, Site, Alive, Dead, Open)
  
  # Table using gt
  table <- format_x %>% 
    gt(rowname_col = "SY", groupname_col = "Site") %>%
    cols_add(empty = NA_character_, .before = "SY") %>%
    sub_missing(columns = empty, missing_text = "     ") %>%
    fmt_markdown() %>% 
    cols_label(
      empty = md('&emsp;&emsp;&emsp;'),
      SY = md("Year <br/> &emsp; "),
      Alive = md("Alive"),
      Dead = md("Dead"),
      Open = md("Open")
    ) %>% 
    cols_align(align = "center", everything()) %>%
    cols_width(empty ~ px(30),
               everything() ~ px(100)) %>% 
    tab_options(
      table.font.size = px(14),  
      row_group.as_column = FALSE,  
      data_row.padding = px(5)
    )
  # opt_interactive(
  # use_search = TRUE,
  # use_filters = FALSE,
  # use_resizers = TRUE,
  # use_highlight = TRUE,
  # use_compact_mode = FALSE,
  # use_text_wrapping = FALSE,
  # use_page_size_select = TRUE
  # )

  return(table)
}

# ------------------------------------------------------------------------------
# Function: densiometer_figure
# Purpose: Visualize changes in canopy structure (vegetation, wood, and open
#          canopy cover) across years using a diverging bar and line chart.
#
# Parameters:
#   - densiometer_data: Data frame of densiometer observations including
#                       vegetation, wood, and open readings by plot and year.
#
# Returns:
#   - A ggplot object showing:
#       - Stacked diverging bars (vegetation: positive, wood: negative)
#       - Line and point overlay for open canopy %
#     Faceted by site for comparison over time.
# ------------------------------------------------------------------------------

densiometer_figure <- function(densiometer_data) {
  a <- densiometer_data %>% select(Site, Plot, SY, starts_with("Calc"))
  b <- df_long <- a %>%
    mutate(across(starts_with("Calc_"), ~ suppressWarnings(as.numeric(.)))) %>% 
    pivot_longer(cols = starts_with("Calc_"),
                 names_to = c("Direction", "Type"),
                 names_pattern = "Calc_([NSEW])_([A-Za-z]+)")
  
  c <- b %>% group_by(SY, Site, Plot, Type) %>% 
    summarise(perc_plot = mean(value, na.rm = T), .groups = "drop")
  d <- c %>% group_by(SY, Site, Type) %>% 
    summarise(perc = mean(perc_plot, na.rm = T), perc_SE = standard_error(perc_plot, na.rm = T), .groups = "drop")
  
  
  d2 <- d %>% 
    filter(Type %in% c("Veg", "Wood", "Open")) %>% 
    select(SY, Site, Type, perc) %>% 
    pivot_wider(names_from = Type, values_from = perc, values_fill = 0)
  
  d2 %>%
    ggplot(aes(x = factor(SY))) +
    geom_bar(aes(y = Veg, fill = "Alive"), stat = "identity") +
    geom_bar(aes(y = -Wood, fill = "Dead"), stat = "identity") +
    geom_line(aes(y = Open, group = 1, color = "Open Area"), 
              size = 1.5, alpha = .85, linetype = "solid") +
    geom_point(aes(y = Open, group = 1, color = "Open Area"), 
               size = 1.5, alpha = .50) +
    facet_wrap(~Site) +
    scale_y_continuous(limits = c(-50,100), labels = abs) +
    labs(y = "Canopy Cover (%)", x = "Year", fill = "Cover Type", color = element_blank()) +
    theme_Publication() +
    scale_fill_manual(values = c("Alive" = "darkolivegreen3", "Dead" = "burlywood4")) +
    scale_color_manual(values = c("Open Area" = "deepskyblue"))
}

# ------------------------------------------------------------------------------
# Function: seedling_density
# Purpose: Create a bar plot of seedling stem density by site and year, with
#          error bars, and overlay average seedling height using a secondary axis.
#
# Parameters:
#   - regen: Data frame of regeneration data (seedlings/saplings) with species,
#            plot, year, and height class.
#   - densio: Data frame used to extract island-site relationships.
#   - breaks: Optional. Custom breaks for the y-axis scale.
#
# Returns:
#   - A ggplot object displaying:
#       - Seedling stem density per hectare (bars with SE)
#       - Mean seedling height (asterisk markers on a secondary y-axis)
#     Grouped by site and year, and faceted by island.
# ------------------------------------------------------------------------------

seedling_density <- function(regen, densio, breaks = NULL) {
  
  site <- densio %>% 
    group_by(Island) %>% 
    reframe(Site = unique(Site))
  
  a <- regen %>%
    pivot_longer(cols = c(RHMA_seedlings, RHMA_saplings, 
                          LARA_seedlings, LARA_saplings, 
                          AVGE_seedlings, AVGE_saplings),
                 names_to = c("Species", "Stage"),
                 names_pattern = "([A-Z]+)_(seedlings|saplings)") %>%
    select(SY, Site, Plot, Species, Stage, Count = value) %>%
    filter(Stage == "seedlings") %>% 
    group_by(SY, Site, Plot, Stage) %>% 
    summarise(mean_plot = mean(Count), .groups = "drop") %>% 
    arrange(SY, Site, Plot, Stage) %>% 
    group_by(SY, Site, Stage) %>% 
    summarise(mean_site = mean(mean_plot), count_SE = standard_error(mean_plot), .groups = "drop")
  
  b <- regen %>% 
    select(SY, Site, Plot, Tall_seedling_cm) %>% 
    group_by(SY, Site, Plot) %>% 
    summarise(height_plot = mean(Tall_seedling_cm, na.rm = T ), .groups = "drop") %>% 
    group_by(SY, Site) %>% 
    summarise(height_site = mean(height_plot, na.rm = T), .groups = "drop")
  
  c <- a %>% 
    full_join(b) %>% 
    left_join(site) %>% 
    mutate(Island = factor(Island, levels = c("St Thomas", "St John", "St Croix")),
           SY = as.factor(SY))
  
  
  
  height_rescale <- function(y) {
    min_h <- min(c$height_site, na.rm = TRUE)
    max_h <- max(c$height_site, na.rm = TRUE)
    min_m <- min(c$mean_site, na.rm = TRUE)
    max_m <- max(c$mean_site, na.rm = TRUE) + max(c$count_SE, na.rm = T)
    
    # Scale height_site to match the primary y-axis range
    (y - min_h) / (max_h - min_h) * (max_m - min_m) + min_m
  }
  
  ggplot(c, aes(x = Site, y = mean_site, fill = SY)) + 
    geom_bar(stat = "identity", position = "dodge", width = 0.75) +
    geom_errorbar(aes(ymax = mean_site + count_SE, ymin = mean_site, color = SY), 
                  position = position_dodge(width = 0.75), width = 0.25, show.legend = FALSE) +
    guides(fill = guide_legend(title = "Year")) +
    geom_point(aes(y = height_rescale(height_site), color = SY), 
               shape = 8,  # Asterisk shape
               size = 4,
               position = position_dodge(width = 0.75),
               show.legend = FALSE) +
    scale_y_continuous(
      name = expression("Seedlings/m"^2), 
      sec.axis = sec_axis(~ (.- min(c$mean_site, na.rm = TRUE)) / 
                            ((max(c$mean_site, na.rm = TRUE) + max(c$count_SE, na.rm = T)) - min(c$mean_site, na.rm = TRUE)) *
                            (max(c$height_site, na.rm = TRUE) - min(c$height_site, na.rm = TRUE)) + 
                            min(c$height_site, na.rm = TRUE),
                            name = "Mean seedling height(cm)") 
    ) + {
      # if breaks parameter is not null add this to ggplot
      if (!is.null(breaks)) scale_y_break(c(breaks, breaks + 1), space = .025, scales = "free")
    } +
    labs(x = element_blank(), fill = "SY", color = "SY") +
    theme_Publication() +
    scale_color_manual(values = c("#CC3333", "#006600", "#3366CC")) +
    theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1),
          axis.text = element_text(size = 12),
          axis.title.y.right = element_text(family = "Calibri", size = 12),
          axis.text.y.right = element_text(family = "Calibri", size = 12),
          strip.placement='outside',
          strip.background.x=element_blank(),
          strip.text=element_text(size=12,color="black",face="bold"),
          panel.spacing.x=unit(0,"pt")) +
    facet_grid(cols=vars(Island),scales="free_x",space="free_x",switch="x")
}

# ------------------------------------------------------------------------------
# Function: sapling_density
# Purpose: Plot sapling stem density by species, site, and year with error bars,
#          and overlay mean sapling height on a secondary y-axis.
#
# Parameters:
#   - regen: Data frame of regeneration data containing sapling counts.
#   - sapling: Data frame with sapling height measurements.
#   - species: String. Target species to filter for (e.g., "PIMA").
#   - densio: Data frame used to extract island-site mapping.
#   - breaks: Optional. Custom breaks for the y-axis scale.
#
# Returns:
#   - A ggplot object showing:
#       - Sapling density (bars with SE) by species and site
#       - Average sapling height (asterisks) overlaid on a secondary axis
#     Faceted by island for geographic comparison.
# ------------------------------------------------------------------------------

sapling_density <- function(regen, sapling, species, densio, breaks = NULL) {

  pretty_name <- function(species) {
    
    a <- case_when(
      species == "RHMA" ~ "R. mangle",
      species == "AVGE" ~ "A. germinians",
      species == "LARA" ~ "L. racemosa"
    )  
  return(a)
  }
  
  site <- densio %>% 
    group_by(Island) %>% 
    reframe(Site = unique(Site))
  
  a <- regen %>%
    pivot_longer(cols = c(RHMA_seedlings, RHMA_saplings,
                          LARA_seedlings, LARA_saplings,
                          AVGE_seedlings, AVGE_saplings),
                 names_to = c("Species", "Stage"),
                 names_pattern = "([A-Z]+)_(seedlings|saplings)") %>%
    select(SY, Site, Plot, Quadrat, Species, Stage, Count = value) %>%
    filter(Stage == "saplings") %>%
    group_by(SY, Site, Plot, Quadrat, Species) %>%
    summarise(regen_sum = sum(Count, na.rm = T), .groups = "drop") %>% 
    group_by(SY, Site, Plot, Species) %>% 
    summarise(plot_dens = mean(regen_sum, na.rm=TRUE), plot_dens_SE = standard_error(regen_sum, na.rm=TRUE), .groups = "drop")
  
  b <- sapling %>%
    select(SY, Site, Plot, Quadrat, Species, Height_cm) %>%
    group_by(SY, Site, Plot, Quadrat, Species) %>%
    summarise(quad_ht = mean(Height_cm, na.rm=TRUE), .groups = "drop") %>%
    group_by(SY, Site, Plot, Species) %>%
    summarise(plot_ht = mean(quad_ht, na.rm=TRUE), .groups = "drop")
  
  c <- a %>% left_join(b) %>%
    group_by(SY, Site, Species) %>% 
    summarise(mean_dens = mean(plot_dens, na.rm=TRUE), mean_dens_SE = standard_error(plot_dens), mean_ht = mean(plot_ht, na.rm=TRUE), .groups = "drop") %>% 
    left_join(site) %>% 
    mutate(Island = factor(Island, levels = c("St Thomas", "St John", "St Croix")),
           SY = as.factor(SY))
  
  
  height_rescale <- function(y) {
    min_h <- min(c$mean_ht, na.rm = TRUE)
    max_h <- max(c$mean_ht, na.rm = TRUE)
    min_m <- min(c$mean_dens, na.rm = TRUE)
    max_m <- max(c$mean_dens, na.rm = TRUE) + max(c$mean_dens_SE, na.rm = T)
    
    # Scale height_site to match the primary y-axis range
    (y - min_h) / (max_h - min_h) * (max_m - min_m) + min_m
  }
  
  c %>% 
    filter(Species == species) %>% 
    ggplot(aes(x = Site, y = mean_dens, fill = SY)) + 
    geom_bar(stat = "identity", position = "dodge", width = 0.75) +
    geom_errorbar(aes(ymax = mean_dens + mean_dens_SE, ymin = mean_dens, color = SY), 
                  position = position_dodge(width = 0.75), width = 0.25, show.legend = FALSE) +
    guides(fill = guide_legend(title = "Year")) +
    geom_point(aes(y = height_rescale(mean_ht), color = SY), 
               shape = 8,  # Asterisk shape
               size = 4,
               position = position_dodge(width = 0.75),
               show.legend = FALSE) +
    scale_y_continuous(
      name = expression("Saplings/m"^2), 
      sec.axis = sec_axis(~ (.- min(c$mean_dens, na.rm = TRUE)) / 
                            ((max(c$mean_dens, na.rm = TRUE) + max(c$mean_dens_SE, na.rm = T)) - min(c$mean_dens, na.rm = TRUE)) *
                            (max(c$mean_ht, na.rm = TRUE) - min(c$mean_ht, na.rm = TRUE)) + 
                            min(c$mean_ht, na.rm = TRUE),
                          name = "Mean sapling height(cm)")
    ) + {
      # if breaks parameter is not null add this to ggplot
      if (!is.null(breaks)) scale_y_break(c(breaks, breaks), space = .025, scales = "free")
    } +
    labs(x = element_blank(), fill = "SY", color = "SY") +
    theme_Publication() +
    scale_color_manual(values = c("#CC3333", "#006600", "#3366CC")) +
    theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1),
          axis.text = element_text(size = 12),
          axis.title.y.right = element_text(family = "Calibri", size = 12),
          axis.text.y.right = element_text(family = "Calibri", size = 12),
          strip.placement='outside',
          strip.background.x=element_blank(),
          strip.text=element_text(size=12,color="black",face="bold"),
          panel.spacing.x=unit(0,"pt")) +
    facet_grid(cols=vars(Island),scales="free_x",space="free_x",switch="x") +
    ggtitle(pretty_name(species)) +
    theme(plot.title = element_text(face = "italic"))

}

# ------------------------------------------------------------------------------
# Function: sapling_rel_abundance_table
# Purpose: Generate a formatted table of sapling stem density by species, site,
#          and year. Displays mean sapling density with associated standard error
#          for each species at each site-year combination.
#
# Parameters:
#   - regen: Data frame containing regeneration counts (saplings and seedlings).
#   - sapling: Unused in current function but reserved for compatibility/extension.
#   - species: String. Species code (e.g., "RHMA", "AVGE", "LARA").
#   - densio: Data frame containing site-to-island mapping.
#
# Returns:
#   - A gt table object with site- and year-grouped rows showing sapling density
#     estimates (mean ± SE) for each species, formatted for clear presentation.
# ------------------------------------------------------------------------------

sapling_rel_abundance_table <- function(regen, sapling, species, densio) {
  
  pretty_name <- function(species) {
    
    a <- case_when(
      species == "RHMA" ~ "R. mangle",
      species == "AVGE" ~ "A. germinians",
      species == "LARA" ~ "L. racemosa"
    )  
    return(a)
  }
  
  site <- densio %>% 
    group_by(Island) %>% 
    reframe(Site = unique(Site))
  
  a <- regen %>%
    pivot_longer(cols = c(RHMA_seedlings, RHMA_saplings,
                          LARA_seedlings, LARA_saplings,
                          AVGE_seedlings, AVGE_saplings),
                 names_to = c("Species", "Stage"),
                 names_pattern = "([A-Z]+)_(seedlings|saplings)") %>%
    select(SY, Site, Plot, Quadrat, Species, Stage, Count = value) %>%
    filter(Stage == "saplings") %>%
    group_by(SY, Site, Plot, Quadrat, Species) %>%
    summarise(regen_sum = sum(Count, na.rm = T), .groups = "drop") %>% 
    group_by(SY, Site, Plot, Species) %>% 
    summarise(plot_dens = mean(regen_sum, na.rm=TRUE), plot_dens_SE = standard_error(regen_sum, na.rm=TRUE), .groups = "drop")
  
  c <- a %>%
    group_by(SY, Site, Species) %>% 
    summarise(mean_dens = mean(plot_dens, na.rm=TRUE), mean_dens_SE = standard_error(plot_dens), .groups = "drop") %>% 
    pivot_wider(
      names_from = Species,
      values_from = c(mean_dens, mean_dens_SE),
      names_glue = "{Species}_{.value}"
    ) %>%
    rename_with(~ gsub("site_", "", .), starts_with(c("AVGE", "LARA", "RHMA"))) %>% 
    arrange(Site, SY)
  
  format <- c %>%
  mutate(across(where(is.numeric), round, 1)) %>%
  mutate(
    AVGE = paste(AVGE_mean_dens, "±", AVGE_mean_dens_SE),
    LARA = paste(LARA_mean_dens, "±", LARA_mean_dens_SE),
    RHMA = paste(RHMA_mean_dens, "±", RHMA_mean_dens_SE)
  ) %>%
  select(SY, Site, RHMA, AVGE, LARA)
  
  # Table using gt
  table <- format %>% 
    gt(rowname_col = "SY", groupname_col = "Site") %>% 
    cols_add(empty = NA_character_, .before = "RHMA") %>%
    sub_missing(columns = empty, missing_text = "     ") %>%
    fmt_markdown() %>% 
    cols_label(
      empty = md('&emsp;&emsp;&emsp;'),
      SY = "",
      RHMA = md("*R. mangle*"),
      AVGE = md("*A. germinians*"),
      LARA = md("*L. racemosa*")
    ) %>% 
    cols_align(align = "center", everything()) %>%
    cols_width(empty ~ px(30),
               everything() ~ px(100)) %>% 
    tab_options(
      table.font.size = px(14),  
      row_group.as_column = FALSE,  
      data_row.padding = px(5)
    )
  
  return(table)
  
}

# ------------------------------------------------------------------------------
# Function: seedling_rel_abundance_table
# Purpose: Create a summary table of seedling density (mean ± SE) by species, 
#          site, and year, based on regeneration data. Only seedling stage is used.
#          The table is formatted using the 'gt' package for presentation.
#
# Parameters:
#   - regen_data: Data frame containing regeneration survey data with seedling and sapling counts.
#
# Returns:
#   - A formatted 'gt' table showing mean seedling density ± standard error for
#     each species (*R. mangle*, *A. germinans*, *L. racemosa*) by site and year.
# ------------------------------------------------------------------------------
seedling_rel_abundance_table <- function(regen_data) {
  
  a <- regen_data %>% 
    pivot_longer(cols = c(RHMA_seedlings, RHMA_saplings, 
                          LARA_seedlings, LARA_saplings, 
                          AVGE_seedlings, AVGE_saplings),
                 names_to = c("Species", "Stage"),
                 names_pattern = "([A-Z]+)_(seedlings|saplings)") %>%
    select(SY, Site, Plot, Species, Stage, Count = value) %>%
    filter(Stage == "seedlings") %>% 
    group_by(SY, Site, Plot, Species) %>% 
    summarise(plot_mean = mean(Count, na.rm = T), .groups = "drop") %>% 
    group_by(SY, Site, Species) %>% 
    summarise(site_mean = mean(plot_mean, na.rm = T), site_SE = standard_error(plot_mean, na.rm = T), .groups = "drop") %>% 
    pivot_wider(
      names_from = Species,
      values_from = c(site_mean, site_SE),
      names_glue = "{Species}_{.value}"
    ) %>%
    rename_with(~ gsub("site_", "", .), starts_with(c("AVGE", "LARA", "RHMA"))) %>% 
    arrange(Site, SY)
  
  format <- a %>% 
    mutate(across(where(is.numeric), round, 1)) %>% 
    mutate(
      AVGE = paste(AVGE_mean, "±", AVGE_SE),
      LARA = paste(LARA_mean, "±", LARA_SE),
      RHMA = paste(RHMA_mean, "±", RHMA_SE)
    ) %>% 
    select(SY, Site, RHMA, AVGE, LARA)
  
  # Table using gt
  table <- format %>% 
    gt(rowname_col = "SY", groupname_col = "Site") %>% 
    cols_add(empty = NA_character_, .before = "RHMA") %>%
    sub_missing(columns = empty, missing_text = "     ") %>%
    fmt_markdown() %>% 
    cols_label(
      empty = md('&emsp;&emsp;&emsp;'),
      SY = "",
      RHMA = md("*R. mangle*"),
      AVGE = md("*A. germinians*"),
      LARA = md("*L. racemosa*")
    ) %>% 
    cols_align(align = "center", everything()) %>%
    cols_width(empty ~ px(30),
               everything() ~ px(100)) %>% 
    tab_options(
      table.font.size = px(14),  
      row_group.as_column = FALSE,  
      data_row.padding = px(5)
    )
  
  return(table)
  
}

# ------------------------------------------------------------------------------
# Function: water_qual_table
# Purpose: Generate a formatted table summarizing key water quality parameters 
#          (e.g., depth, temperature, dissolved oxygen, salinity, total dissolved solids) 
#          by site and year.
#
# Parameters:
#   - YSI: Data frame of water quality measurements collected in the field.
#
# Returns:
#   - A 'gt' table displaying mean ± standard error for each metric, grouped by site 
#     and year. Values are rounded and labeled with appropriate units for clarity.
# ------------------------------------------------------------------------------

water_qual_table <- function(YSI) {
  
  a <- mean_water_quality(YSI)
  
  format_x <- a %>% 
    mutate(across(where(is.numeric), round, 1)) %>% 
    mutate(mean_SAL = round(mean_SAL, digits = 0),
           mean_SAL_SE = round(mean_SAL_SE, digits = 0),
           mean_TDS = round(mean_TDS, digits = 0),
           mean_TDS_SE = round(mean_TDS_SE, digits = 0)) %>%
    mutate(
      Depth = paste(mean_depth, "±", mean_depth_SE),
      Temp = paste(mean_temp, "±", mean_temp_SE),
      DO = paste(mean_DO, "±", mean_DO_SE),
      SAL = paste(mean_SAL, "±", mean_SAL_SE),
      TDS = paste(mean_TDS, "±", mean_TDS_SE)
    ) %>% 
    select(SY, Site, Depth, Temp, DO, SAL, TDS)
  
  # Table using gt
  table <- format_x %>% 
    gt(rowname_col = "SY", groupname_col = "Site") %>%
    # cols_add(empty = NA_character_, .before = "SY") %>%
    # sub_missing(columns = empty, missing_text = "     ") %>%
    fmt_markdown() %>% 
    cols_label(
      # empty = md('&emsp;&emsp;&emsp;'),
      SY = md("Year <br/> &emsp; "),
      Depth = md("Depth <br/> (cm)"),
      Temp = md("Temperature <br/> (&deg;C)"),
      DO = md("Dissolved Oxygen <br/> (mg/L)"),
      SAL = md("Salinity <br/> (ppt)"),
      TDS = md("Total Dissolved Solids <br/> (mg/L)")
    ) %>% 
    cols_align(align = "center", everything()) %>%
    # cols_width(empty ~ px(30),
    #            Density ~ px(120),
    #            everything() ~ px(100)) %>%
    cols_width(everything() ~ px(120)) %>%
    tab_options(
      table.font.size = px(14),
      row_group.as_column = FALSE,
      data_row.padding = px(5)
    )
  # opt_interactive(
  # use_search = TRUE,
  # use_filters = FALSE,
  # use_resizers = TRUE,
  # use_highlight = TRUE,
  # use_compact_mode = FALSE,
  # use_text_wrapping = FALSE,
  # use_page_size_select = TRUE
  # )
  
  return(table)
  
}

################################################################################
## Dataframe Functions
################################################################################

# ------------------------------------------------------------------------------
# Function: site_coordinates
# Description:
#   Calculates the average latitude and longitude for each site by averaging the 
#   coordinates of its three plots.
#
# Parameters:
#   Coordinates - Dataframe containing plot-level coordinates and metadata.
#
# Returns:
#   Dataframe with site-level coordinates and forest type.
# ------------------------------------------------------------------------------
site_coordinates <- function(Coordinates) {
  a <- Coordinates %>% 
    group_by(Island, Site) %>% 
    summarise(
      Typology = first(`Forest type`),
      Latitude = mean(Latitude),
      Longitude = mean(Longitude),
      .groups = "drop"
    )
  return(a)
}

# ------------------------------------------------------------------------------
# Function: mean_basal_area
# Description:
#   Calculates the mean and standard error of basal area by site and year 
#   for living trees (Alive or Dying).
#
# Parameters:
#   tree_measurements - Dataframe with tree DBH and metadata.
#
# Returns:
#   Dataframe with mean and SE of basal area per site and year.
# ------------------------------------------------------------------------------
mean_basal_area <- function(tree_measurements) {
  a <- tree_measurements %>% 
    filter(Species %in% c("RHMA", "AVGE", "LARA") & Mortality %in% c("Alive", "Dying")) %>% 
    select(SY, Site, Plot, Species, DBH_cm) %>% 
    mutate(basal_area = basal_area(DBH_cm)) %>% 
    group_by(SY, Site, Plot) %>% 
    summarise(basal_plot = sum(basal_area, na.rm = TRUE), .groups = "drop") %>% 
    group_by(SY, Site) %>% 
    summarise(
      mean_basal_area = mean(basal_plot, na.rm = TRUE),
      mean_basal_area_SE = standard_error(basal_plot),
      .groups = "drop"
    ) %>% 
    mutate(across(everything(), ~replace_na(., 0)))
  return(a)
}

# ------------------------------------------------------------------------------
# Function: mean_basal_area_dead
# Description:
#   Calculates the mean and SE of basal area by site and year for dead trees only.
#
# Parameters:
#   tree_measurements - Dataframe with tree DBH and metadata.
#
# Returns:
#   Dataframe with mean and SE of basal area for dead trees per site and year.
# ------------------------------------------------------------------------------
mean_basal_area_dead <- function(tree_measurements) {
  a <- tree_measurements %>% 
    filter(Species %in% c("RHMA", "AVGE", "LARA") & Mortality == "Dead") %>% 
    select(SY, Site, Plot, Species, DBH_cm) %>% 
    mutate(basal_area = basal_area(DBH_cm)) %>% 
    group_by(SY, Site, Plot) %>% 
    summarise(basal_plot = sum(basal_area, na.rm = TRUE), .groups = "drop") %>% 
    group_by(SY, Site) %>% 
    summarise(
      mean_basal_area = mean(basal_plot, na.rm = TRUE),
      mean_basal_area_SE = standard_error(basal_plot),
      .groups = "drop"
    ) %>% 
    filter(!is.na(SY)) %>% 
    mutate(across(everything(), ~replace_na(., 0)))
  return(a)
}

# ------------------------------------------------------------------------------
# Function: spp_contibution
# Description:
#   Calculates percent contribution of each species to total basal area per 
#   site and year. Supports either alive or dead trees.
#
# Parameters:
#   tree_measurements - Dataframe with tree DBH and metadata.
#   mortality         - String ("alive" or "dead") to filter tree status.
#
# Returns:
#   Dataframe with percent basal area contribution per species, per site and year.
# ------------------------------------------------------------------------------
spp_contibution <- function(tree_measurements, mortality) {
  mortality_filter <- if (tolower(mortality) == "alive") {
    c("Alive", "Dying")
  } else {
    c("Dead")
  }
  suffix <- if (mortality == "alive") "_alive" else "_dead"
  
  x <- tree_measurements %>%
    filter(Species %in% c("RHMA", "AVGE", "LARA") & Mortality %in% mortality_filter) %>%
    select(SY, Site, Plot, Species, DBH_cm) %>%
    group_by(SY, Site, Plot, Species) %>%
    summarise(basal_area_sum = sum(basal_area(DBH_cm), na.rm = TRUE), .groups = "drop") %>% 
    pivot_wider(names_from = Species, values_from = basal_area_sum, values_fill = 0) %>%
    pivot_longer(cols = LARA:AVGE, names_to = "Species", values_to = "basal_area_sum") %>% 
    group_by(SY, Site, Species) %>% 
    summarise(mean_basal_area = mean(basal_area_sum), .groups = "drop") %>%
    group_by(SY, Site) %>% 
    mutate(contibution = (mean_basal_area / sum(mean_basal_area)) * 100) %>% 
    ungroup() %>% 
    pivot_wider(id_cols = !mean_basal_area, names_from = Species, values_from = contibution) %>% 
    rename_with(.cols = c(RHMA, AVGE, LARA), ~ paste0(., suffix))
  return(x)
}

# ------------------------------------------------------------------------------
# Function: mean_stem_density
# Description:
#   Calculates mean and SE stem density (stems/ha) by site and year for all 
#   mangrove species (alive or dying).
#
# Parameters:
#   tree_measurements - Dataframe with tree DBH and metadata.
#
# Returns:
#   Dataframe with mean and SE of stem density.
# ------------------------------------------------------------------------------
mean_stem_density <- function(tree_measurements) {
  a <- tree_measurements %>% 
    filter(Species %in% c("RHMA", "AVGE", "LARA") & Mortality %in% c("Alive", "Dying")) %>% 
    select(SY, Site, Plot, Species, DBH_cm) %>% 
    group_by(SY, Site, Plot) %>% 
    summarise(stem_plot_density = length(DBH_cm) / 0.01, .groups = "drop") %>% 
    group_by(SY, Site) %>% 
    summarise(
      mean_stem_density = mean(stem_plot_density, na.rm = TRUE),
      mean_stem_density_SE = standard_error(stem_plot_density),
      .groups = "drop"
    )
  return(a)
}

# ------------------------------------------------------------------------------
# Function: mean_tree_height
# Description:
#   Calculates mean and SE of tree height per site and year for mangroves 
#   (including NA and UNK).
#
# Parameters:
#   tree_heights - Dataframe with tree height and species info.
#
# Returns:
#   Dataframe with mean and SE of tree height.
# ------------------------------------------------------------------------------
mean_tree_height <- function(tree_heights) {
  a <- tree_heights %>% 
    filter(Species %in% c("RHMA", "AVGE", "LARA", NA, "UNK")) %>% 
    select(SY, Site, Plot, Species, Tree_Height_m) %>% 
    group_by(SY, Site, Plot) %>% 
    summarise(mean_height_plot = mean(Tree_Height_m, na.rm = TRUE), .groups = "drop") %>% 
    group_by(SY, Site) %>% 
    summarise(
      mean_height = mean(mean_height_plot, na.rm = TRUE),
      mean_height_SE = standard_error(mean_height_plot),
      .groups = "drop"
    ) %>% 
    mutate(across(everything(), ~replace_na(., 0)))
  return(a)
}

# ------------------------------------------------------------------------------
# Function: percent_basal_area_change
# Description:
#   Calculates the percentage change in mean basal area between 2022 and 2024.
#
# Parameters:
#   tree_measurements - Dataframe with tree DBH and metadata.
#
# Returns:
#   Dataframe with percent change in basal area per site.
# ------------------------------------------------------------------------------
percent_basal_area_change <- function(tree_measurements) {
  a <- mean_basal_area(tree_measurements) %>% 
    select(SY, Site, mean_basal_area) %>%
    pivot_wider(names_from = SY, values_from = mean_basal_area, values_fill = 0) %>%
    mutate(percentage_change = ((`2024` - `2022`) / `2022`) * 100)
  return(a)
}

# ------------------------------------------------------------------------------
# Function: total_basal_area_change
# Description:
#   Calculates absolute change in mean basal area (m²) from 2022 to 2024.
#
# Parameters:
#   tree_measurements - Dataframe with tree DBH and metadata.
#
# Returns:
#   Dataframe with absolute change in basal area.
# ------------------------------------------------------------------------------
total_basal_area_change <- function(tree_measurements) {
  a <- mean_basal_area(tree_measurements) %>% 
    select(SY, Site, mean_basal_area) %>%
    pivot_wider(names_from = SY, values_from = mean_basal_area, values_fill = 0) %>%
    mutate(total_change = (`2024` - `2022`)) %>% 
    select(Site, total_change) %>% 
    mutate(across(everything(), ~replace_na(., 0)))
  return(a)
}

# ------------------------------------------------------------------------------
# Function: total_basal_area_change_dead
# Description:
#   Calculates absolute change in mean basal area (m²) for dead trees from 
#   2022 to 2024.
#
# Parameters:
#   tree_measurements - Dataframe with tree DBH and metadata.
#
# Returns:
#   Dataframe with absolute change in basal area for dead trees.
# ------------------------------------------------------------------------------
total_basal_area_change_dead <- function(tree_measurements) {
  a <- mean_basal_area_dead(tree_measurements) %>% 
    select(SY, Site, mean_basal_area) %>%
    pivot_wider(names_from = SY, values_from = mean_basal_area, values_fill = 0) %>%
    mutate(total_change_dead = (`2024` - `2022`)) %>% 
    select(Site, total_change_dead) %>% 
    mutate(across(everything(), ~replace_na(., 0)))
  return(a)
}

# ------------------------------------------------------------------------------
# Function: total_tree_height_change
# Description:
#   Calculates change in average tree height (m) from 2022 to 2024.
#
# Parameters:
#   tree_heights - Dataframe with tree height and species info.
#
# Returns:
#   Dataframe with change in mean tree height.
# ------------------------------------------------------------------------------
total_tree_height_change <- function(tree_heights) {
  a <- mean_tree_height(tree_heights) %>% 
    select(SY, Site, mean_height) %>%
    pivot_wider(names_from = SY, values_from = mean_height, values_fill = 0) %>%
    mutate(height_change = (`2024` - `2022`)) %>% 
    select(Site, height_change) %>% 
    mutate(across(everything(), ~replace_na(., 0)))
  return(a)
}

# ------------------------------------------------------------------------------
# Function: total_stem_density_change
# Description:
#   Calculates change in stem density from 2022 to 2024.
#
# Parameters:
#   tree_measurements - Dataframe with tree DBH and metadata.
#
# Returns:
#   Dataframe with stem density change per site.
# ------------------------------------------------------------------------------
total_stem_density_change <- function(tree_measurements) {
  a <- mean_stem_density(tree_measurements) %>% 
    select(SY, Site, mean_stem_density) %>%
    pivot_wider(names_from = SY, values_from = mean_stem_density, values_fill = 0) %>%
    mutate(stem_density_change = (`2024` - `2022`)) %>% 
    select(Site, stem_density_change) %>% 
    mutate(across(everything(), ~replace_na(., 0)))
  return(a)
}

# ------------------------------------------------------------------------------
# Function: mean_water_quality
# Description:
#   Calculates mean and SE for various water quality parameters (e.g., DO, 
#   salinity, temperature) by site and year.
#
# Parameters:
#   YSI - Dataframe containing YSI water quality measurements.
#
# Returns:
#   Dataframe with mean and SE of each water quality variable.
# ------------------------------------------------------------------------------
mean_water_quality <- function(YSI) {
  a <- YSI %>% 
    select(SY, Site, Plot, Water_depth, Temp, DO_mg_L, Salinity_ppt, TDS) %>% 
    group_by(SY, Site, Plot) %>% 
    summarise(
      plot_depth = mean(Water_depth, na.rm = TRUE),
      plot_temp = mean(Temp, na.rm = TRUE),
      plot_DO = mean(DO_mg_L, na.rm = TRUE),
      plot_SAL = mean(Salinity_ppt, na.rm = TRUE),
      plot_TDS = mean(TDS, na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    group_by(SY, Site) %>% 
    summarise(
      mean_depth = mean(plot_depth, na.rm = TRUE), mean_depth_SE = standard_error(plot_depth),
      mean_temp = mean(plot_temp, na.rm = TRUE), mean_temp_SE = standard_error(plot_temp),
      mean_DO = mean(plot_DO, na.rm = TRUE), mean_DO_SE = standard_error(plot_DO),
      mean_SAL = mean(plot_SAL, na.rm = TRUE), mean_SAL_SE = standard_error(plot_SAL),
      mean_TDS = mean(plot_TDS, na.rm = TRUE), mean_TDS_SE = standard_error(plot_TDS),
      .groups = "drop"
    )
  return(a)
}


################################################################################
## Helper Functions
################################################################################

# ------------------------------------------------------------------------------
# Function: basal_area
# Description:
#   Calculates the basal area of a tree from its diameter at breast height (DBH).
#   The output is in square meters per hectare (m²/ha), assuming a fixed plot area.
#
# Parameters:
#   DBH - Numeric vector. Diameter at breast height in centimeters (cm).
#
# Returns:
#   Numeric vector. Basal area in square meters per hectare (m²/ha).
# ------------------------------------------------------------------------------
basal_area <- function(DBH) {
  x <- (pi * ((DBH / 100) / 2)^2) / 0.01
  return(x)
}


# ------------------------------------------------------------------------------
# Function: standard_error
# Description:
#   Calculates the standard error of a numeric vector.
#   Standard error is computed as the standard deviation divided by the 
#   square root of the sample size.
#
# Parameters:
#   x     - Numeric vector of values.
#   na.rm - Logical. If TRUE (default), missing values are removed before 
#           calculation.
#
# Returns:
#   Numeric. The standard error of the input vector.
# ------------------------------------------------------------------------------
standard_error <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- na.omit(x)
  }
  sd_x <- sd(x)
  n <- length(x)
  se_x <- sd_x / sqrt(n)
  return(se_x)
}

library(dplyr)
library(readxl)
library(tidyr)
library(purrr)
library(gt)
library(ggbreak)
library(ggplot2)

# calculate basal area 
# DBH: field measurement (cm)
# returns basal area in m^2/ha
basal_area <- function(DBH) {
  
  x <- (pi*((DBH/100)/2)^2) / .01
  return(x)
}

# Function to calculate standard error
standard_error <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- na.omit(x)
  }
  sd_x <- sd(x)
  n <- length(x)
  se_x <- sd_x / sqrt(n)
  return(se_x)
}

# Calculate mean +- SE basal area by year and site
# Hard coded filters: red, white, black mangroves & Alive or Dying
mean_basal_area <- function(tree_measurements) {
  
  a <- tree_measurements %>% 
    filter(Species %in% c("RHMA", "AVGE", "LARA") & Mortality %in% c("Alive", "Dying")) %>% 
    select(SY, Site, Plot, Species, DBH_cm) %>% 
    mutate(basal_area = basal_area(DBH_cm)) %>% 
    group_by(SY, Site, Plot) %>% 
    summarise(basal_plot = sum(basal_area, na.rm = T), .groups = "drop") %>% 
    group_by(SY, Site) %>% 
    summarise(mean_basal_area = mean(basal_plot, na.rm = T), mean_basal_area_SE = standard_error(basal_plot), .groups = "drop") %>% 
    mutate(across(.cols = everything(), \(x) replace_na(x, 0)))
    
  
  return(a)
  
}

# Calculate % species contribution to mean basal area by year and site
# Hard coded filters: red, white, black mangroves & Alive or Dead & Years 1,3
spp_contibution <- function(tree_measurements) {
  
  x <- tree_measurements %>%
    filter(Species %in% c("RHMA", "AVGE", "LARA")) %>%
    select(SY, Site, Plot, Species, DBH_cm) %>%
    group_by(SY, Site, Plot, Species) %>%
    summarise(basal_area_sum = sum(basal_area(DBH_cm), na.rm = T), .groups = "drop") %>% 
    pivot_wider(names_from = Species, values_from = basal_area_sum, values_fill = 0) %>%
    pivot_longer(cols = LARA:AVGE, names_to = "Species", values_to = "basal_area_sum") %>% 
    group_by(SY, Site, Species) %>% 
    summarise(mean_basal_area = mean(basal_area_sum), .groups = "drop") %>%
    group_by(SY, Site) %>% 
    mutate(contibution = (mean_basal_area / sum(mean_basal_area))*100) %>% 
    ungroup() %>% 
    pivot_wider(, id_cols = !mean_basal_area, names_from = Species, values_from = contibution)
  
  return(x)
  
}

# Calculate mean tree height(m) by year and site
# Hard coded filters: red, white, black, NA, UNK mangroves & Alive or Dead & Years 1,3
mean_tree_height <- function(tree_heights) {
  
  a <- tree_heights %>% 
    filter(Species %in% c("RHMA", "AVGE", "LARA", NA, "UNK")) %>% 
    select(SY, Site, Plot, Species, Tree_Height_m) %>% 
    group_by(SY, Site, Plot) %>% 
    summarise(mean_height_plot = mean(Tree_Height_m, na.rm = T), .groups = "drop") %>% 
    group_by(SY, Site) %>% 
    summarise(mean_height = mean(mean_height_plot, na.rm = T), mean_height_SE = standard_error(mean_height_plot), .groups = "drop") %>% 
    mutate(across(.cols = everything(), \(x) replace_na(x, 0)))
  
  return(a)
  
}

mean_stem_density <- function(tree_measurements) {
  
  a <- tree_measurements %>% 
    filter(Species %in% c("RHMA", "AVGE", "LARA") & Mortality %in% c("Alive", "Dying")) %>% 
    select(SY, Site, Plot, Species, DBH_cm) %>% 
    group_by(SY, Site, Plot) %>% 
    summarise(stem_plot_density = (length(DBH_cm))/.01, .groups = "drop") %>% 
    group_by(SY, Site) %>% 
    summarise(mean_stem_density = mean(stem_plot_density, na.rm = T), mean_stem_density_SE = standard_error(stem_plot_density), .groups = "drop")
  
  return(a)
}

mean_water_quality <- function(YSI) {
  
  a <- YSI %>% 
    select(SY, Site, Plot, Water_depth_cm, Temp, DO_mg, SPC, TDS) %>% 
    group_by(SY, Site, Plot) %>% 
    summarise(plot_depth = mean(Water_depth_cm, na.rm = TRUE),
              plot_temp = mean(Temp, na.rm = TRUE),
              plot_DO = mean(DO_mg, na.rm = TRUE),
              plot_SPC = mean(SPC, na.rm = TRUE),
              plot_TDS = mean(TDS, na.rm = TRUE),
              .groups = "drop") %>% 
    group_by(SY, Site) %>% 
    summarise(mean_depth = mean(plot_depth, na.rm = TRUE), mean_depth_SE = standard_error(plot_depth),
              mean_temp = mean(plot_temp, na.rm = TRUE), mean_temp_SE = standard_error(plot_temp),
              mean_DO = mean(plot_DO, na.rm = TRUE), mean_DO_SE = standard_error(plot_DO),
              mean_SPC = mean(plot_SPC, na.rm = TRUE), mean_SPC_SE = standard_error(plot_SPC),
              mean_TDS = mean(plot_TDS, na.rm = TRUE), mean_TDS_SE = standard_error(plot_TDS),
              .groups = "drop")
  
  return(a)
  
}

create_tree_measurement_table <- function(tree_measurements, tree_heights) {
  
  a <- mean_basal_area(tree_measurements)
  b <- spp_contibution(tree_measurements)
  c <- mean_tree_height(tree_heights)
  d <- mean_stem_density(tree_measurements)
  
  x <- reduce(.x = list(a,b,c,d), .f = full_join)
  
  format_x <- x %>% 
    mutate(across(where(is.numeric), round, 1)) %>% 
    mutate(mean_stem_density    = round(mean_stem_density, digits = 0),
           mean_stem_density_SE = round(mean_stem_density_SE, digits = 0)) %>% 
    mutate(
      Height = paste(mean_height, "±", mean_height_SE),
      Basal_Area = paste(mean_basal_area, "±", mean_basal_area_SE),
      Density = paste(mean_stem_density, "±", mean_stem_density_SE)
    ) %>% 
    select(SY, Site, Height, Basal_Area, Density, RHMA, AVGE, LARA)
  
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
      RHMA = md("*R. mangle* <br/> (%)"),
      AVGE = md("*A. germinians* <br/> (%)"),
      LARA = md("*L. racemosa* <br/> (%)"),
      Height = md("Height <br/> (m)"),
      Basal_Area = md("Basal area <br/> (m²/ha)"),
      Density = md("Density <br/> (stems/ha)")
    ) %>% 
    tab_spanner(label = "Relative Distribution of Species", columns = RHMA:LARA) %>% 
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

percent_basal_area_change <- function(tree_measurements) {
  
  a <- mean_basal_area(tree_measurements) %>% 
    select(SY, Site, mean_basal_area) %>%
    pivot_wider(names_from = SY, values_from = mean_basal_area, values_fill = 0) %>%
    mutate(percentage_change = ((`2024` - `2022`) / `2022` * 100))
  
  return(a)
}

site_coordinates <- function(Coordinates) {
  
  a <- Coordinates %>% 
    group_by(Island, Site) %>% 
    summarise(Typology = first(`Forest type`),Latitude = mean(Latitude), Longitude = mean(Longitude))
  
  return(a)
  
}

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
    xlab(label = "DBH Size Class (cm)") + 
    ylab(label = "Relative Frequency") + 
    ggtitle(site) +
    theme_Publication() +
    facet_wrap(~SY)
  
}

# Work in progress
# Currently a diverging bar chart with open plotted as line over top
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
    labs(y = "Forest Cover (%)", x = "Year", fill = "Cover Type", color = element_blank()) +
    theme_Publication() +
    scale_fill_manual(values = c("Alive" = "darkolivegreen3", "Dead" = "burlywood4")) +
    scale_color_manual(values = c("Open Area" = "deepskyblue"))
}

# This is working but needs some refactoring
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
    mutate(Island = factor(Island, levels = c("St Thomas", "St John", "St Croix")))
  
  
  
  height_rescale <- function(y) {
    min_h <- min(c$height_site, na.rm = TRUE)
    max_h <- max(c$height_site, na.rm = TRUE)
    min_m <- min(c$mean_site, na.rm = TRUE)
    max_m <- max(c$mean_site, na.rm = TRUE) + max(c$count_SE, na.rm = T)
    
    # Scale height_site to match the primary y-axis range
    (y - min_h) / (max_h - min_h) * (max_m - min_m) + min_m
  }
  
  ggplot(c, aes(x = Site, y = mean_site, fill = factor(SY))) + 
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
      if (!is.null(breaks)) scale_y_break(c(breaks, breaks), space = .025, scales = "free")
    } +
    labs(x = element_blank(), fill = "SY", color = "SY") +
    theme_Publication() +
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

sapling_density <- function(df) {
  
  a <- regen_data %>% 
    select(SY, Site, Plot) %>% 
    group_by(SY, Site) %>% 
    summarise(cnt = n_distinct(Plot))
  
  b <- sapling_data %>% 
    select(!c(Date, Notes)) %>% 
    group_by(SY, Site, Plot, Species) %>% 
    summarise(count = length(Height_cm), .groups = "drop") %>% 
    pivot_wider(names_from = Species, values_from = count, values_fill = 0)
  
}

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

water_qual_table <- function(YSI) {
  
  a <- mean_water_quality(YSI)
  
  format_x <- a %>% 
    mutate(across(where(is.numeric), round, 1)) %>% 
    mutate(mean_SPC = round(mean_SPC, digits = 0),
           mean_SPC_SE = round(mean_SPC_SE, digits = 0),
           mean_TDS = round(mean_TDS, digits = 0),
           mean_TDS_SE = round(mean_TDS_SE, digits = 0)) %>%
    mutate(
      Depth = paste(mean_depth, "±", mean_depth_SE),
      Temp = paste(mean_temp, "±", mean_temp_SE),
      DO = paste(mean_DO, "±", mean_DO_SE),
      SPC = paste(mean_SPC, "±", mean_SPC_SE),
      TDS = paste(mean_TDS, "±", mean_TDS_SE)
    ) %>% 
    select(SY, Site, Depth, Temp, DO, SPC, TDS)
  
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
      SPC = md("Specific Conductance <br/> (\u00b5S/cm)"),
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
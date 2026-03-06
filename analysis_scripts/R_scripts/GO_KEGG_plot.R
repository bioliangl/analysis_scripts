save_GO_plot <- function(go_data, save_name){
  library(tidyverse)
  library(scales)
  ## load GO data
  go <- read_csv(go_data)
  ## arrange raw data
  go_top <- go %>%
    group_by(ONTOLOGY) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    mutate(
      ONTOLOGY = recode(ONTOLOGY,
                        BP = "Biological Process",
                        CC = "Cellular Component",
                        MF = "Molecular Function"),
      Description_wrapped = str_wrap(Description, width = 50),
      Description_wrapped = factor(
        Description_wrapped,
        levels = rev(Description_wrapped[order(p.adjust)])
      ),
      logP = -log10(p.adjust)
    ) %>%
    ungroup()
  ## set figure size
  row_height <- 0.25
  extra_height <- 1.8
  plot_height <- nrow(go_top) * row_height + extra_height
  ## get figure
  p <- ggplot(go_top,
         aes(x = Count,
             y = Description_wrapped,
             fill = logP)) +
    geom_bar(stat = "identity", width = 0.8,) +
    facet_wrap(~ ONTOLOGY,
               ncol = 1,
               scales = "free_y",
               space = "free_y") + 
    scale_x_continuous(
      breaks = pretty_breaks(n = 5),
      labels = label_number(accuracy = 1),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_fill_gradient(
      low = "#4DBBD5",
      high = "#E64B35",
      name = "-log10(adj.P)"
    ) +
    guides(fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
    )) +
    labs(x = "Gene Count",
         y = NULL) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(
        face = "bold",
        size = 12
      ),
      axis.text.y = element_text(size = 10, colour = "black", face = "bold"),
      axis.text.x = element_text(size = 10, colour = "black", face = "bold"),
      axis.title.x = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 8, colour = "black"),
      legend.title = element_text(size = 10, colour = "black"),
      plot.margin = margin(5, 40, 5, 5)
    )
  ## save figure
  ggsave(
    filename = save_name,
    plot = p,
    width = 8,
    height = plot_height,
    dpi = 600,
    limitsize = FALSE
  )
}

save_KEGG_plot <- function(kegg_data, save_name){
  library(tidyverse)
  library(scales)
  ## load KEGG data
  kegg <- read_csv(kegg_data)
  ## arrange raw data
  kegg_top <- kegg %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    mutate(
      Description_wrapped = str_wrap(Description, width = 50),
      Description_wrapped = factor(
        Description_wrapped,
        levels = rev(Description_wrapped[order(p.adjust)])
      ),
      logP = -log10(p.adjust),
    ) %>%
    ungroup()
  ## set figure size
  row_height <- 0.25
  extra_height <- 1.8
  plot_height <- nrow(kegg_top) * row_height + extra_height
  ## get figure
  p <- ggplot(kegg_top,
              aes(x = Count,
                  y = Description_wrapped,
                  fill = logP)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_x_continuous(
      breaks = pretty_breaks(n = 5),
      labels = label_number(accuracy = 1),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_fill_gradient(
      low = "#4DBBD5",
      high = "#E64B35",
      name = "-log10(adj.P)"
    ) +
    guides(fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5
    )) +
    labs(x = "Gene Count",
         y = NULL) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(
        face = "bold",
        size = 12
      ),
      axis.text.y = element_text(size = 10, colour = "black", face = "bold"),
      axis.text.x = element_text(size = 10, colour = "black", face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 8, colour = "black"),
      legend.title = element_text(size = 10, colour = "black"),
      plot.margin = margin(5, 40, 5, 5)
    )
  ## save figure
  ggsave(
    filename = save_name,
    plot = p,
    width = 8,
    height = plot_height,
    dpi = 600,
    limitsize = FALSE
  )
}
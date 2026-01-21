#Model summary figures
#March 3, 2025
#Ana Miller-ter Kuile

#exploreing results of the Ebird model without lags

# Load packages -----------------------------------------------------------

package.list <- c('tidyverse',
                  'here', 'patchwork', 'coda', 
                  "FNN", 'ggtext',
                  'grid') 

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())
# Load model summary ------------------------------------------------------
# 
betas <- readRDS(here('monsoon',
                      'outputs',
                      'ebird_abund_model2_summary.RDS'))

# Weights -----------------------------------------------------------------

weights_df <- as.data.frame(betas$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "w")) %>%
  mutate(covariate = case_when(str_detect(parm, "wA") ~ "Cones",
                               str_detect(parm, "wB") ~ "Temperature",
                               str_detect(parm, 'wC') ~ "Precipitation")) %>%
  mutate(lag = str_sub(parm, 4, (nchar(parm)-1)),
         lag = as.numeric(lag)) 

# Figure 3: cone wts ------------------------------------------------------

start <- c(0.5,2.55, 3.55)
end <- c(2.45, 3.45, 5.5)
times <- c("Anticipatory", 
           "Immediate", 
           "Delayed")

box_df <- as.data.frame(start) %>%
  bind_cols(end = end, times = times) %>%
  mutate(times = factor(times, levels = c("Anticipatory", 
                                          "Immediate", 
                                          "Delayed")))

#get the figure icons to all be square
draw_square <- function(data, params, size) {
  if (is.null(data$size)) {
    data$size <- 0.5
  }
  lwd <- min(data$size, min(size) /4)
  grid::rectGrob(
    width  = unit(1, "snpc") - unit(lwd, "mm"),
    height = unit(1, "snpc") - unit(lwd, "mm"),
    gp = gpar(
      col = data$colour %||% NA,
      fill = alpha(data$fill %||% "grey20", data$alpha),
      lty = data$linetype %||% 1,
      lwd = lwd * .pt,
      linejoin = params$linejoin %||% "mitre",
      lineend = if (identical(params$linejoin, "round")) "round" else "square"
    )
  )
}
 
(cone_weight_plot <- weights_df %>%
  filter(covariate == "Cones") %>%
ggplot(aes(x = lag, y = `50%`)) +
  geom_rect(data = box_df, aes(xmin = start,
                               xmax = end,
                               group = times,
                               ymin = 0,
                               ymax = 1,
                               fill = times),
            inherit.aes = F, 
            alpha = 0.6,
            key_glyph = draw_square) +
  scale_fill_manual(values = c('#02818A',
                               '#67A9CF',
                               "#A6BDDB"),
                    labels = c("H1: **Anticipatory response** <br> (jays eat acorns and juniper <br> berries or see pinyon <br> cones develop)",
                               "H2: **Immediate response** <br> (jays eat seeds cached <br> the prior fall)",
                               "H3: **Delayed response** <br> (population response to <br> prior seed availability<br> and maturation of fledglings)")) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`,
                    ymax = `97.5%`),
                width = 0.2) +
    scale_x_reverse(breaks = c(1:5),
                    labels = c( "+2", "+1", "-1",
                                "-2", "-3")) +
  # scale_x_continuous(breaks = c(1:5),
  #                    labels = c( "+2", "+1", "-1", 
  #                               "-2", "-3")) +
  #annotate(geom = "text", x = 1.5, y = 1.1, label = "Cone years before \n birds surveyed") +
  #annotate(geom = "text", x = 4.5, y = 1.1, label = "Cone years after \n birds surveyed") +
  geom_vline(xintercept = 2.5, linetype = 2) +
    #geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 1)) +
  labs(x = "Cone abundance year (fall) \n (relative to spring of bird survey)",
       y = "Importance weight \n (median and 95% CI)",
       subtitle = c("Cone years before \nbirds surveyed", "Cone years after\n birds surveyed")) +
    theme(plot.margin = unit(c(0.5,0.25,0.25,0.25), "cm"), #starts at top and clockwise
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          legend.text = element_markdown(size = 12, 
                                         hjust = 0.5,
                                         margin = margin(t = 10)),
          legend.title = element_blank(),
          legend.key.width = unit(0.75, 'cm'), #change legend key width
          legend.key.spacing.y = unit(0.5, 'cm'),
          plot.subtitle = element_text(size = 12, 
                                       hjust = c(0, 1))) +
  guides(fill = guide_legend(title = "",
                             label.position = "right")))


ggsave(here('pictures',
            'R',
            'cone_weights.pdf'),
       width = 7,
       height = 4,
       units = 'in',
       dpi = 300)

ggsave(here('pictures',
            'final',
            'cone_weights.jpg'),
       width = 7,
       height = 4,
       units = 'in',
       dpi = 300)
#green - lag
#ccebc5
#blue - predictive
#b3cde3
#yellow - immediate
#ffffcc

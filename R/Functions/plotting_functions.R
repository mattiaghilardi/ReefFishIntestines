# Functions for plotting each panel of figure 3

# Standard length
#plot_fig3_sl <- function(model, data, ndraw, xlab, ylab, lim, col1, col2, trait) {
#  model %>% 
#    tidybayes::spread_draws(b_Intercept, b_scalespec_mean_sl_log, n = ndraw, seed = 29) %>% # extract ndraw random draws
#    mutate(spec_mean_sl_log = list(seq(min(data$spec_mean_sl_log), max(data$spec_mean_sl_log), 0.01))) %>% # the observed value range of body size
#    unnest(spec_mean_sl_log) %>%
#    mutate(pred = b_Intercept + (b_scalespec_mean_sl_log/sd(data$spec_mean_sl_log))*(spec_mean_sl_log - mean(data$spec_mean_sl_log)), # predict using unscaled slope and centred body size
#           isometry = b_Intercept + isometry[isometry$trait == trait, 2, drop = TRUE]*(spec_mean_sl_log - mean(data$spec_mean_sl_log))) %>% # isometric relationship with same intercept
#    group_by(spec_mean_sl_log) %>%
#    mutate(pred_m = mean(pred, na.rm = TRUE), # mean prediction
#           iso_m = mean(isometry, na.rm = TRUE)) %>% # mean isometry
#    ggplot(aes(x = exp(spec_mean_sl_log))) +
#    geom_line(aes(y = exp(pred), group = .draw), color = col2, alpha = 0.2) + # add ndraw predictions to show uncertainty around the mean
#    geom_line(aes(y = exp(iso_m)), color = "black", linetype = "longdash") + # add isometric relationship
#    geom_line(aes(y = exp(pred_m)), color = col1) + # add mean prediction
#    geom_rug(data = data, aes(x = exp(spec_mean_sl_log)), sides = "b", color = "black") + # add raw data as rug plot on the x axis
#    ylab(ylab) + xlab(xlab) +
#   scale_y_continuous(limits = lim)
#}

# Trophic level - also used for sensitivity plot in appendix S2
plot_fig3_troph <- function(model, data, ndraw, xlab, ylab, lim, col1, col2) {
  model %>% 
    tidybayes::spread_draws(b_Intercept, b_stomachPresent, b_scaletroph, n = ndraw, seed = 29) %>% # extract ndraw random draws
    mutate(troph = list(seq(min(data$troph), max(data$troph), 0.01))) %>% # the observed value range of trophic level
    unnest(troph) %>%
    mutate(pred = b_Intercept + b_stomachPresent + (b_scaletroph/sd(data$troph))*(troph - mean(data$troph))) %>% # predict using unscaled slope and centred trophic level
    group_by(troph) %>%
    mutate(pred_m = mean(pred, na.rm = TRUE)) %>% # mean prediction
    ggplot(aes(x = troph)) +
    geom_line(aes(y = exp(pred), group = .draw), color = col2, alpha = 0.2) + # add ndraw predictions to show uncertainty around the mean
    geom_line(aes(y = exp(pred_m)), color = col1) + # add mean prediction
    geom_rug(data = data, aes(x = troph), sides = "b", color = "black") + # add raw data as rug plot on the x axis
    ylab(ylab) + xlab(xlab) +
    scale_y_continuous(limits = lim)
}

# ELongation
plot_fig3_elon <- function(model, data, ndraw, xlab, ylab, lim, col1, col2) {
  model %>% 
    tidybayes::spread_draws(b_Intercept, b_stomachPresent, b_scaleelon_log, n = ndraw, seed = 29) %>% # extract ndraw random draws
    mutate(elon_log = list(seq(min(data$elon_log), max(data$elon_log), 0.01))) %>% # the observed value range of elongation
    unnest(elon_log) %>%
    mutate(pred = b_Intercept + b_stomachPresent + (b_scaleelon_log/sd(data$elon_log))*(elon_log - mean(data$elon_log))) %>% # predict using unscaled slope and centred elongation
    group_by(elon_log) %>%
    mutate(pred_m = mean(pred, na.rm = TRUE)) %>% # mean prediction
    ggplot(aes(x = exp(elon_log))) +
    geom_line(aes(y = exp(pred), group = .draw), color = col2, alpha = 0.2) + # add ndraw predictions to show uncertainty around the mean
    geom_line(aes(y = exp(pred_m)), color = col1) + # add mean prediction
    geom_rug(data = data, aes(x = exp(elon_log)), sides = "b", color = "black") + # add raw data as rug plot on the x axis
    ylab(ylab) + xlab(xlab) +
    scale_y_continuous(limits = lim)
}

#######################################################

# Function for plotting each panel of figure 4 and sensitivity plot in appendix S2
# xlab in html format
plot_fig4 <- function(data, model, xlab) {
  data %>%
    select(stomach, durophagy, troph, spec_mean_sl_log, elon_log, within_spec_sl_log) %>% 
    mutate(spec_mean_sl_log = mean(spec_mean_sl_log), 
           troph = mean(troph),
           elon_log = mean(elon_log),
           within_spec_sl_log = mean(within_spec_sl_log)) %>%
    unique() %>%
    tidybayes::add_fitted_draws(model, re_formula = NA) %>%
    mutate(stomach = factor(stomach, levels = c('Absent', 'Present'), labels = c('Stomach\nabsent', 'Stomach\npresent'))) %>%
    ggplot(aes(x = .value, y = stomach, color = durophagy, fill = durophagy)) +
    ggridges::geom_density_ridges(scale = 0.5, alpha = 0.4, rel_min_height = 0.01) + # add density plots filled by durophagy
    stat_pointinterval(aes(y = as.numeric(stomach) - 0.2), position = position_dodge(width = 0.3), 
                       .width = c(0.5, 0.95), scale = 0.5, interval_size_range = c(0.3, 0.8),
                       show.legend = FALSE) + # add intervals colored by durophagy
    xlab(xlab) + ylab("") +
    scale_color_viridis_d(option = "C", end = .8) +
    scale_fill_viridis_d(option = "C", end = .8) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm"),
          axis.title.x = element_markdown(size = 10))
}

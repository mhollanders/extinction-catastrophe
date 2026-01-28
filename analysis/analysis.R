# load packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(here, readxl, tidyverse, janitor, cmdstanr, tidybayes, posterior)
options(mc.cores = 8)

# read and organise data
files <- list.files(here("data"), ".csv")
dat_raw <- map(files,
               ~read_csv(here(str_c("data/", .x))) |> clean_names())
names(dat_raw) <- files

# process and filter data
dat <- dat_raw$`Table S1 - timeline.csv` |> 
  filter(!(record_id %in% str_c("pop_", c(2, 20, 31, 32, 36, 47, 55, 57))),
         record_certainty == "Confirmed") |> 
  mutate(across(c(scientific_name, predator), factor)) |> 
  select(prey = scientific_name, predator, contains("years_since"))

# prepare Stan data
scale_years <- dat |> 
  select(contains("years")) |> 
  as.matrix() |>
  c() |> 
  sd(na.rm = T)
stan_data <- dat |> 
  drop_na() |> 
  mutate(min_years_since = (min_years_since - 0.5) / scale_years,
         max_years_since = (max_years_since + 0.5) / scale_years) |> 
  compose_data() |> 
  glimpse()

# fit
mod <- cmdstan_model(here("analysis/mod.stan"))
fit <- mod$sample(stan_data, chains = 8, iter_warmup = 200, iter_sampling = 500)

# summarise
fit |> 
  spread_rvars(alpha[predator], epsilon[prey, predator]) |> 
  summarise(years = rvar_mean(scale_years * (alpha + epsilon)), .by = predator) |> 
  median_hdci(years)

# ggplot theme
my_theme <- function(base_size = 10,
                     base_family = "", 
                     base_line_size = base_size / 20 / 2, 
                     base_rect_size = base_size / 20) {
  my_black <- "#333333"
  half_line <- base_size / 2
  theme_grey(base_size = base_size,
             base_family = base_family, 
             base_line_size = base_line_size,
             base_rect_size = base_rect_size) %+replace%
    theme(axis.line = element_blank(), 
          axis.text = element_text(colour = my_black, size = rel(0.9)),
          axis.title = element_text(colour = my_black, size = rel(1)),
          axis.ticks =  element_line(colour = my_black),
          legend.key = element_rect(fill = "white", colour = NA),
          legend.text = element_text(size = rel(0.9)),
          panel.background = element_rect(fill = NA, colour = NA),
          panel.border = element_rect(fill = NA, colour = my_black),
          panel.grid = element_blank(),
          plot.margin = margin(10, 10, 10, 10),
          plot.title = element_text(size = rel(1.1), hjust = 0, vjust = 1, margin = margin(b = half_line)),
          strip.background = element_rect(fill = my_black, colour = my_black, linewidth = base_size / 2), 
          strip.text = element_text(colour = "white", size = rel(0.9), margin = margin(rep(0.8 * half_line, 4))),
          complete = TRUE)
}
theme_set(my_theme())

# plot
fit |> 
  spread_rvars(alpha[predator], epsilon[prey, predator]) |> 
  mutate(predator = factor(predator, labels = levels(dat$predator)),
         prey = factor(prey, labels = levels(dat$prey))) |> 
  left_join(dat |> 
              mutate(years = 0.5 * (min_years_since + max_years_since))) |> 
  ggplot(aes(x = years, 
             xdist = scale_years * (alpha + epsilon), 
             y = fct_rev(prey))) +
  facet_wrap(~ predator) + 
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red4") +
  geom_point(size = 1, alpha = 0.6, colour = "red4", shape = 16, 
             position = position_nudge(y = -0.3),
             show.legend = F) +
  stat_pointinterval(point_interval = median_hdi,
                     .width = 0.95,
                     linewidth = 0.5, size = 0.5) + 
  theme(axis.text.y = element_text(face = "italic")) + 
  labs(x = "Years of last sighting since predator arrival",
       y = NULL)
ggsave(here("figs/fig-years.png"), width = 8, height = 7, dpi = 600)

## ---- cluster-selection-figure ----
crit_mat <- tryCatch(
  as(kml_obj, "ListPartition")["criterionValuesAsMatrix"],
  error = function(e) kml_obj["criterionValuesAsMatrix"]
) %>%
  as.matrix()

# Ensure rows correspond to k = 2, 3, 4, and 5
if (nrow(crit_mat) == 20 && ncol(crit_mat) == 4) {
  crit_mat <- t(crit_mat)
}

# kml may return a list-matrix; convert every cell to numeric
crit_dim <- dim(crit_mat)
crit_dimnames <- dimnames(crit_mat)

crit_mat <- matrix(
  as.numeric(unlist(crit_mat, use.names = FALSE)),
  nrow = crit_dim[1],
  ncol = crit_dim[2],
  dimnames = crit_dimnames
)

# Extract k from row names such as c2, c3, c4, c5
k_values <- suppressWarnings(
  as.integer(gsub("[^0-9]", "", rownames(crit_mat)))
)

if (is.null(rownames(crit_mat)) ||
    length(k_values) != nrow(crit_mat) ||
    anyNA(k_values)) {
  k_values <- 2:(nrow(crit_mat) + 1)
}

crit_long <- as.data.frame(crit_mat, check.names = FALSE) %>%
  mutate(k = k_values) %>%
  filter(k %in% 2:5) %>%
  pivot_longer(
    cols = -k,
    names_to = "redrawing",
    values_to = "CH"
  ) %>%
  filter(is.finite(CH)) %>%
  group_by(k) %>%
  arrange(desc(CH), .by_group = TRUE) %>%
  mutate(
    redrawing_rank = row_number(),
    solution = factor(
      paste0(k, "-group solution"),
      levels = paste0(2:5, "-group solution")
    ),
    selected = factor(k == 2, levels = c(FALSE, TRUE))
  ) %>%
  ungroup()

if (nrow(crit_long) == 0) {
  stop("No finite Calinski–Harabasz values were extracted.")
}

ch_cols <- c(
  "2-group solution" = "#0072B2",
  "3-group solution" = "#D55E00",
  "4-group solution" = "#009E73",
  "5-group solution" = "grey55"
)

p_ch <- ggplot(
  crit_long,
  aes(
    x = redrawing_rank,
    y = CH,
    colour = solution,
    group = solution
  )
) +
  geom_line(
    aes(linewidth = selected, alpha = selected),
    lineend = "round"
  ) +
  geom_point(
    aes(size = selected, alpha = selected)
  ) +
  scale_colour_manual(values = ch_cols, name = NULL) +
  scale_linewidth_manual(
    values = c("FALSE" = 0.7, "TRUE" = 1.35),
    guide = "none"
  ) +
  scale_size_manual(
    values = c("FALSE" = 1.8, "TRUE" = 2.5),
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c("FALSE" = 0.65, "TRUE" = 1),
    guide = "none"
  ) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  labs(
    x = "Rank across repeated random initializations",
    y = "Calinski–Harabasz criterion"
  ) +
  theme_manuscript()

print(p_ch)

save_plot(
  p_ch,
  "Figure_4_cluster_selection_CH-2.png",
  width = 7.4,
  height = 4.6
)
save_plot(
  p_ch,
  "Figure_4_cluster_selection_CH-2.pdf",
  width = 7.4,
  height = 4.6
)

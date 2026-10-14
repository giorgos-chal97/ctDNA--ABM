# ============================================================
# REQUIRED PACKAGES
# ============================================================

# Run this ONCE if the packages are not already installed:

#install.packages(c("dplyr", "ggplot2", "tidyr", "scales"))
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)


# ============================================================
# SIMPLE INDIVIDUAL-CELL AGENT-BASED MODEL OF ctDNA
# ============================================================

set.seed(123)


# ============================================================
# 1. GENERAL SETTINGS
# ============================================================

simulation_days <- 120

initial_cells <- 1000

initial_resistant_fraction <- 0.05

treatment_days <- c(30, 60, 90)

treatment_duration <- 3


# ============================================================
# 2. BIOLOGICAL PARAMETERS
# ============================================================

# Sensitive cells
division_prob_S <- 0.025
apoptosis_prob_S <- 0.012
necrosis_prob_S <- 0.003

# Resistant cells
division_prob_R <- 0.022
apoptosis_prob_R <- 0.009
necrosis_prob_R <- 0.002


# ============================================================
# 3. TREATMENT EFFECT
# ============================================================

treatment_apoptosis_multiplier_S <- 6
treatment_necrosis_multiplier_S <- 3

treatment_apoptosis_multiplier_R <- 1.5
treatment_necrosis_multiplier_R <- 1.2


# ============================================================
# 4. ctDNA RELEASE
# ============================================================

active_release_S <- 0.002
active_release_R <- 0.003

apoptosis_release_S <- 1.0
apoptosis_release_R <- 1.2

necrosis_release_S <- 2.0
necrosis_release_R <- 2.4

ctDNA_clearance_fraction <- 0.80


# ============================================================
# 5. CREATE INITIAL AGENTS
# ============================================================

n_resistant <- round(
  initial_cells * initial_resistant_fraction
)

n_sensitive <- initial_cells - n_resistant

agents <- data.frame(
  id = 1:initial_cells,
  type = c(
    rep("Sensitive", n_sensitive),
    rep("Resistant", n_resistant)
  ),
  alive = TRUE,
  stringsAsFactors = FALSE
)

next_id <- initial_cells + 1


# ============================================================
# 6. STORAGE
# ============================================================

results <- data.frame()

ctDNA_pool <- 0


# ============================================================
# 7. SIMULATION LOOP
# ============================================================

for (day in 0:simulation_days) {
  
  # ----------------------------------------------------------
  # Determine whether treatment is active
  # ----------------------------------------------------------
  
  treatment_active <- any(
    day >= treatment_days &
      day < treatment_days + treatment_duration
  )
  
  
  # ----------------------------------------------------------
  # Identify living agents
  # ----------------------------------------------------------
  
  living_ids <- agents$id[
    agents$alive
  ]
  
  
  # ----------------------------------------------------------
  # Daily counters
  # ----------------------------------------------------------
  
  apoptosis_today <- 0
  necrosis_today <- 0
  
  new_cells <- list()
  
  
  # ==========================================================
  # EACH AGENT MAKES ITS OWN STOCHASTIC DECISION
  # ==========================================================
  
  for (cell_id in living_ids) {
    
    idx <- which(
      agents$id == cell_id
    )
    
    cell_type <- agents$type[idx]
    
    
    # --------------------------------------------------------
    # Assign biological probabilities
    # --------------------------------------------------------
    
    if (cell_type == "Sensitive") {
      
      division_prob <- division_prob_S
      apoptosis_prob <- apoptosis_prob_S
      necrosis_prob <- necrosis_prob_S
      
      if (treatment_active) {
        
        apoptosis_prob <-
          apoptosis_prob *
          treatment_apoptosis_multiplier_S
        
        necrosis_prob <-
          necrosis_prob *
          treatment_necrosis_multiplier_S
      }
      
    } else {
      
      division_prob <- division_prob_R
      apoptosis_prob <- apoptosis_prob_R
      necrosis_prob <- necrosis_prob_R
      
      if (treatment_active) {
        
        apoptosis_prob <-
          apoptosis_prob *
          treatment_apoptosis_multiplier_R
        
        necrosis_prob <-
          necrosis_prob *
          treatment_necrosis_multiplier_R
      }
    }
    
    
    # --------------------------------------------------------
    # Ensure probabilities do not exceed 1
    # --------------------------------------------------------
    
    apoptosis_prob <- min(
      apoptosis_prob,
      1
    )
    
    necrosis_prob <- min(
      necrosis_prob,
      1
    )
    
    
    # --------------------------------------------------------
    # APOPTOSIS
    # --------------------------------------------------------
    
    if (runif(1) < apoptosis_prob) {
      
      agents$alive[idx] <- FALSE
      apoptosis_today <- apoptosis_today + 1
      
      next
    }
    
    
    # --------------------------------------------------------
    # NECROSIS
    # --------------------------------------------------------
    
    if (runif(1) < necrosis_prob) {
      
      agents$alive[idx] <- FALSE
      necrosis_today <- necrosis_today + 1
      
      next
    }
    
    
    # --------------------------------------------------------
    # CELL DIVISION
    # --------------------------------------------------------
    
    if (runif(1) < division_prob) {
      
      new_cells[[length(new_cells) + 1]] <-
        data.frame(
          id = next_id,
          type = cell_type,
          alive = TRUE,
          stringsAsFactors = FALSE
        )
      
      next_id <- next_id + 1
    }
  }
  
  
  # ==========================================================
  # ADD DAUGHTER CELLS
  # ==========================================================
  
  if (length(new_cells) > 0) {
    
    agents <- bind_rows(
      agents,
      bind_rows(new_cells)
    )
  }
  
  
  # ==========================================================
  # CURRENT POPULATION
  # ==========================================================
  
  sensitive_alive <- sum(
    agents$alive &
      agents$type == "Sensitive"
  )
  
  resistant_alive <- sum(
    agents$alive &
      agents$type == "Resistant"
  )
  
  total_alive <-
    sensitive_alive +
    resistant_alive
  
  
  # ==========================================================
  # ctDNA RELEASE FROM LIVING CELLS
  # ==========================================================
  
  active_ctDNA <-
    sensitive_alive * active_release_S +
    resistant_alive * active_release_R
  
  
  # ==========================================================
  # ctDNA FROM CELL DEATH
  # ==========================================================
  
  death_ctDNA <-
    apoptosis_today *
    mean(
      c(
        apoptosis_release_S,
        apoptosis_release_R
      )
    ) +
    necrosis_today *
    mean(
      c(
        necrosis_release_S,
        necrosis_release_R
      )
    )
  
  
  # ==========================================================
  # ctDNA POOL
  # ==========================================================
  
  ctDNA_pool <-
    ctDNA_pool *
    (
      1 -
        ctDNA_clearance_fraction
    ) +
    active_ctDNA +
    death_ctDNA
  
  
  # ==========================================================
  # SAVE DAILY RESULTS
  # ==========================================================
  
  results <- bind_rows(
    results,
    data.frame(
      day = day,
      sensitive = sensitive_alive,
      resistant = resistant_alive,
      total_cells = total_alive,
      apoptosis = apoptosis_today,
      necrosis = necrosis_today,
      ctDNA = ctDNA_pool,
      treatment = treatment_active
    )
  )
}


# ============================================================
# 8. PRINT FINAL RESULTS
# ============================================================

cat(
  "\n====================================\n"
)

cat(
  "SIMULATION COMPLETE\n"
)

cat(
  "====================================\n"
)

cat(
  "Final sensitive cells:",
  tail(results$sensitive, 1),
  "\n"
)

cat(
  "Final resistant cells:",
  tail(results$resistant, 1),
  "\n"
)

cat(
  "Final total cells:",
  tail(results$total_cells, 1),
  "\n"
)

cat(
  "Final ctDNA:",
  round(
    tail(results$ctDNA, 1),
    2
  ),
  "\n"
)


# ============================================================
# 9. PLOT 1 - TOTAL TUMOR CELLS
# ============================================================

p1 <- ggplot(
  results,
  aes(
    x = day,
    y = total_cells
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_vline(
    xintercept = treatment_days,
    linetype = "dashed"
  ) +
  labs(
    title = "Tumor population over time",
    subtitle = "Dashed lines indicate treatment periods",
    x = "Day",
    y = "Living tumor cells"
  ) +
  theme_minimal(
    base_size = 13
  )

print(p1)


# ============================================================
# 10. PLOT 2 - SENSITIVE VS RESISTANT CELLS
# ============================================================

population_long <- results %>%
  select(
    day,
    sensitive,
    resistant
  ) %>%
  pivot_longer(
    cols = c(
      sensitive,
      resistant
    ),
    names_to = "cell_type",
    values_to = "cells"
  )

p2 <- ggplot(
  population_long,
  aes(
    x = day,
    y = cells,
    linetype = cell_type
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_vline(
    xintercept = treatment_days,
    linetype = "dotted"
  ) +
  labs(
    title = "Sensitive and resistant tumor-cell agents",
    x = "Day",
    y = "Number of living agents",
    linetype = "Cell type"
  ) +
  theme_minimal(
    base_size = 13
  )

print(p2)


# ============================================================
# 11. PLOT 3 - ctDNA
# ============================================================

p3 <- ggplot(
  results,
  aes(
    x = day,
    y = ctDNA
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_vline(
    xintercept = treatment_days,
    linetype = "dashed"
  ) +
  labs(
    title = "Simulated ctDNA dynamics",
    subtitle = "ctDNA released by viable and dying tumor-cell agents",
    x = "Day",
    y = "Simulated ctDNA units"
  ) +
  theme_minimal(
    base_size = 13
  )

print(p3)


# ============================================================
# 12. PLOT 4 - RESISTANT FRACTION
# ============================================================

results_fraction <- results %>%
  mutate(
    resistant_fraction =
      resistant / total_cells
  )

p4 <- ggplot(
  results_fraction,
  aes(
    x = day,
    y = resistant_fraction
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_vline(
    xintercept = treatment_days,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    labels = percent_format(
      accuracy = 1
    )
  ) +
  labs(
    title = "Resistant tumor-cell fraction over time",
    subtitle = "Dashed lines indicate treatment periods",
    x = "Day",
    y = "Resistant fraction"
  ) +
  theme_minimal(
    base_size = 13
  )

print(p4)


# ============================================================
# 13. SAVE OUTPUTS
# ============================================================

dir.create(
  "D:/ctdna/model2/plots",
  showWarnings = FALSE,
  recursive = TRUE
)

if (!dir.exists("D:/ctdna/model2/output_data")) {
  dir.create(
    "D:/ctdna/model2/output_data",
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# Save simulation results
# ------------------------------------------------------------

write.csv(
  results,
  "D:/ctdna/model2/output_data/Simple_ABM_results.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# Save plots
# ------------------------------------------------------------

ggsave(
  filename = "D:/ctdna/model2/plots/Simple_ABM_total_cells.png",
  plot = p1,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "D:/ctdna/model2/plots/Simple_ABM_sensitive_resistant.png",
  plot = p2,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "D:/ctdna/model2/plots/Simple_ABM_ctDNA.png",
  plot = p3,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "D:/ctdna/model2/plots/Simple_ABM_resistant_fraction.png",
  plot = p4,
  width = 8,
  height = 5,
  dpi = 300
)


# ============================================================
# END
# ============================================================

cat(
  "\nFiles saved successfully.\n"
)
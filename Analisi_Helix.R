library(PReMiuM)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, ggplot2)

load("exposome.RData")

air_pollution <- c("h_no2_ratio_preg_Log", "h_pm25_ratio_preg_None", 
                   "h_pm10_ratio_preg_None", "h_abs_ratio_preg_Log")
metals <- c("hs_pb_m_Log2", "hs_hg_m_Log2", "hs_cd_m_Log2", "hs_as_m_Log2")
pops <- c("hs_dde_madj_Log2", "hs_hcb_madj_Log2", "hs_pcb153_madj_Log2", "hs_pfos_m_Log2")
lifestyle <- c("e3_alcpreg_yn_None", "h_fish_preg_Ter", "h_fruit_preg_Ter", "h_pamod_t3_None")
environment <- c("h_popdens_preg_Sqrt", "h_ndvi100_preg_None", "h_greenyn300_preg_None", "h_lden_cat_preg_None")

selected_vars <- c(air_pollution, metals, pops, lifestyle, environment)

vars_to_tertiles <- c("h_no2_ratio_preg_Log", "h_pm25_ratio_preg_None", "h_pm10_ratio_preg_None", 
                      "h_abs_ratio_preg_Log", "hs_pb_m_Log2", "hs_hg_m_Log2", "hs_cd_m_Log2", 
                      "hs_as_m_Log2", "hs_dde_madj_Log2", "hs_hcb_madj_Log2", "hs_pcb153_madj_Log2", 
                      "hs_pfos_m_Log2", "h_popdens_preg_Sqrt", "h_ndvi100_preg_None", "h_lden_cat_preg_None")
vars_tertiles_recode <- c("h_fish_preg_Ter", "h_fruit_preg_Ter")
vars_binary <- c("e3_alcpreg_yn_None", "h_greenyn300_preg_None")
vars_categorical_ok <- c("h_pamod_t3_None")

categorize_tertiles <- function(x, var_name = "var") {
  q <- quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  if(length(unique(q)) < 4) {
    x_cat <- cut(x, breaks = 3, labels = c("Low", "Medium", "High"), include.lowest = TRUE)
  } else {
    x_cat <- cut(x, breaks = q, labels = c("Low", "Medium", "High"), include.lowest = TRUE)
  }
  return(x_cat)
}

expo_cat <- data.frame(row.names = rownames(exposome))

for(var in selected_vars) {
  if(var %in% vars_to_tertiles) {
    expo_cat[[var]] <- categorize_tertiles(exposome[[var]], var)
  } else if(var %in% vars_tertiles_recode) {
    old_levels <- levels(exposome[[var]])
    expo_cat[[var]] <- factor(exposome[[var]], levels = old_levels, labels = c("Low", "Medium", "High"))
  } else if(var %in% vars_binary) {
    expo_cat[[var]] <- factor(exposome[[var]], levels = c("0", "1"), labels = c("No", "Yes"))
  } else if(var %in% vars_categorical_ok) {
    expo_cat[[var]] <- exposome[[var]]
  } else {
    expo_cat[[var]] <- as.factor(exposome[[var]])
  }
}

outcome <- ifelse(as.numeric(phenotype$hs_bmi_c_cat) >= 3, 1, 0)

short_names <- c("NO2", "PM25", "PM10", "PMabs", "Pb", "Hg", "Cd", "As",
                 "DDE", "HCB", "PCB153", "PFOS", "Alcol", "Pesce", "Frutta", "Attivita",
                 "Densita", "NDVI", "Verde300", "Rumore")

data_premium <- expo_cat
names(data_premium) <- short_names
data_premium$outcome <- as.factor(outcome)

save(data_premium, file = "data_premium_ready.RData")

load("data_premium_ready.RData")

levels(data_premium$Frutta) <- list("Low-Med" = c("Low", "Medium"), "High" = "High")

lev_att <- levels(data_premium$Attivita)
levels(data_premium$Attivita) <- list(
  "Sedentary/Low" = c(lev_att[1], lev_att[2]),
  "Medium" = lev_att[3],
  "High" = lev_att[4]
)



save(data_premium, file = "data_premium_ready.RData")

out_file <- "exposome_premium"


levels(phenotype$hs_bmi_c_cat) 
round(prop.table(table(phenotype$hs_bmi_c_cat)) * 100, 2)

MCMC_CONFIG <- list(nSweeps = 20000, nBurn = 30000, nClusInit = 20, seed = 12345)

run_mcmc <- function() {
  load("data_premium_ready.RData")
  covariates_names <- names(data_premium)[1:20]
  data_input <- data.frame(outcome = as.numeric(as.character(data_premium$outcome)),
                           data_premium[, covariates_names])
  save(data_input, covariates_names, file = "data_input_premium.RData")
  
  old_files <- list.files(pattern = paste0("^", out_file, "_"))
  if(length(old_files) > 0) file.remove(old_files)
  
  set.seed(MCMC_CONFIG$seed)
  runInfoObj <- profRegr(yModel = "Bernoulli", xModel = "Discrete",
                         nSweeps = MCMC_CONFIG$nSweeps, nBurn = MCMC_CONFIG$nBurn, nFilter = 2,
                         data = data_input, outcome = "outcome", covNames = covariates_names,
                         nClusInit = MCMC_CONFIG$nClusInit, seed = MCMC_CONFIG$seed,
                         output = out_file, reportBurnIn = FALSE)
  save(runInfoObj, file = "runInfoObj.RData")
  return(runInfoObj)
}

run_postprocessing <- function() {
  load("runInfoObj.RData")
  dissimObj <- calcDissimilarityMatrix(runInfoObj)
  clusObj <- calcOptimalClustering(dissimObj)
  riskProfileObj <- calcAvgRiskAndProfile(clusObj)
  save(runInfoObj, dissimObj, clusObj, riskProfileObj, file = "premium_postprocessing.RData")
  return(list(dissimObj = dissimObj, clusObj = clusObj, riskProfileObj = riskProfileObj))
}

runInfoObj <- run_mcmc()
postproc <- run_postprocessing()

load("premium_postprocessing.RData")
load("data_premium_ready.RData")

dirichlet_phi <- function(var_data, all_levels, alpha_prior = 1) {
  counts <- table(factor(var_data, levels = all_levels))
  alpha_post <- alpha_prior + as.numeric(counts)
  phi_hat <- alpha_post / sum(alpha_post)
  names(phi_hat) <- all_levels
  return(phi_hat)
}

risk_matrix <- riskProfileObj$risk[, , 1]
n_clusters <- ncol(risk_matrix)

risk_mean <- colMeans(risk_matrix)
risk_sd <- apply(risk_matrix, 2, sd)
risk_q025 <- apply(risk_matrix, 2, quantile, probs = 0.025)
risk_q500 <- apply(risk_matrix, 2, quantile, probs = 0.50)
risk_q975 <- apply(risk_matrix, 2, quantile, probs = 0.975)

cluster_assignment <- clusObj$clustering
cluster_sizes <- as.numeric(table(cluster_assignment))
overall_prevalence <- mean(as.numeric(as.character(data_premium$outcome)))

risk_summary <- data.frame(
  Cluster = 1:n_clusters, N = cluster_sizes,
  Percent = round(cluster_sizes / sum(cluster_sizes) * 100, 1),
  Risk_Mean = round(risk_mean, 3), Risk_Median = round(risk_q500, 3),
  Risk_SD = round(risk_sd, 3), CI_2.5 = round(risk_q025, 3), CI_97.5 = round(risk_q975, 3)
)

write.csv(risk_summary, "Tabella_Rischio_Cluster.csv", row.names = FALSE)
save(risk_summary, risk_matrix, risk_mean, risk_sd, risk_q025, risk_q975,
     overall_prevalence, cluster_sizes, n_clusters, file = "risk_analysis_results.RData")

pdf("Rischio_Cluster.pdf", width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))

colors <- ifelse(risk_mean > overall_prevalence + 0.05, "firebrick3",
                 ifelse(risk_mean < overall_prevalence - 0.05, "forestgreen", "gray50"))
ord <- order(risk_mean)

bp <- barplot(risk_mean[ord] * 100, names.arg = paste("C", (1:n_clusters)[ord], sep = ""),
              col = colors[ord], ylim = c(0, max(risk_q975) * 100 + 10),
              ylab = "Probabilità Sovrappeso/Obesità (%)", xlab = "Cluster",
              main = "A. Rischio per Cluster", cex.names = 0.8)
arrows(bp, risk_q025[ord] * 100, bp, risk_q975[ord] * 100, angle = 90, code = 3, length = 0.05)
abline(h = overall_prevalence * 100, lty = 2, col = "blue", lwd = 2)
text(bp, 2, paste("n=", cluster_sizes[ord], sep = ""), cex = 0.6)

ord_desc <- order(risk_mean, decreasing = TRUE)
plot(risk_mean[ord_desc] * 100, 1:n_clusters, xlim = c(0, max(risk_q975) * 100 + 5),
     ylim = c(0.5, n_clusters + 0.5), pch = 19, cex = 1.5, col = colors[ord_desc],
     xlab = "Probabilità Sovrappeso/Obesità (%)", ylab = "", yaxt = "n", 
     main = "B. Forest Plot (IC 95% Credibilità)")
axis(2, at = 1:n_clusters, labels = paste("Cluster ", (1:n_clusters)[ord_desc], " (n=", cluster_sizes[ord_desc], ")", sep = ""), las = 1, cex.axis = 0.7)
segments(risk_q025[ord_desc] * 100, 1:n_clusters, risk_q975[ord_desc] * 100, 1:n_clusters, col = colors[ord_desc], lwd = 2)
abline(v = overall_prevalence * 100, lty = 2, col = "blue", lwd = 2)

barplot(cluster_sizes[ord], names.arg = paste("C", (1:n_clusters)[ord], sep = ""),
        col = colors[ord], ylab = "Numero Soggetti", xlab = "Cluster", 
        main = "C. Dimensione Cluster", cex.names = 0.8)

boxplot(risk_matrix * 100, names = paste("C", 1:n_clusters, sep = ""), col = colors,
        ylab = "Probabilità (%)", xlab = "Cluster", 
        main = "D. Distribuzione a Posteriori del Rischio", cex.axis = 0.8, outline = FALSE)
abline(h = overall_prevalence * 100, lty = 2, col = "blue", lwd = 2)
dev.off()

cov_names_original <- names(data_premium)[1:20]
cov_names_short <- c("NO2", "PM2.5", "PM10", "PM_abs", "Pb", "Hg", "Cd", "As",
                     "DDE", "HCB", "PCB153", "PFOS", "Alcol", "Pesce", "Frutta", "Attività",
                     "Densità", "NDVI", "Verde300", "Rumore")
families <- c(rep("Air Pollution", 4), rep("Metals", 4), rep("POPs", 4),
              rep("Lifestyle", 4), rep("Environment", 4))

high_risk_clusters <- c(3, 9, 10)
low_risk_clusters <- c(4, 7, 8)

calc_dirichlet_distribution <- function(data, cluster_ids, cov_index, cluster_assign) {
  subjects <- which(cluster_assign %in% cluster_ids)
  var_data <- data[subjects, cov_index]
  all_levels <- levels(data[, cov_index])
  phi_hat <- dirichlet_phi(var_data, all_levels)
  return(phi_hat)
}

comparison_results <- data.frame(Covariata = cov_names_short, Famiglia = families, stringsAsFactors = FALSE)

for(j in 1:20) {
  phi_high <- calc_dirichlet_distribution(data_premium, high_risk_clusters, j, cluster_assignment)
  phi_low <- calc_dirichlet_distribution(data_premium, low_risk_clusters, j, cluster_assignment)
  n_levels <- length(phi_high)
  
  if(n_levels == 3) {
    first_level <- names(phi_high)[1]
    comparison_results$High_pct_HighRisk[j] <- round(phi_high["High"] * 100, 1)
    comparison_results$High_pct_LowRisk[j] <- round(phi_low["High"] * 100, 1)
    comparison_results$Low_pct_HighRisk[j] <- round(phi_high[first_level] * 100, 1)
    comparison_results$Low_pct_LowRisk[j] <- round(phi_low[first_level] * 100, 1)
  } else if(n_levels == 2) {
    comparison_results$High_pct_HighRisk[j] <- round(phi_high[2] * 100, 1)
    comparison_results$High_pct_LowRisk[j] <- round(phi_low[2] * 100, 1)
    comparison_results$Low_pct_HighRisk[j] <- NA
    comparison_results$Low_pct_LowRisk[j] <- NA
  } else {
    comparison_results$High_pct_HighRisk[j] <- round(phi_high[n_levels] * 100, 1)
    comparison_results$High_pct_LowRisk[j] <- round(phi_low[n_levels] * 100, 1)
    comparison_results$Low_pct_HighRisk[j] <- round(phi_high[1] * 100, 1)
    comparison_results$Low_pct_LowRisk[j] <- round(phi_low[1] * 100, 1)
  }
}

comparison_results$Diff_High <- comparison_results$High_pct_HighRisk - comparison_results$High_pct_LowRisk
comparison_results <- comparison_results[order(abs(comparison_results$Diff_High), decreasing = TRUE), ]
comparison_results$Interpretazione <- ifelse(comparison_results$Diff_High > 5, "↑ Più alto nei cluster a RISCHIO",
                                             ifelse(comparison_results$Diff_High < -5, "↓ Più basso nei cluster a RISCHIO", "≈ Simile"))

write.csv(comparison_results, "Confronto_Profili_Corretto.csv", row.names = FALSE)

top_cov <- comparison_results[1:10, ]
png("Covariate_Discriminanti_Top10.png", width = 1000, height = 600, res = 120)
par(mar = c(8, 4, 4, 2))
bar_data <- rbind(top_cov$High_pct_HighRisk, top_cov$High_pct_LowRisk)
colors_bar <- c("firebrick3", "forestgreen")
bp <- barplot(bar_data, beside = TRUE, col = colors_bar, names.arg = top_cov$Covariata, las = 2,
              ylab = "φ(High) - Media a Posteriori (%)", ylim = c(0, max(bar_data, na.rm = TRUE) + 10),
              main = "")
text(bp[1,], bar_data[1,] + 2, labels = round(bar_data[1,], 0), cex = 0.7, col = "firebrick4")
text(bp[2,], bar_data[2,] + 2, labels = round(bar_data[2,], 0), cex = 0.7, col = "darkgreen")
legend("topright", legend = c("Cluster Alto Rischio", "Cluster Basso Rischio"), fill = colors_bar, bty = "n")
dev.off()

png("Differenze_Covariate.png", width = 900, height = 700, res = 120)
par(mar = c(5, 12, 4, 2))
comp_ordered <- comparison_results[order(comparison_results$Diff_High), ]
cols <- ifelse(comp_ordered$Diff_High > 0, "firebrick3", "forestgreen")
dotchart(comp_ordered$Diff_High, labels = paste(comp_ordered$Covariata, " (", comp_ordered$Famiglia, ")", sep = ""),
         pch = 19, color = cols, cex = 0.9,
         xlim = c(min(comp_ordered$Diff_High, na.rm = TRUE) - 5, max(comp_ordered$Diff_High, na.rm = TRUE) + 5),
         main = "",
         xlab = "Differenza (punti percentuali)")
abline(v = 0, lty = 2, col = "gray50", lwd = 2)
dev.off()

png("Heatmap_Cluster_Chiave.png", width = 1100, height = 600, res = 120)
key_clusters <- c(3, 9, 10, 8, 7, 4)
heatmap_data <- matrix(NA, nrow = length(key_clusters), ncol = 20)
rownames(heatmap_data) <- paste("C", key_clusters, " (", round(risk_mean[key_clusters] * 100, 0), "%)", sep = "")
colnames(heatmap_data) <- cov_names_short
for(i in 1:length(key_clusters)) {
  c_id <- key_clusters[i]
  subjects <- which(cluster_assignment == c_id)
  for(j in 1:20) {
    var_data <- data_premium[subjects, j]
    all_levels <- levels(data_premium[, j])
    phi_hat <- dirichlet_phi(var_data, all_levels)
    heatmap_data[i, j] <- phi_hat[length(phi_hat)] * 100
  }
}

# Layout
layout(matrix(c(1, 2), nrow = 1), widths = c(9, 1))

# --- Heatmap ---
par(mar = c(8, 10, 4, 1))
col_palette <- colorRampPalette(c("lightyellow", "orange", "orangered"))(100)
zlim <- c(0, max(heatmap_data, na.rm = TRUE))

image(1:20, 1:length(key_clusters), t(heatmap_data), col = col_palette,
      zlim = zlim, axes = FALSE, xlab = "", ylab = "", main = "")
axis(1, at = 1:20, labels = cov_names_short, las = 2, cex.axis = 0.7)
axis(2, at = 1:length(key_clusters), labels = rownames(heatmap_data), las = 1, cex.axis = 0.8)
for(i in 1:length(key_clusters)) {
  for(j in 1:20) {
    if(!is.na(heatmap_data[i, j])) text(j, i, round(heatmap_data[i, j], 0), cex = 0.6)
  }
}
abline(v = c(4.5, 8.5, 12.5, 16.5), col = "black", lwd = 1.5)
n_high_risk <- sum(key_clusters %in% high_risk_clusters)
abline(h = n_high_risk + 0.5, col = "blue", lwd = 2, lty = 2)

# --- Color bar ---
par(mar = c(8, 0.5, 4, 3))
legend_vals <- seq(zlim[1], zlim[2], length.out = 100)
image(1, legend_vals, t(matrix(legend_vals, ncol = 1)), col = col_palette,
      zlim = zlim, axes = FALSE, xlab = "", ylab = "")
axis(4, at = pretty(zlim), las = 1, cex.axis = 0.75)
mtext(expression(phi(High) ~ "%"), side = 4, line = 2, cex = 0.8)

dev.off()

psi_hat <- cluster_sizes / sum(cluster_sizes)

png("Psi_Barplot.png", width = 900, height = 500, res = 120)
par(mar = c(5, 5, 4, 2))
ord <- order(psi_hat, decreasing = TRUE)
bp <- barplot(psi_hat[ord] * 100, names.arg = paste("C", (1:n_clusters)[ord], sep = ""),
              col = "darkblue", ylim = c(0, max(psi_hat) * 100 + 5),
              ylab = expression(hat(psi)[c] ~ "(%)"), xlab = "Cluster",
              main = "")
text(bp, psi_hat[ord] * 100 + 1.5, labels = paste0(round(psi_hat[ord] * 100, 1), "%"), cex = 0.8)
text(bp, 0.5, labels = paste0("n=", cluster_sizes[ord]), cex = 0.7, col = "white")
dev.off()

logit <- function(p) { p <- pmax(pmin(p, 0.999), 0.001); log(p / (1 - p)) }
theta_matrix <- logit(risk_matrix)
theta_mean <- colMeans(theta_matrix)
ref_cluster <- which.min(risk_mean)

or_matrix <- matrix(NA, nrow = nrow(risk_matrix), ncol = n_clusters)
for(c in 1:n_clusters) or_matrix[, c] <- exp(theta_matrix[, c] - theta_matrix[, ref_cluster])

or_mean <- colMeans(or_matrix)
or_median <- apply(or_matrix, 2, median)
or_q025 <- apply(or_matrix, 2, quantile, probs = 0.025)
or_q975 <- apply(or_matrix, 2, quantile, probs = 0.975)

significativo <- ifelse(or_q025 > 1, "↑ Rischio aumentato",
                        ifelse(or_q975 < 1, "↓ Rischio ridotto", "Non significativo"))

or_table <- data.frame(
  Cluster = 1:n_clusters, N = cluster_sizes,
  Pct_Camp = round(psi_hat * 100, 1),
  Prob_pct = round(risk_mean * 100, 1),
  Prob_CI95 = paste0("[", round(risk_q025 * 100, 1), "; ", round(risk_q975 * 100, 1), "]"),
  Theta = round(theta_mean, 2),
  OR = round(or_mean, 2),
  OR_Median = round(or_median, 2),
  OR_CI95 = paste0("[", round(or_q025, 2), "; ", round(or_q975, 2), "]"),
  Significativo = significativo
)
write.csv(or_table[order(or_table$OR, decreasing = TRUE), ], "Odds_Ratio_Cluster.csv", row.names = FALSE)

png("Forest_Plot_OddsRatio.png", width = 1000, height = 700, res = 120)
par(mar = c(5, 12, 4, 2))
plot_data <- or_table[or_table$Cluster != ref_cluster, ]
plot_data <- plot_data[order(plot_data$OR), ]
n_plot <- nrow(plot_data)
colors_or <- ifelse(grepl("aumentato", plot_data$Significativo), "firebrick3",
                    ifelse(grepl("ridotto", plot_data$Significativo), "forestgreen", "gray50"))
x_min <- min(0.3, min(or_q025)); x_max <- max(5, max(or_q975))

plot(plot_data$OR, 1:n_plot, xlim = c(x_min, x_max), ylim = c(0.5, n_plot + 0.5),
     pch = 19, cex = 1.5, col = colors_or, xlab = "Odds Ratio (IC 95% Credibilità)", ylab = "", yaxt = "n", log = "x",
     main = paste("Forest Plot - Odds Ratio dalla Distribuzione anPosteriori\n(Riferimento: Cluster", ref_cluster, ")"))
axis(2, at = 1:n_plot, labels = paste("Cluster ", plot_data$Cluster, " (n=", plot_data$N, ")", sep = ""), las = 1, cex.axis = 0.75)
segments(or_q025[plot_data$Cluster], 1:n_plot, or_q975[plot_data$Cluster], 1:n_plot, col = colors_or, lwd = 2)
abline(v = 1, lty = 2, col = "blue", lwd = 2)
text(plot_data$OR, 1:n_plot, labels = round(plot_data$OR, 2), pos = 4, cex = 0.7, offset = 0.3)
legend("bottomright", legend = c("Rischio aumentato (p<0.05)", "Rischio ridotto (p<0.05)", "Non significativo", "OR = 1"),
       col = c("firebrick3", "forestgreen", "gray50", "blue"), pch = c(19, 19, 19, NA), lty = c(NA, NA, NA, 2), bty = "n", cex = 0.8)
dev.off()

impatto <- psi_hat * risk_mean
impatto_norm <- impatto / sum(impatto) * 100
ord_imp <- order(impatto_norm, decreasing = TRUE)

png("Rilevanza_Salute_Pubblica.png", width = 1100, height = 700, res = 120)
par(mar = c(7, 5, 4, 2))
cols_imp <- ifelse(risk_mean[ord_imp] > 0.35, "firebrick3",
                   ifelse(risk_mean[ord_imp] > 0.25, "darkorange", "forestgreen"))
y_max <- max(impatto_norm) + 8
bp <- barplot(impatto_norm[ord_imp], names.arg = paste("Cluster", (1:n_clusters)[ord_imp]),
              col = cols_imp, ylim = c(0, y_max), ylab = "Contributo relativo al totale obesità (%)", las = 2, cex.names = 0.9,
              main = "")
text(bp, impatto_norm[ord_imp] + 1.5, labels = paste0(round(impatto_norm[ord_imp], 1), "%"), cex = 0.85, font = 2)
mtext(side = 1, at = bp, line = 4.5, text = paste0("ψ=", round(psi_hat[ord_imp] * 100, 1), "%"), cex = 0.7, col = "steelblue")
mtext(side = 1, at = bp, line = 5.5, text = paste0("p=", round(risk_mean[ord_imp] * 100, 1), "%"), cex = 0.7, col = "firebrick")
legend("topright", legend = c("Rischio Alto (p > 35%)", "Rischio Medio (25-35%)", "Rischio Basso (p < 25%)"),
       fill = c("firebrick3", "darkorange", "forestgreen"), bty = "n", cex = 0.9)
contributo_uniforme <- 100 / n_clusters
abline(h = contributo_uniforme, lty = 2, col = "gray50", lwd = 1.5)
dev.off()

p_ref <- risk_mean[ref_cluster]
casi_totali <- sum(as.numeric(as.character(data_premium$outcome)))

paf_table <- data.frame(
  Cluster = 1:n_clusters, N = cluster_sizes,
  Psi_pct = round(psi_hat * 100, 1),
  Rischio_pct = round(risk_mean * 100, 1),
  Excess_Risk_pct = round((risk_mean - p_ref) * 100, 1),
  Casi_Eccesso = round(cluster_sizes * (risk_mean - p_ref), 1),
  PAF_pct = round(cluster_sizes * (risk_mean - p_ref) / casi_totali * 100, 1)
)
paf_table <- paf_table[order(paf_table$PAF_pct, decreasing = TRUE), ]
paf_table$PAF_Cumulativo <- cumsum(paf_table$PAF_pct)
write.csv(paf_table, "PAF_Cluster.csv", row.names = FALSE)

png("PAF_Bubble_Chart.png", width = 1000, height = 700, res = 120)
par(mar = c(5, 5, 4, 8))
bubble_size <- sqrt(paf_table$N) * 1.5
paf_colors <- colorRampPalette(c("forestgreen", "gold", "firebrick"))(100)
paf_scaled <- pmax(1, pmin(100, round((paf_table$PAF_pct + 5) / 25 * 99) + 1))

plot(paf_table$Rischio_pct, paf_table$PAF_pct, cex = bubble_size / 5, pch = 21,
     bg = paf_colors[paf_scaled], col = "black", xlim = c(10, 50),
     ylim = c(min(paf_table$PAF_pct) - 3, max(paf_table$PAF_pct) + 5),
     xlab = "Rischio nel Cluster - Media a Posteriori (%)", ylab = "PAF (%)",
     main = "")
abline(h = 0, lty = 2, col = "gray50")
abline(v = p_ref * 100, lty = 2, col = "blue")
text(paf_table$Rischio_pct, paf_table$PAF_pct, labels = paste("C", paf_table$Cluster, sep = ""), pos = 3, cex = 0.7, offset = 0.8)
legend("bottomright", legend = c("N = 50", "N = 100", "N = 200"), pt.cex = sqrt(c(50, 100, 200)) * 1.5 / 5,
       pch = 21, pt.bg = "gray80", title = "Dimensione Cluster", bty = "n", cex = 0.8, inset = c(-0.15, 0), xpd = TRUE)
dev.off()

subj_c3 <- which(cluster_assignment == 3)
subj_c8 <- which(cluster_assignment == 8)

cov_labels <- c("NO2", "PM2.5", "PM10", "PM abs", "Piombo (Pb)", "Mercurio (Hg)", "Cadmio (Cd)", "Arsenico (As)",
                "DDE", "HCB", "PCB153", "PFOS", "Alcol in gravidanza", "Consumo pesce", "Consumo frutta", "Attività fisica",
                "Densità popolazione", "NDVI (verde)", "Verde entro 300m", "Rumore traffico")

summary_table <- data.frame(Covariata = cov_labels, Famiglia = families, stringsAsFactors = FALSE)

for(i in 1:20) {
  all_levels <- levels(data_premium[, i])
  phi_c3 <- dirichlet_phi(data_premium[subj_c3, i], all_levels)
  phi_c8 <- dirichlet_phi(data_premium[subj_c8, i], all_levels)
  n_lev <- length(phi_c3)
  
  summary_table$Pct_High_C3[i] <- round(phi_c3[n_lev] * 100, 1)
  summary_table$Pct_High_C8[i] <- round(phi_c8[n_lev] * 100, 1)
  summary_table$Differenza[i] <- round(summary_table$Pct_High_C3[i] - summary_table$Pct_High_C8[i], 1)
  
  combined_tab <- rbind(
    table(factor(data_premium[subj_c3, i], levels = all_levels)),
    table(factor(data_premium[subj_c8, i], levels = all_levels))
  )
  chi_test <- tryCatch(chisq.test(combined_tab), error = function(e) list(p.value = NA))
  summary_table$P_value[i] <- round(chi_test$p.value, 4)
}
summary_table <- summary_table[order(abs(summary_table$Differenza), decreasing = TRUE), ]
write.csv(summary_table, "Confronto_Cluster3_vs_Cluster8.csv", row.names = FALSE)

png("Confronto_C3_vs_C8_Lollipop.png", width = 900, height = 800, res = 120)
par(mar = c(5, 14, 4, 2))
cols <- ifelse(summary_table$Differenza > 0, "firebrick3", "forestgreen")
dotchart(summary_table$Differenza, labels = paste(summary_table$Covariata, " (", summary_table$Famiglia, ")", sep = ""),
         pch = 19, color = cols, cex = 0.9, xlim = c(min(summary_table$Differenza) - 10, max(summary_table$Differenza) + 10),
         main = "", 
         xlab = "Differenza (punti percentuali)")
abline(v = 0, lty = 2, col = "gray50", lwd = 2)
mtext("← Più alto in C8 (protettivo)", side = 1, line = 3, adj = 0, cex = 0.8, col = "forestgreen")
mtext("Più alto in C3 (rischio) →", side = 1, line = 3, adj = 1, cex = 0.8, col = "firebrick3")
dev.off()

nomi_display <- c("NO2", "PM2.5", "PM10", "PM abs", "Piombo", "Mercurio", "Cadmio", "Arsenico",
                  "DDE", "HCB", "PCB153", "PFOS", "Alcol", "Pesce", "Frutta", "Attività",
                  "Densità", "NDVI", "Verde 300m", "Rumore")

dannoso <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
             TRUE, FALSE, FALSE, FALSE, NA, FALSE, FALSE, TRUE)

pct_alto_c3 <- numeric(20)
pct_alto_c8 <- numeric(20)

for(j in 1:20) {
  all_levels <- levels(data_premium[, j])
  phi_c3 <- dirichlet_phi(data_premium[subj_c3, j], all_levels)
  phi_c8 <- dirichlet_phi(data_premium[subj_c8, j], all_levels)
  pct_alto_c3[j] <- phi_c3[length(phi_c3)] * 100
  pct_alto_c8[j] <- phi_c8[length(phi_c8)] * 100
}

diff <- pct_alto_c3 - pct_alto_c8

coerenza <- character(20)
for(j in 1:20) {
  if(is.na(dannoso[j])) {
    coerenza[j] <- "Neutro"
  } else if(abs(diff[j]) < 5) {
    coerenza[j] <- "Simile"
  } else if(dannoso[j] && diff[j] > 0) {
    coerenza[j] <- "Coerente"
  } else if(!dannoso[j] && diff[j] < 0) {
    coerenza[j] <- "Coerente"
  } else {
    coerenza[j] <- "Incongruente"
  }
}

png("Lollipop_C3_C8_Corretto.png", width = 1100, height = 900, res = 120)
par(mar = c(5, 11, 4, 3))
ord <- order(abs(diff), decreasing = FALSE)
cols <- ifelse(coerenza == "Coerente", "steelblue",
               ifelse(coerenza == "Incongruente", "firebrick3",
                      ifelse(coerenza == "Simile", "gray50", "gray70")))
y_pos <- 1:20

plot(diff[ord], y_pos, xlim = c(-70, 110), ylim = c(0.5, 20.5), yaxt = "n", pch = 19, cex = 1.5,
     col = cols[ord], xlab = "Differenza φ(High) in punti percentuali", ylab = "",
     main = "Cluster 3 vs Cluster 8")

axis(2, at = y_pos, labels = FALSE, tick = TRUE, col = "gray30")
for(i in 1:length(y_pos)) {
  axis(2, at = y_pos[i], labels = nomi_display[ord][i], col.axis = cols[ord][i], las = 1, cex.axis = 0.85, tick = FALSE, font = 2)
}

abline(v = 0, lty = 2, col = "black", lwd = 2)
segments(0, y_pos, diff[ord], y_pos, col = cols[ord], lwd = 2.5)

for(i in 1:20) {
  d <- diff[ord[i]]
  label_text <- sprintf("%+.0f", d)
  if(d >= 0) text(d + 3, y_pos[i], labels = label_text, cex = 0.7, pos = 4, font = 2)
  else text(d - 3, y_pos[i], labels = label_text, cex = 0.7, pos = 2, font = 2)
}

text(55, 20.8, "Più alto in C3 →", col = "gray40", cex = 0.85, font = 3, xpd = TRUE)
text(-35, 20.8, "← Più alto in C8", col = "gray40", cex = 0.85, font = 3, xpd = TRUE)
abline(v = seq(-60, 100, 20), col = "gray90", lty = 1)
points(diff[ord], y_pos, pch = 19, cex = 1.5, col = cols[ord])
legend("bottomright", legend = c("Coerente con aspettative", "Incongruente", "Differenza minima (<5 pp)", "Neutro"),
       pch = 19, col = c("steelblue", "firebrick3", "gray50", "gray70"), pt.cex = 1.5, bty = "n", cex = 0.85)
dev.off()

save(theta_matrix, or_matrix, risk_matrix, or_table, paf_table, ref_cluster, file = "model_parameters.RData")

dirichlet_example <- function(cluster_id, var_name, data, cluster_assign) {
  subj <- which(cluster_assign == cluster_id)
  n_c <- length(subj)
  var_data <- data[[var_name]][subj]
  all_levels <- levels(data[[var_name]])
  counts <- table(factor(var_data, levels = all_levels))
  n_levels <- length(counts)
  alpha_prior <- rep(1, n_levels)
  alpha_post <- alpha_prior + as.numeric(counts)
  alpha_sum <- sum(alpha_post)
  phi_hat <- alpha_post / alpha_sum
  
  cat(sprintf("\n══════════════════════════════════════════════════════════════\n"))
  cat(sprintf("CLUSTER %d - Variabile: %s (n = %d)\n", cluster_id, var_name, n_c))
  cat(sprintf("══════════════════════════════════════════════════════════════\n\n"))
  
  cat("CONTEGGI OSSERVATI:\n")
  for(i in 1:n_levels) cat(sprintf("  n_{%d,%s,%s} = %d\n", cluster_id, var_name, all_levels[i], counts[i]))
  
  cat(sprintf("\nPRIOR: Dirichlet(%s)\n", paste(alpha_prior, collapse = ", ")))
  cat(sprintf("POSTERIOR: Dirichlet(%s)\n", paste(alpha_post, collapse = ", ")))
  
  cat("\nSTIME PUNTUALI:\n")
  for(i in 1:n_levels) cat(sprintf("  φ(%s) = %d/%d = %.1f%%\n", all_levels[i], alpha_post[i], alpha_sum, phi_hat[i] * 100))
  
  cat("\n--- LaTeX ---\n")
  cat("\\textit{Dati osservati:}\n\\begin{itemize}\n")
  for(i in 1:n_levels) cat(sprintf("    \\item $n_{%d,\\text{%s},\\text{%s}} = %d$\n", cluster_id, var_name, all_levels[i], counts[i]))
  cat("\\end{itemize}\n\n")
  cat(sprintf("\\textit{Posterior:} $\\boldsymbol{\\phi}_{%d,\\text{%s}} \\mid \\text{Dati} \\sim \\text{Dirichlet}(%s)$\n\n",
              cluster_id, var_name, paste(alpha_post, collapse = ", ")))
  cat("\\textit{Stime puntuali:}\n\\begin{align*}\n")
  for(i in 1:n_levels) cat(sprintf("\\hat{\\phi}_{%d,\\text{%s},\\text{%s}} &= \\frac{%d}{%d} \\approx %.1f\\%% \\\\\n",
                                   cluster_id, var_name, all_levels[i], alpha_post[i], alpha_sum, phi_hat[i] * 100))
  cat("\\end{align*}\n")
  
  return(invisible(list(counts = counts, alpha_post = alpha_post, phi_hat = phi_hat)))
}

cat("\n\n")
cat("████████████████████████████████████████████████████████████████████████████\n")
cat("██                    ESEMPI DIRICHLET PER TESI                           ██\n")
cat("████████████████████████████████████████████████████████████████████████████\n")

dirichlet_example(3, "PM10", data_premium, cluster_assignment)
dirichlet_example(3, "Frutta", data_premium, cluster_assignment)
dirichlet_example(3, "Attivita", data_premium, cluster_assignment)
dirichlet_example(8, "PM10", data_premium, cluster_assignment)

dirichlet_example(8, "Frutta", data_premium, cluster_assignment)











###############################################################################
##  DIAGNOSTIC MCMC – Profile Regression (PReMiuM)
##  Parameter: alpha, K, MMP (marginal model posterior)
###############################################################################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(coda)

# ── 0. Configuration ────────────────────────────────────────────────────────
out_file   <- "exposome_premium"          
n_sweeps   <- 10000                       

# ── 1. MMP ──────────────────────────────────

load("runInfoObj.RData")

mmp_file <- paste0(out_file, "_margModPost.txt")
if (!file.exists(mmp_file)) {
  cat("Calcolo Marginal Model Posterior (MMP)...\n")
  margModelPosterior(runInfoObj)
  cat("✓ MMP calcolata e salvata in:", mmp_file, "\n")
} else {
  cat("✓ File MMP già presente:", mmp_file, "\n")
}

# ── 2. ────────────────────────

alpha_file <- paste0(out_file, "_alpha.txt")
if (!file.exists(alpha_file)) stop("File non trovato: ", alpha_file)
alpha_chain <- scan(alpha_file, quiet = TRUE)

nclus_file <- paste0(out_file, "_nClusters.txt")
if (!file.exists(nclus_file)) stop("File non trovato: ", nclus_file)
K_chain <- scan(nclus_file, quiet = TRUE)

mmp_chain <- tryCatch({
  tmp <- read.table(mmp_file, header = FALSE)[, 1]
  tmp
}, error = function(e) {
  tryCatch({
    tmp <- read.table(mmp_file, header = TRUE)[, 1]
    tmp
  }, error = function(e2) {
    scan(mmp_file, quiet = TRUE)
  })
})

len <- min(length(alpha_chain), length(K_chain), length(mmp_chain))
alpha_chain <- alpha_chain[1:len]
K_chain     <- K_chain[1:len]
mmp_chain   <- mmp_chain[1:len]

cat(sprintf("\nLunghezza catene (allineate): %d iterazioni post burn-in\n", len))
cat(sprintf("  Alpha: media = %.3f, sd = %.3f\n", mean(alpha_chain), sd(alpha_chain)))
cat(sprintf("  K:     media = %.2f, sd = %.2f\n", mean(K_chain), sd(K_chain)))
cat(sprintf("  MMP:   media = %.2f, sd = %.2f\n", mean(mmp_chain), sd(mmp_chain)))

# ── 3. MCMC conversion ─────────────────────────────────────────────
alpha_mcmc <- mcmc(alpha_chain)
K_mcmc     <- mcmc(K_chain)
mmp_mcmc   <- mcmc(mmp_chain)

# ── 4. GEWEKE TEST ────────────────────────────────────────────────────────

geweke_alpha <- geweke.diag(alpha_mcmc, frac1 = 0.1, frac2 = 0.5)
geweke_K     <- geweke.diag(K_mcmc,     frac1 = 0.1, frac2 = 0.5)
geweke_mmp   <- geweke.diag(mmp_mcmc,   frac1 = 0.1, frac2 = 0.5)

geweke_z <- c(geweke_alpha$z, geweke_K$z, geweke_mmp$z)
geweke_p <- 2 * pnorm(-abs(geweke_z))

interpret_geweke <- function(z, alpha_level = 0.05) {
  p_val <- 2 * pnorm(-abs(z))
  if (abs(z) <= 1.96) {
    return(sprintf("z = %.3f, p = %.4f, |z| ≤ 1.96 → CONVERGENZA OK", z, p_val))
  } else {
    return(sprintf("z = %.3f, p = %.4f, |z| > 1.96 → ATTENZIONE: possibile non-convergenza", z, p_val))
  }
}

cat("\n")
cat("══════════════════════════════════════════════════════════════\n")
cat("              TEST DI GEWEKE (α = 0.05, soglia |z| = 1.96)  \n")
cat("══════════════════════════════════════════════════════════════\n\n")
cat("Alpha: ", interpret_geweke(geweke_alpha$z), "\n")
cat("K:     ", interpret_geweke(geweke_K$z), "\n")
cat("MMP:   ", interpret_geweke(geweke_mmp$z), "\n")

# ── 5. EFFECTIVE SAMPLE SIZE (ESS) ──────────────────────────────────────────
ess_alpha <- effectiveSize(alpha_mcmc)
ess_K     <- effectiveSize(K_mcmc)
ess_mmp   <- effectiveSize(mmp_mcmc)

cat("\n")
cat("══════════════════════════════════════════════════════════════\n")
cat("           EFFECTIVE SAMPLE SIZE (ESS)                       \n")
cat("══════════════════════════════════════════════════════════════\n\n")
cat(sprintf("Alpha:  ESS = %.1f  (su %d iterazioni → efficienza = %.1f%%)\n",
            ess_alpha, len, ess_alpha / len * 100))
cat(sprintf("K:      ESS = %.1f  (su %d iterazioni → efficienza = %.1f%%)\n",
            ess_K, len, ess_K / len * 100))
cat(sprintf("MMP:    ESS = %.1f  (su %d iterazioni → efficienza = %.1f%%)\n",
            ess_mmp, len, ess_mmp / len * 100))

for (param in c("Alpha", "K", "MMP")) {
  ess_val <- switch(param, Alpha = ess_alpha, K = ess_K, MMP = ess_mmp)
  if (ess_val >= 400) {
    cat(sprintf("  %s: ESS = %.0f ≥ 400 → ECCELLENTE\n", param, ess_val))
  } else if (ess_val >= 200) {
    cat(sprintf("  %s: ESS = %.0f ≥ 200 → ACCETTABILE\n", param, ess_val))
  } else if (ess_val >= 100) {
    cat(sprintf("  %s: ESS = %.0f ≥ 100 → MARGINALE (considerare più iterazioni)\n", param, ess_val))
  } else {
    cat(sprintf("  %s: ESS = %.0f < 100 → INSUFFICIENTE (aumentare nSweeps o thinning)\n", param, ess_val))
  }
}

# ── 6. TAB ─────────────────────────────────────────────────
ess_vals <- c(ess_alpha, ess_K, ess_mmp)

diag_table <- data.frame(
  Parametro   = c("Alpha (concentrazione)", "K (n. cluster)", "MMP (marginal model posterior)"),
  Media       = c(mean(alpha_chain), mean(K_chain), mean(mmp_chain)),
  SD          = c(sd(alpha_chain), sd(K_chain), sd(mmp_chain)),
  ESS         = round(ess_vals, 1),
  Efficienza  = paste0(round(ess_vals / len * 100, 1), "%"),
  Geweke_z    = round(geweke_z, 3),
  Geweke_p    = round(geweke_p, 4),
  Convergenza = ifelse(abs(geweke_z) <= 1.96, "OK", "Attenzione")
)

write.csv(diag_table, "Tabella_Diagnostica_MCMC.csv", row.names = FALSE)
cat("\n✓ Tabella diagnostica salvata: Tabella_Diagnostica_MCMC.csv\n")

# ── 7. TRACE PLOT + CORRELOGRAM ───────────────────────
pdf("Diagnostica_Panel_Completo.pdf", width = 14, height = 9)
layout(matrix(1:6, nrow = 2, byrow = TRUE))
par(mar = c(4, 5, 3, 1), cex.main = 1.1)

# --- Trace plots ---

# Alpha
plot(as.numeric(alpha_mcmc), type = "l", col = "purple", lwd = 0.3,
     main = expression("A. Trace plot – " * alpha * " (concentrazione)"),
     xlab = "Iterazione (post burn-in)", ylab = expression(alpha))
abline(h = mean(alpha_chain), col = "firebrick", lty = 2, lwd = 1.5)
legend("topright", legend = c(
  sprintf("Media = %.3f", mean(alpha_chain)),
  sprintf("ESS = %.0f", ess_alpha),
  sprintf("Geweke z = %.3f (p = %.3f)", geweke_alpha$z, geweke_p[1])
), bty = "n", cex = 0.8)

# K
plot(as.numeric(K_mcmc), type = "l", col = "forestgreen", lwd = 0.3,
     main = "B. Trace plot – K (n. cluster attivi)",
     xlab = "Iterazione (post burn-in)", ylab = "K")
abline(h = mean(K_chain), col = "firebrick", lty = 2, lwd = 1.5)
legend("topright", legend = c(
  sprintf("Media = %.2f", mean(K_chain)),
  sprintf("ESS = %.0f", ess_K),
  sprintf("Geweke z = %.3f (p = %.3f)", geweke_K$z, geweke_p[2])
), bty = "n", cex = 0.8)

# MMP
plot(as.numeric(mmp_mcmc), type = "l", col = "darkblue", lwd = 0.3,
     main = "C. Trace plot – MMP (Marginal Model Posterior)",
     xlab = "Iterazione (post burn-in)", ylab = "log MMP")
abline(h = mean(mmp_chain), col = "firebrick", lty = 2, lwd = 1.5)
legend("topright", legend = c(
  sprintf("Media = %.2f", mean(mmp_chain)),
  sprintf("ESS = %.0f", ess_mmp),
  sprintf("Geweke z = %.3f (p = %.3f)", geweke_mmp$z, geweke_p[3])
), bty = "n", cex = 0.8)

# --- Riga 2: Correlogram ---
max_lag <- min(200, len - 1)

acf(alpha_chain, lag.max = max_lag,
    main = expression("D. Autocorrelazione – " * alpha),
    col = "purple", lwd = 1.2)

acf(K_chain, lag.max = max_lag,
    main = "E. Autocorrelazione – K",
    col = "forestgreen", lwd = 1.2)

acf(mmp_chain, lag.max = max_lag,
    main = "F. Autocorrelazione – MMP",
    col = "darkblue", lwd = 1.2)

dev.off()
cat("✓ Panel diagnostico salvato: Diagnostica_Panel_Completo.pdf\n")

# ── 8. GEWEKE PLOT ─────────────────────────────────────
pdf("Diagnostica_Geweke.pdf", width = 12, height = 10)
par(mfrow = c(3, 1), mar = c(4, 5, 3, 2))

n_windows <- 20
frac2 <- 0.5
start_fracs <- seq(0.05, 0.45, length.out = n_windows)

geweke_windows <- function(chain_mcmc, start_fracs, frac2) {
  z_scores <- numeric(length(start_fracs))
  for (i in seq_along(start_fracs)) {
    gd <- geweke.diag(chain_mcmc, frac1 = start_fracs[i], frac2 = frac2)
    z_scores[i] <- gd$z
  }
  return(z_scores)
}

# Alpha
z_alpha <- geweke_windows(alpha_mcmc, start_fracs, frac2)
plot(start_fracs * 100, z_alpha, type = "b", pch = 19, col = "purple",
     xlab = "Inizio prima finestra (%)", ylab = "z-score di Geweke",
     main = expression("Geweke diagnostic – " * alpha),
     ylim = c(min(c(z_alpha, -3)), max(c(z_alpha, 3))))
abline(h = c(-1.96, 1.96), lty = 2, col = "firebrick", lwd = 1.5)
abline(h = 0, lty = 1, col = "gray50")
legend("topright", legend = "Soglia |z| = 1.96", col = "firebrick", lty = 2, bty = "n")

# K
z_K <- geweke_windows(K_mcmc, start_fracs, frac2)
plot(start_fracs * 100, z_K, type = "b", pch = 19, col = "forestgreen",
     xlab = "Inizio prima finestra (%)", ylab = "z-score di Geweke",
     main = "Geweke diagnostic – K",
     ylim = c(min(c(z_K, -3)), max(c(z_K, 3))))
abline(h = c(-1.96, 1.96), lty = 2, col = "firebrick", lwd = 1.5)
abline(h = 0, lty = 1, col = "gray50")

# MMP
z_mmp <- geweke_windows(mmp_mcmc, start_fracs, frac2)
plot(start_fracs * 100, z_mmp, type = "b", pch = 19, col = "darkblue",
     xlab = "Inizio prima finestra (%)", ylab = "z-score di Geweke",
     main = "Geweke diagnostic – MMP",
     ylim = c(min(c(z_mmp, -3)), max(c(z_mmp, 3))))
abline(h = c(-1.96, 1.96), lty = 2, col = "firebrick", lwd = 1.5)
abline(h = 0, lty = 1, col = "gray50")

dev.off()
cat("✓ Geweke plot salvati: Diagnostica_Geweke.pdf\n")

# ── 9. Save ─────────────────────────────────────────────
save(alpha_chain, K_chain, mmp_chain,
     alpha_mcmc, K_mcmc, mmp_mcmc,
     geweke_alpha, geweke_K, geweke_mmp,
     ess_alpha, ess_K, ess_mmp,
     diag_table,
     file = "diagnostica_mcmc_results.RData")

cat("\n")
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat("                     RIEPILOGO DIAGNOSTICA MCMC                             \n")
cat("══════════════════════════════════════════════════════════════════════════════\n")
print(diag_table, row.names = FALSE)
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat("  DIAGNOSTICA COMPLETATA\n")
cat("══════════════════════════════════════════════════════════════════════════════\n")


rm(list = ls())

library(nlme)
library(splines)
library(ggplot2)
library(dplyr)
library(readxl)
library(effsize)
library(car)
library(tidyr)

# =====================================================
# LABEL MAPS — đổi tên chỉ số cho hình ảnh
# =====================================================

label_en <- c(
  EW     = "Edge Width",
  GP     = "Gradient Peak",
  Blur   = "Edge Degradation",
  Weight = "Weight",
  WL     = "Weight Loss (%)"
)

label_vi <- c(
  EW     = "Độ rộng cạnh (Edge Width)",
  GP     = "Đỉnh gradient (Gradient Peak)",
  Blur   = "Suy giảm cạnh (Edge Degradation)",
  Weight = "Khối lượng (Weight)",
  WL     = "Hao hụt khối lượng (%)"
)

save_bilingual <- function(plot_en, plot_vi, base_name) {
  ggsave(paste0("Figures/", base_name, "_EN.png"),
         plot_en, width = 10, height = 8, dpi = 400)
  ggsave(paste0("Figures/", base_name, "_VI.png"),
         plot_vi, width = 10, height = 8, dpi = 400)
  cat("[Saved]", base_name, "_EN.png &", base_name, "_VI.png\n")
}

# =====================================================
# 1. LOAD DATA
# =====================================================

file_path <- file.choose()
data <- read_excel(file_path)

cat("=== Cấu trúc file gốc ===\n")
cat("Số cột:", ncol(data), "\n")
cat("Số dòng:", nrow(data), "\n")

colnames(data) <- c("STT_mau", "SampleID", "Mark", "Condition", "Day",
                    "EW_value", "GP_value", "Blur_value",
                    "EW", "GP", "Blur", "Weight", "WL")

data$SampleID  <- as.factor(data$SampleID)
data$Day       <- as.numeric(data$Day)
data$Condition <- factor(data$Condition, levels = c("Room", "Cold"))
data$EW        <- as.numeric(data$EW)
data$GP        <- as.numeric(data$GP)
data$Blur      <- as.numeric(data$Blur)
data$Weight    <- as.numeric(data$Weight)
data$WL        <- as.numeric(data$WL)

data$Sample_Group <- as.factor(paste0(data$Mark, data$STT_mau))
data$IndivID      <- as.character(paste0(data$Mark, data$STT_mau))

cat("\n=== Kiểm tra cấu trúc ===\n")
cat("Số cá thể unique:", length(unique(data$Sample_Group)), "(kỳ vọng: 60)\n")
data %>% select(Sample_Group, Condition) %>% distinct() %>% count(Condition) %>% print()

overlap <- data %>%
  select(Sample_Group, Condition) %>% distinct() %>%
  group_by(Sample_Group) %>% filter(n() > 1)
cat("SampleID ở 2 conditions:", nrow(overlap), "(kỳ vọng: 0)\n")

cat("\nWL range:", range(data$WL, na.rm = TRUE), "\n")
cat("WL mean:",  mean(data$WL,  na.rm = TRUE), "\n")

dir.create("Figures", showWarnings = FALSE)
dir.create("Tables",  showWarnings = FALSE)

# =====================================================
# 2. T-TEST NGÀY ĐẦU TIÊN
# =====================================================

day1_data <- data %>% filter(Day == 1)

perform_ttest <- function(data_day1, var_name, var_key) {

  lbl_en <- label_en[var_key]
  lbl_vi <- label_vi[var_key]

  cat("\n=== T-TEST:", lbl_en, "===\n")

  room_vals <- data_day1 %>% filter(Condition == "Room") %>% pull(!!sym(var_name))
  cold_vals <- data_day1 %>% filter(Condition == "Cold") %>% pull(!!sym(var_name))
  room_vals <- room_vals[!is.na(room_vals)]
  cold_vals <- cold_vals[!is.na(cold_vals)]

  if (length(room_vals) < 2 || length(cold_vals) < 2) {
    cat("Không đủ dữ liệu!\n"); return(NULL)
  }

  desc_stats <- data_day1 %>%
    group_by(Condition) %>%
    summarise(n    = sum(!is.na(!!sym(var_name))),
              Mean = round(mean(!!sym(var_name), na.rm = TRUE), 3),
              SD   = round(sd(!!sym(var_name),   na.rm = TRUE), 3),
              .groups = "drop")

  sw_room_p <- tryCatch(shapiro.test(room_vals)$p.value, error = function(e) NA_real_)
  sw_cold_p <- tryCatch(shapiro.test(cold_vals)$p.value, error = function(e) NA_real_)

  test_df <- data.frame(
    Value     = c(room_vals, cold_vals),
    Condition = factor(c(rep("Room", length(room_vals)), rep("Cold", length(cold_vals))),
                       levels = c("Room", "Cold"))
  )

  levene_result <- leveneTest(Value ~ Condition, data = test_df)
  levene_p      <- round(levene_result$`Pr(>F)`[1], 4)
  use_equal_var <- levene_p > 0.05
  ttest_result  <- t.test(room_vals, cold_vals, var.equal = use_equal_var)
  cohen_d       <- cohen.d(room_vals, cold_vals)

  result_table <- data.frame(
    Parameter = c("n (Room)", "n (Cold)", "Mean +/- SD (Room)", "Mean +/- SD (Cold)",
                  "Shapiro p (Room)", "Shapiro p (Cold)", "Levene p", "Test used",
                  "t-statistic", "df", "p-value", "CI95 lower", "CI95 upper",
                  "Cohen's d", "Effect size"),
    Value = c(
      desc_stats$n[1], desc_stats$n[2],
      paste0(desc_stats$Mean[1], " +/- ", desc_stats$SD[1]),
      paste0(desc_stats$Mean[2], " +/- ", desc_stats$SD[2]),
      ifelse(is.na(sw_room_p), "NA", round(sw_room_p, 4)),
      ifelse(is.na(sw_cold_p), "NA", round(sw_cold_p, 4)),
      levene_p,
      ifelse(use_equal_var, "Student's t", "Welch's t"),
      round(ttest_result$statistic, 4),
      round(ttest_result$parameter, 2),
      round(ttest_result$p.value, 4),
      round(ttest_result$conf.int[1], 4),
      round(ttest_result$conf.int[2], 4),
      round(cohen_d$estimate, 4),
      as.character(cohen_d$magnitude)
    )
  )

  write.csv(result_table, paste0("Tables/Ttest_", var_key, "_Day1.csv"), row.names = FALSE)

  subtitle_text <- paste0("p = ", round(ttest_result$p.value, 4),
                          "  |  Cohen's d = ", round(cohen_d$estimate, 3))

  base_plot <- ggplot(test_df, aes(x = Condition, y = Value, fill = Condition)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
    scale_fill_manual(values = c("Room" = "#d95f5f", "Cold" = "#2b6cb0")) +
    theme_classic(base_size = 18) +
    theme(plot.background  = element_rect(fill = "grey97"),
          panel.background = element_rect(fill = "grey97"),
          legend.position  = "none")

  # --- Ban tieng Anh ---
  p_en <- base_plot +
    labs(title    = paste0("T-test: ", lbl_en, " at Day 1"),
         subtitle = subtitle_text,
         x = "Condition", y = lbl_en)

  # --- Ban tieng Viet ---
  cond_vi <- c("Room" = "Nhiet do phong", "Cold" = "Lanh")
  test_df_vi <- test_df
  test_df_vi$Condition <- factor(cond_vi[as.character(test_df$Condition)],
                                  levels = c("Nhiet do phong", "Lanh"))
  p_vi <- ggplot(test_df_vi, aes(x = Condition, y = Value, fill = Condition)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
    scale_fill_manual(values = c("Nhiet do phong" = "#d95f5f", "Lanh" = "#2b6cb0")) +
    theme_classic(base_size = 18) +
    labs(title    = paste0(lbl_vi, " tai Ngay 1"),
         subtitle = subtitle_text,
         x = "Dieu kien bao quan", y = lbl_vi) +
    theme(plot.background  = element_rect(fill = "grey97"),
          panel.background = element_rect(fill = "grey97"),
          legend.position  = "none")

  save_bilingual(p_en, p_vi, paste0("Ttest_", var_key, "_Day1"))
  print(p_en)
}

perform_ttest(day1_data, "Weight", "Weight")
perform_ttest(day1_data, "EW",     "EW")
perform_ttest(day1_data, "GP",     "GP")
perform_ttest(day1_data, "Blur",   "Blur")

# =====================================================
# 3. SPLINE MIXED MODELS
# =====================================================

cat("\n=== FITTING SPLINE LMM ===\n")

data_complete <- data %>%
  filter(!is.na(Day), !is.na(Condition),
         !is.na(EW), !is.na(GP), !is.na(Blur), !is.na(Weight), !is.na(WL))

model_EW     <- lme(EW     ~ bs(Day, df = 4) * Condition, random = ~1 | Sample_Group,
                    data = data_complete, na.action = na.omit)
model_GP     <- lme(GP     ~ bs(Day, df = 4) * Condition, random = ~1 | Sample_Group,
                    data = data_complete, na.action = na.omit)
model_Blur   <- lme(Blur   ~ bs(Day, df = 4) * Condition, random = ~1 | Sample_Group,
                    data = data_complete, na.action = na.omit)
model_Weight <- lme(Weight ~ bs(Day, df = 4) * Condition, random = ~1 | Sample_Group,
                    data = data_complete, na.action = na.omit)
model_WL     <- lme(WL     ~ bs(Day, df = 4) * Condition, random = ~1 | Sample_Group,
                    data = data_complete, na.action = na.omit)

cat("All models fitted.\n")

# =====================================================
# R2 MARGINAL/CONDITIONAL -- Nakagawa & Schielzeth 2013
# =====================================================

cat("\n=== R2 MARGINAL/CONDITIONAL ===\n")

compute_r2_nlme <- function(model, var_name) {
  fitted_fixed <- predict(model, level = 0)
  var_fixed    <- var(fitted_fixed)
  rand_var     <- as.numeric(VarCorr(model)[1, "Variance"])
  resid_var    <- model$sigma^2
  total_var    <- var_fixed + rand_var + resid_var

  R2m <- round(var_fixed / total_var, 4)
  R2c <- round((var_fixed + rand_var) / total_var, 4)

  cat(sprintf("  %-8s R2m = %.4f | R2c = %.4f | Var_fixed=%.4f | Var_rand=%.4f | Var_resid=%.4f\n",
              var_name, R2m, R2c, var_fixed, rand_var, resid_var))

  data.frame(
    Variable       = var_name,
    R2_marginal    = R2m,
    R2_conditional = R2c,
    Var_fixed      = round(var_fixed,  4),
    Var_random     = round(rand_var,   4),
    Var_residual   = round(resid_var,  4),
    R2m_interpret  = ifelse(R2m >= 0.75, "Excellent",
                     ifelse(R2m >= 0.50, "Good",
                     ifelse(R2m >= 0.25, "Moderate", "Weak"))),
    R2c_interpret  = ifelse(R2c >= 0.75, "Excellent",
                     ifelse(R2c >= 0.50, "Good",
                     ifelse(R2c >= 0.25, "Moderate", "Weak")))
  )
}

r2_EW     <- compute_r2_nlme(model_EW,     "EW")
r2_GP     <- compute_r2_nlme(model_GP,     "GP")
r2_Blur   <- compute_r2_nlme(model_Blur,   "Blur")
r2_Weight <- compute_r2_nlme(model_Weight, "Weight")
r2_WL     <- compute_r2_nlme(model_WL,     "WL")

r2_table <- bind_rows(r2_EW, r2_GP, r2_Blur, r2_Weight, r2_WL)

write.csv(r2_table, "Tables/Spline_R2_Nakagawa.csv", row.names = FALSE)
cat("[Saved] Tables/Spline_R2_Nakagawa.csv\n")

cat("\n=== TONG HOP CHAT LUONG SPLINE MODEL ===\n")
print(r2_table %>% select(Variable, R2_marginal, R2m_interpret,
                           R2_conditional, R2c_interpret))

# =====================================================
# MODEL DIAGNOSTICS
#   - RMSE
#   - Residual SD (sigma)
#   - Residual normality: Shapiro-Wilk
#   - Homoscedasticity: Spearman(|residuals|, fitted)
# =====================================================

cat("\n=== MODEL DIAGNOSTICS ===\n")

compute_diagnostics <- function(model, var_name) {

  resid_vals  <- residuals(model, type = "response")
  fitted_vals <- fitted(model)

  # --- RMSE ---
  rmse_val <- round(sqrt(mean(resid_vals^2)), 4)

  # --- Residual SD (sigma tu model) ---
  resid_sd <- round(model$sigma, 4)

  # --- Shapiro-Wilk: kiem tra phan phoi chuan cua residuals ---
  # Lay mau toi da 5000 neu n lon (gioi han cua shapiro.test la 5000)
  set.seed(42)
  resid_sample <- if (length(resid_vals) > 5000)
                    sample(resid_vals, 5000) else resid_vals
  sw        <- shapiro.test(resid_sample)
  sw_stat   <- round(sw$statistic, 4)
  sw_p      <- round(sw$p.value,   4)
  sw_normal <- ifelse(sw$p.value > 0.05, "Yes (p>0.05)", "No (p<=0.05)")

  # --- Homoscedasticity: Spearman(|residuals|, fitted values) ---
  # p > 0.05 va |rho| nho: khong co pattern he thong => phuong sai dong deu
  cor_res <- cor.test(fitted_vals, abs(resid_vals),
                      method = "spearman", exact = FALSE)
  rf_rho  <- round(cor_res$estimate, 4)
  rf_p    <- round(cor_res$p.value,  4)
  rf_homo <- ifelse(cor_res$p.value > 0.05,
                    "Homoscedastic (p>0.05)",
                    "Heteroscedastic (p<=0.05)")

  cat(sprintf(
    "  %-8s RMSE=%.4f | sigma=%.4f | SW W=%.4f p=%.4f [%s] | Res~Fit rho=%.4f p=%.4f [%s]\n",
    var_name, rmse_val, resid_sd,
    sw_stat, sw_p, sw_normal,
    rf_rho, rf_p, rf_homo
  ))

  data.frame(
    Variable             = var_name,
    RMSE                 = rmse_val,
    Residual_SD          = resid_sd,
    Shapiro_W            = sw_stat,
    Shapiro_p            = sw_p,
    Normality            = sw_normal,
    ResFit_Spearman_rho  = rf_rho,
    ResFit_Spearman_p    = rf_p,
    Homoscedasticity     = rf_homo,
    N_obs                = length(resid_vals)
  )
}

diag_EW     <- compute_diagnostics(model_EW,     "EW")
diag_GP     <- compute_diagnostics(model_GP,     "GP")
diag_Blur   <- compute_diagnostics(model_Blur,   "Blur")
diag_Weight <- compute_diagnostics(model_Weight, "Weight")
diag_WL     <- compute_diagnostics(model_WL,     "WL")

diag_table <- bind_rows(diag_EW, diag_GP, diag_Blur, diag_Weight, diag_WL)

# Luu bang diagnostics rieng
write.csv(diag_table, "Tables/Spline_Model_Diagnostics.csv", row.names = FALSE)
cat("[Saved] Tables/Spline_Model_Diagnostics.csv\n")

# Luu bang tong hop diagnostics + R2 (de dung trong bai bao)
diag_full <- diag_table %>%
  left_join(r2_table %>% select(Variable, R2_marginal, R2_conditional,
                                 R2m_interpret, R2c_interpret),
            by = "Variable")

write.csv(diag_full, "Tables/Spline_Model_Full_Evaluation.csv", row.names = FALSE)
cat("[Saved] Tables/Spline_Model_Full_Evaluation.csv\n")

cat("\n=== TONG HOP DIAGNOSTICS ===\n")
print(diag_table %>% select(Variable, RMSE, Residual_SD,
                              Shapiro_p, Normality,
                              ResFit_Spearman_rho, ResFit_Spearman_p,
                              Homoscedasticity))

# =====================================================
# PREDICTION GRID & SPLINE PLOTS
# =====================================================

day_range <- range(data_complete$Day, na.rm = TRUE)
newdata <- expand.grid(
  Day       = seq(day_range[1], day_range[2], length.out = 200),
  Condition = levels(data$Condition)
)

newdata$EW_pred     <- predict(model_EW,     newdata, level = 0)
newdata$GP_pred     <- predict(model_GP,     newdata, level = 0)
newdata$Blur_pred   <- predict(model_Blur,   newdata, level = 0)
newdata$Weight_pred <- predict(model_Weight, newdata, level = 0)
newdata$WL_pred     <- predict(model_WL,     newdata, level = 0)

plot_dual_spline <- function(original_data, pred_data, yvar, ypred,
                              title_text, subtitle_text,
                              xlab, ylab,
                              cond_levels, cond_labels, cond_colors) {

  od <- original_data
  pd <- pred_data
  od$Condition_disp <- factor(cond_labels[as.character(od$Condition)],
                               levels = cond_labels)
  pd$Condition_disp <- factor(cond_labels[as.character(pd$Condition)],
                               levels = cond_labels)

  ggplot() +
    geom_point(data = od,
               aes(x = Day, y = .data[[yvar]], color = Condition_disp),
               size = 3, alpha = 0.8) +
    geom_line(data = pd,
              aes(x = Day, y = .data[[ypred]], color = Condition_disp),
              linewidth = 1.6) +
    scale_color_manual(values = cond_colors) +
    theme_classic(base_size = 18) +
    labs(title    = title_text,
         subtitle = subtitle_text,
         x = xlab, y = ylab) +
    theme(plot.background  = element_rect(fill = "grey97"),
          panel.background = element_rect(fill = "grey97"),
          legend.title     = element_blank())
}

make_spline_bilingual <- function(yvar, ypred, var_key, r2_obj) {

  sub_txt <- paste0("R2m = ", r2_obj$R2_marginal,
                    " | R2c = ", r2_obj$R2_conditional)

  cond_en <- c("Room" = "Room", "Cold" = "Cold")
  col_en  <- c("Room" = "#d95f5f", "Cold" = "#2b6cb0")

  p_en <- plot_dual_spline(
    data_complete, newdata, yvar, ypred,
    title_text    = paste0("Spline LMM - ", label_en[var_key]),
    subtitle_text = sub_txt,
    xlab = "Day", ylab = label_en[var_key],
    cond_levels = c("Room", "Cold"),
    cond_labels = cond_en,
    cond_colors = col_en
  )

  cond_vi <- c("Room" = "Nhiet do phong", "Cold" = "Lanh")
  col_vi  <- c("Nhiet do phong" = "#d95f5f", "Lanh" = "#2b6cb0")

  p_vi <- plot_dual_spline(
    data_complete, newdata, yvar, ypred,
    title_text    = paste0("Spline LMM - ", label_vi[var_key]),
    subtitle_text = sub_txt,
    xlab = "Ngay", ylab = label_vi[var_key],
    cond_levels = c("Room", "Cold"),
    cond_labels = cond_vi,
    cond_colors = col_vi
  )

  save_bilingual(p_en, p_vi, paste0("Spline_", var_key))
  print(p_en)
}

make_spline_bilingual("EW",     "EW_pred",     "EW",     r2_EW)
make_spline_bilingual("GP",     "GP_pred",     "GP",     r2_GP)
make_spline_bilingual("Blur",   "Blur_pred",   "Blur",   r2_Blur)
make_spline_bilingual("Weight", "Weight_pred", "Weight", r2_Weight)
make_spline_bilingual("WL",     "WL_pred",     "WL",     r2_WL)

# =====================================================
# 4. COUPLING ANALYSIS -- LME (VAR ~ WL)
# =====================================================

cat("\n=== COUPLING ANALYSIS (LME) ===\n")

run_lme_wl <- function(response_var, var_key, data_lme) {
  cat("\n---", label_en[var_key], "vs Weight Loss ---\n")
  results_list <- list()

  for (cond in c("Overall", "Room", "Cold")) {
    sub <- if (cond == "Overall") data_lme else
           data_lme %>% filter(Condition == cond)
    sub <- sub %>% filter(!is.na(WL), !is.na(.data[[response_var]]))
    if (nrow(sub) < 10) next

    fit <- tryCatch(
      lme(as.formula(paste(response_var, "~ WL")),
          random = ~1 | Sample_Group, data = sub, na.action = na.omit),
      error = function(e) NULL
    )
    if (is.null(fit)) { cat("LME khong hoi tu:", var_key, cond, "\n"); next }

    fit_sum   <- summary(fit)
    fix_tab   <- fit_sum$tTable
    r_var     <- as.numeric(VarCorr(fit)[1, "Variance"])
    res_var   <- fit$sigma^2
    icc_val   <- round(r_var / (r_var + res_var), 4)
    cor_res   <- cor.test(sub$WL, sub[[response_var]], method = "spearman", exact = FALSE)
    ci_slope  <- intervals(fit, which = "fixed")$fixed
    ci_lower  <- round(ci_slope["WL", "lower"], 4)
    ci_upper  <- round(ci_slope["WL", "upper"], 4)
    slope_val <- fix_tab["WL", "Value"]
    slope_p   <- fix_tab["WL", "p-value"]

    cat(sprintf("  [%s] slope=%.4f [CI: %.4f, %.4f] | p=%s | rho=%.4f | ICC=%.4f\n",
                cond, slope_val, ci_lower, ci_upper,
                formatC(slope_p, format = "e", digits = 2),
                cor_res$estimate, icc_val))

    results_list[[length(results_list) + 1]] <- data.frame(
      Condition     = cond,
      Response      = var_key,
      n_obs         = nrow(sub),
      n_individuals = length(unique(sub$Sample_Group)),
      Spearman_rho  = round(cor_res$estimate, 4),
      Spearman_p    = cor_res$p.value,
      Spearman_sig  = ifelse(cor_res$p.value < 0.001, "***",
                     ifelse(cor_res$p.value < 0.01,   "**",
                     ifelse(cor_res$p.value < 0.05,   "*", "ns"))),
      LME_intercept = round(fix_tab["(Intercept)", "Value"], 4),
      LME_slope_WL  = round(slope_val, 4),
      CI95_lower    = ci_lower,
      CI95_upper    = ci_upper,
      LME_SE_slope  = round(fix_tab["WL", "Std.Error"], 4),
      LME_t_slope   = round(fix_tab["WL", "t-value"], 4),
      LME_p_slope   = slope_p,
      LME_sig       = ifelse(slope_p < 0.001, "***",
                     ifelse(slope_p < 0.01,   "**",
                     ifelse(slope_p < 0.05,   "*", "ns"))),
      ICC           = icc_val,
      AIC           = round(AIC(fit), 2)
    )
  }
  bind_rows(results_list)
}

lme_EW   <- run_lme_wl("EW",   "EW",   data_complete)
lme_GP   <- run_lme_wl("GP",   "GP",   data_complete)
lme_Blur <- run_lme_wl("Blur", "Blur", data_complete)

all_lme_results <- bind_rows(lme_EW, lme_GP, lme_Blur)

# =====================================================
# COUPLING PLOTS (EN + VI)
# =====================================================

cat("\n=== Coupling Plots ===\n")

plot_coupling <- function(data_raw, response_var, var_key, fig_num, results,
                           lang = "EN") {

  lbl       <- if (lang == "EN") label_en[var_key] else label_vi[var_key]
  wl_lbl    <- if (lang == "EN") "Weight Loss (%)" else "Hao hut khoi luong (%)"
  title_txt <- if (lang == "EN") {
    paste0(lbl, " vs Weight Loss (%)")
  } else {
    paste0(lbl, " vs Hao hut khoi luong (%)")
  }

  cond_map <- if (lang == "EN") {
    c("Room" = "Room", "Cold" = "Cold")
  } else {
    c("Room" = "Nhiet do phong", "Cold" = "Lanh")
  }
  cond_colors <- if (lang == "EN") {
    c("Room" = "#d95f5f", "Cold" = "#2b6cb0")
  } else {
    c("Nhiet do phong" = "#d95f5f", "Lanh" = "#2b6cb0")
  }

  data_plot <- data_raw
  data_plot$Condition_disp <- factor(cond_map[as.character(data_raw$Condition)],
                                      levels = cond_map)

  stats_ov   <- results %>% filter(Condition == "Overall", Response == var_key)
  stats_room <- results %>% filter(Condition == "Room",    Response == var_key)
  stats_cold <- results %>% filter(Condition == "Cold",    Response == var_key)

  room_lbl <- if (lang == "EN") "Room" else "Nhiet do phong"
  cold_lbl <- if (lang == "EN") "Cold" else "Lanh"

  p <- ggplot(data_plot, aes(x = WL, y = .data[[response_var]], color = Condition_disp)) +
    geom_point(size = 2, alpha = 0.4) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.4, alpha = 0.2) +
    scale_color_manual(values = cond_colors) +
    theme_classic(base_size = 16) +
    labs(
      title    = title_txt,
      subtitle = paste0("Overall: rho=", stats_ov$Spearman_rho,
                        " | slope=", stats_ov$LME_slope_WL,
                        " [", stats_ov$CI95_lower, ", ", stats_ov$CI95_upper, "]",
                        " | p=", formatC(stats_ov$LME_p_slope, format = "e", digits = 2)),
      x = wl_lbl, y = lbl
    ) +
    theme(plot.background  = element_rect(fill = "grey97"),
          panel.background = element_rect(fill = "grey97"),
          legend.title     = element_blank(),
          legend.position  = c(0.15, 0.85))

  if (nrow(stats_room) > 0)
    p <- p + annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 2,
                      label = paste0(room_lbl, ": rho=", stats_room$Spearman_rho,
                                     " slope=", stats_room$LME_slope_WL,
                                     " [", stats_room$CI95_lower, ", ", stats_room$CI95_upper, "]"),
                      color = "#d95f5f", size = 3.5, fontface = "bold")

  if (nrow(stats_cold) > 0)
    p <- p + annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 4,
                      label = paste0(cold_lbl, ": rho=", stats_cold$Spearman_rho,
                                     " slope=", stats_cold$LME_slope_WL,
                                     " [", stats_cold$CI95_lower, ", ", stats_cold$CI95_upper, "]"),
                      color = "#2b6cb0", size = 3.5, fontface = "bold")

  return(p)
}

for (i in seq_along(c("EW", "GP", "Blur"))) {
  var_keys <- c("EW", "GP", "Blur")
  vk <- var_keys[i]
  p_en <- plot_coupling(data_complete, vk, vk, i, all_lme_results, lang = "EN")
  p_vi <- plot_coupling(data_complete, vk, vk, i, all_lme_results, lang = "VI")
  save_bilingual(p_en, p_vi, paste0("Fig", i, "_", vk, "_vs_WL"))
  print(p_en)
}

# =====================================================
# 5. XUAT KET QUA
# =====================================================

write.csv(all_lme_results, "Tables/Coupling_LME_Results.csv", row.names = FALSE)

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list(
    Coupling_Results  = all_lme_results,
    Model_R2          = r2_table,
    Model_Diagnostics = diag_table,
    Model_Full_Eval   = diag_full
  ), "Tables/All_Results.xlsx")
  cat("[Saved] Tables/All_Results.xlsx\n")
} else {
  cat("[Note] install.packages('writexl') de xuat Excel\n")
}

cat("[Saved] Tables/Coupling_LME_Results.csv\n")

# =====================================================
# IN KET QUA TONG HOP CUOI
# =====================================================

cat("\n=============================================\n")
cat("           KET QUA TONG HOP\n")
cat("=============================================\n")

cat("\n[1] CHAT LUONG SPLINE MODEL (Nakagawa & Schielzeth 2013):\n")
print(r2_table %>%
        select(Variable, R2_marginal, R2m_interpret,
               R2_conditional, R2c_interpret))

cat("\n[2] MODEL DIAGNOSTICS:\n")
print(diag_table %>%
        select(Variable, RMSE, Residual_SD,
               Shapiro_p, Normality,
               ResFit_Spearman_rho, ResFit_Spearman_p,
               Homoscedasticity))

cat("\n[3] COUPLING ANALYSIS - slope + 95% CI + Spearman:\n")
print(all_lme_results %>%
        select(Condition, Response,
               Spearman_rho, Spearman_sig,
               LME_slope_WL, CI95_lower, CI95_upper, LME_sig,
               ICC))

cat("\nDONE!\n")
cat("\n=== DANH SACH FILE BANG (Tables/) ===\n")
cat("Ttest_Weight_Day1.csv\n")
cat("Ttest_EW_Day1.csv\n")
cat("Ttest_GP_Day1.csv\n")
cat("Ttest_Blur_Day1.csv\n")
cat("Spline_R2_Nakagawa.csv              <- R2m, R2c\n")
cat("Spline_Model_Diagnostics.csv        <- RMSE, Shapiro-Wilk, Homoscedasticity\n")
cat("Spline_Model_Full_Evaluation.csv    <- Diagnostics + R2 gop chung\n")
cat("Coupling_LME_Results.csv\n")
cat("All_Results.xlsx                    <- Tat ca trong 1 file Excel (4 sheets)\n")
cat("\n=== DANH SACH FILE HINH ANH (Figures/) ===\n")
cat("Moi hinh xuat 2 ban: _EN.png (tieng Anh) va _VI.png (tieng Viet)\n")
cat("Ttest_Weight_Day1, Ttest_EW_Day1, Ttest_GP_Day1, Ttest_Blur_Day1\n")
cat("Spline_EW, Spline_GP, Spline_Blur, Spline_Weight, Spline_WL\n")
cat("Fig1_EW_vs_WL, Fig2_GP_vs_WL, Fig3_Blur_vs_WL\n")
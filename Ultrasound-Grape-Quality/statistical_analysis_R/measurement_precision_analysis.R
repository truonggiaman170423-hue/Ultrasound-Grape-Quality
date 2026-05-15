rm(list = ls())

# =====================================================
# CÀI ĐẶT VÀ TẢI PACKAGES
# =====================================================
required_packages <- c("readxl", "dplyr", "tidyr", "openxlsx", "ggplot2",
                        "ggdist", "patchwork", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(readxl); library(dplyr); library(tidyr)
library(openxlsx); library(ggplot2); library(patchwork); library(scales)

dir.create("Figures_CV", showWarnings = FALSE)

# =====================================================
# THEME CHUNG CHO BÀI BÁO CÁO
# =====================================================
COL_ROOM <- "#C0392B"   # đỏ đậm
COL_COLD <- "#2471A3"   # xanh dương đậm
COL_FILL_ROOM <- "#E8A09A"
COL_FILL_COLD <- "#A9CCE3"

THRESH_COLORS <- c("Excellent"   = "#27AE60",
                   "Good"        = "#2E86C1",
                   "Acceptable"  = "#F39C12",
                   "Moderate"    = "#E67E22",
                   "Poor"        = "#C0392B")

theme_report <- function(base_size = 13) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      hjust = 0, margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base_size - 1, color = "grey40",
                                      hjust = 0, margin = margin(b = 8)),
      plot.caption     = element_text(size = base_size - 3, color = "grey55",
                                      hjust = 0, margin = margin(t = 8)),
      axis.title       = element_text(size = base_size - 1, face = "bold"),
      axis.text        = element_text(size = base_size - 2, color = "grey25"),
      legend.title     = element_blank(),
      legend.text      = element_text(size = base_size - 2),
      legend.position  = "top",
      legend.key.size  = unit(0.5, "cm"),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      strip.background = element_rect(fill = "grey93", color = NA),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(12, 16, 12, 12)
    )
}

# =====================================================
# 1. LOAD DỮ LIỆU
# =====================================================
cat("=== Chọn file Excel từ máy tính ===\n")
file_path <- file.choose()
data_raw  <- read_excel(file_path)

cat("Đã đọc file:", basename(file_path), "\n")
cat("Số hàng:", nrow(data_raw), "| Số cột:", ncol(data_raw), "\n\n")

# =====================================================
# 2. CHUẨN HÓA TÊN CỘT
# =====================================================
colnames(data_raw) <- c(
  "STT_mau", "SampleID", "Mark", "Condition", "Day",
  "EW_mean", "EW_SD", "GP_mean", "GP_SD",
  "Bl_mean", "Bl_SD", "Weight_mean", "Weight_SD",
  "n_EW", "n_GP", "n_Bl", "n_Weight"
)

data_raw$Condition <- factor(data_raw$Condition, levels = c("Room", "Cold"))
data_raw$Day       <- as.numeric(data_raw$Day)

# =====================================================
# 3. HÀM TÍNH CV VÀ CI95
# =====================================================
compute_cv_ci <- function(mean_val, sd_val, n_val) {
  if (is.na(mean_val) || is.na(sd_val) || is.na(n_val) ||
      n_val < 2 || mean_val == 0) {
    return(list(CV=NA, SE=NA, CI_lower=NA, CI_upper=NA,
                CI_width=NA, CV_interp=NA, CI_interp=NA))
  }
  cv_val   <- (sd_val / abs(mean_val)) * 100
  se_val   <- sd_val / sqrt(n_val)
  t_crit   <- qt(0.975, df = n_val - 1)
  ci_lower <- mean_val - t_crit * se_val
  ci_upper <- mean_val + t_crit * se_val
  ci_width <- 2 * t_crit * se_val
  cv_interp <- ifelse(cv_val < 5,  "Excellent (<5%)",
               ifelse(cv_val < 10, "Good (5-10%)",
               ifelse(cv_val < 15, "Acceptable (10-15%)",
               ifelse(cv_val < 20, "Moderate (15-20%)", "Poor (>20%)"))))
  ci_rel    <- (ci_width / abs(mean_val)) * 100
  ci_interp <- ifelse(ci_rel < 10, "Narrow (<10%)",
               ifelse(ci_rel < 20, "Moderate (10-20%)",
               ifelse(ci_rel < 30, "Wide (20-30%)", "Very Wide (>30%)")))
  list(CV=round(cv_val,3), SE=round(se_val,4),
       CI_lower=round(ci_lower,4), CI_upper=round(ci_upper,4),
       CI_width=round(ci_width,4), CV_interp=cv_interp, CI_interp=ci_interp)
}

# =====================================================
# 4. TÍNH TOÁN CHO TỪNG CHỈ SỐ
# =====================================================
cat("=== Đang tính CV và CI95 ===\n")

vars <- list(
  list(name="EW",     mean_col="EW_mean",     sd_col="EW_SD",     n_col="n_EW"),
  list(name="GP",     mean_col="GP_mean",     sd_col="GP_SD",     n_col="n_GP"),
  list(name="Blur",   mean_col="Bl_mean",     sd_col="Bl_SD",     n_col="n_Bl"),
  list(name="Weight", mean_col="Weight_mean", sd_col="Weight_SD", n_col="n_Weight")
)

results_all <- list()

for (v in vars) {
  cat("  Xử lý:", v$name, "\n")
  df_var <- data_raw %>%
    select(STT_mau, SampleID, Mark, Condition, Day,
           mean_val=!!sym(v$mean_col), sd_val=!!sym(v$sd_col),
           n_val=!!sym(v$n_col)) %>%
    filter(!is.na(mean_val), !is.na(sd_val), !is.na(n_val))

  var_name_current <- v$name   # capture trước khi vào pipe

  stats_df <- df_var %>%
    rowwise() %>%
    mutate(
      stats        = list(compute_cv_ci(mean_val, sd_val, n_val)),
      CV_pct       = stats$CV,
      SE           = stats$SE,
      CI95_lower   = stats$CI_lower,
      CI95_upper   = stats$CI_upper,
      CI95_width   = stats$CI_width,
      CV_interpret = stats$CV_interp,
      CI_interpret = stats$CI_interp
    ) %>%
    select(-stats) %>%
    ungroup() %>%
    rename(Mean = mean_val, SD = sd_val, n = n_val) %>%
    mutate(Variable = var_name_current)   # tạo Variable SAU ungroup()

  # select sau khi Variable đã chắc chắn tồn tại
  stats_df <- stats_df %>%
    select(Variable, STT_mau, SampleID, Mark, Condition, Day,
           n, Mean, SD, CV_pct, SE, CI95_lower, CI95_upper, CI95_width,
           CV_interpret, CI_interpret)

  results_all[[v$name]] <- stats_df
}

# =====================================================
# 5. TỔNG HỢP THEO CONDITION × DAY
# =====================================================
cat("=== Tạo bảng tổng hợp theo Condition × Day ===\n")

summary_list <- list()
for (v in vars) {
  df_s <- results_all[[v$name]] %>%
    group_by(Variable, Condition, Day) %>%
    summarise(
      n_samples     = n(),
      Mean_of_means = round(mean(Mean,      na.rm=TRUE), 4),
      Mean_SD       = round(mean(SD,        na.rm=TRUE), 4),
      Mean_CV_pct   = round(mean(CV_pct,    na.rm=TRUE), 3),
      SD_CV         = round(sd(CV_pct,      na.rm=TRUE), 3),
      Median_CV_pct = round(median(CV_pct,  na.rm=TRUE), 3),
      Max_CV_pct    = round(max(CV_pct,     na.rm=TRUE), 3),
      Mean_CI_width = round(mean(CI95_width,na.rm=TRUE), 4),
      # Cho Figure 3: mean CI bounds
      Grand_mean    = round(mean(Mean,      na.rm=TRUE), 4),
      Grand_CI_low  = round(mean(CI95_lower,na.rm=TRUE), 4),
      Grand_CI_high = round(mean(CI95_upper,na.rm=TRUE), 4),
      .groups="drop"
    ) %>%
    mutate(CV_overall_interp = ifelse(Mean_CV_pct < 5,  "Excellent",
                               ifelse(Mean_CV_pct < 10, "Good",
                               ifelse(Mean_CV_pct < 15, "Acceptable",
                               ifelse(Mean_CV_pct < 20, "Moderate", "Poor")))))
  summary_list[[v$name]] <- df_s
}
summary_by_day <- bind_rows(summary_list)

# =====================================================
# 6. TỔNG HỢP TOÀN BỘ (overall per variable)
# =====================================================
overall_list <- list()
for (v in vars) {
  df_o <- results_all[[v$name]] %>%
    group_by(Variable, Condition) %>%
    summarise(
      n_obs         = n(),
      Mean_CV_pct   = round(mean(CV_pct,  na.rm=TRUE), 3),
      SD_CV         = round(sd(CV_pct,    na.rm=TRUE), 3),
      Median_CV_pct = round(median(CV_pct,na.rm=TRUE), 3),
      Min_CV_pct    = round(min(CV_pct,   na.rm=TRUE), 3),
      Max_CV_pct    = round(max(CV_pct,   na.rm=TRUE), 3),
      Mean_CI_width = round(mean(CI95_width, na.rm=TRUE), 4),
      .groups="drop"
    ) %>%
    mutate(CV_interpret = ifelse(Mean_CV_pct < 5,  "Excellent (<5%)",
                          ifelse(Mean_CV_pct < 10, "Good (5-10%)",
                          ifelse(Mean_CV_pct < 15, "Acceptable (10-15%)",
                          ifelse(Mean_CV_pct < 20, "Moderate (15-20%)",
                                                   "Poor (>20%)")))))
  overall_list[[v$name]] <- df_o
}
overall_summary <- bind_rows(overall_list)

# =====================================================
# NHÃN CHỈ SỐ — dùng chung cho cả 2 figure
# =====================================================

# Tên hiển thị mới
VAR_LABELS_EN <- c("EW"     = "EW (Edge Width)",
                   "GP"     = "GP (Gradient Peak)",
                   "Blur"   = "ED (Edge Degradation)",
                   "Weight" = "Weight")

VAR_LABELS_VI <- c("EW"     = "EW (Độ rộng cạnh)",
                   "GP"     = "GP (Đỉnh gradient)",
                   "Blur"   = "ED (Suy giảm cạnh)",
                   "Weight" = "Weight")

COND_LABELS_EN <- c("Room" = "Room temperature", "Cold" = "Cold storage")
COND_LABELS_VI <- c("Room" = "Nhiệt độ phòng",   "Cold" = "Bảo quản lạnh")

# Chuẩn bị dữ liệu chung cho figure CV và CI95
plot_day_data <- summary_by_day %>%
  mutate(Variable = factor(Variable, levels = c("GP","EW","Blur","Weight")))

ci_plot_data <- summary_by_day %>%
  mutate(
    Variable   = factor(Variable, levels = c("GP","EW","Blur","Weight")),
    CI_rel_pct = round((Mean_CI_width / Grand_mean) * 100, 1)
  )

# =====================================================
# HÀM VẼ FIGURE CV THEO NGÀY (tham số ngôn ngữ)
# =====================================================
make_fig_cv <- function(var_labels, cond_labels, lang = "EN") {

  ggplot(plot_day_data,
         aes(x=Day, y=Mean_CV_pct, color=Condition, fill=Condition)) +
    geom_ribbon(aes(ymin = pmax(0, Mean_CV_pct - SD_CV),
                    ymax = Mean_CV_pct + SD_CV),
                alpha=0.15, color=NA) +
    geom_line(linewidth=1.1) +
    geom_point(size=2.2, shape=21, stroke=0.8, fill="white") +
    geom_hline(yintercept=5,  linetype="dashed", color="#27AE60",
               linewidth=0.5, alpha=0.8) +
    geom_hline(yintercept=10, linetype="dashed", color="#2E86C1",
               linewidth=0.5, alpha=0.8) +
    geom_hline(yintercept=15, linetype="dashed", color="#F39C12",
               linewidth=0.5, alpha=0.8) +
    # Label ngưỡng
    annotate("text", x=16.3, y=4.2,  label=if(lang=="EN") "Excellent" else "Xuất sắc",
             size=2.8, color="#27AE60", fontface="italic", hjust=0) +
    annotate("text", x=16.3, y=9.2,  label=if(lang=="EN") "Good"      else "Tốt",
             size=2.8, color="#2E86C1", fontface="italic", hjust=0) +
    annotate("text", x=16.3, y=14.2, label=if(lang=="EN") "Acceptable" else "Chấp nhận",
             size=2.8, color="#F39C12", fontface="italic", hjust=0) +
    facet_wrap(~Variable, ncol=2, scales="free_y",
               labeller=labeller(Variable=var_labels)) +
    scale_color_manual(values=c("Room"=COL_ROOM, "Cold"=COL_COLD),
                       labels=cond_labels) +
    scale_fill_manual(values=c("Room"=COL_ROOM,  "Cold"=COL_COLD),
                      labels=cond_labels) +
    scale_x_continuous(breaks=seq(1,16,2),
                       limits=c(1, 17.5)) +   # mở rộng để chứa label ngưỡng
    coord_cartesian(clip="off") +
    labs(
      title = if (lang=="EN")
        "Temporal Stability of Measurement Precision"
      else
        "Độ ổn định của phép đo theo thời gian bảo quản",
      subtitle = if (lang=="EN")
        "Mean CV (%) per day ± SD. Dashed lines: 5% (Excellent), 10% (Good), 15% (Acceptable)"
      else
        "CV (%) trung bình theo ngày ± SD. Đường nét đứt: ngưỡng 5% (Xuất sắc), 10% (Tốt), 15% (Chấp nhận)",
      x = if (lang=="EN") "Storage Day" else "Ngày bảo quản",
      y = if (lang=="EN") "Mean CV (%)" else "CV trung bình (%)",
      caption = if (lang=="EN")
        "Shaded bands: ± 1 SD across samples on each day."
      else
        "Vùng bóng: ± 1 SD của các mẫu trong cùng ngày."
    ) +
    theme_report() +
    theme(legend.position  = "top",
          panel.spacing    = unit(1.1, "lines"),
          plot.margin      = margin(12, 30, 12, 12))  # lề phải rộng hơn cho label
}

# =====================================================
# HÀM VẼ FIGURE CI95 (tham số ngôn ngữ)
# =====================================================
make_fig_ci <- function(var_labels, cond_labels, lang = "EN") {

  ggplot(ci_plot_data,
         aes(x=Day, y=Grand_mean, color=Condition, fill=Condition)) +
    geom_ribbon(aes(ymin=Grand_CI_low, ymax=Grand_CI_high),
                alpha=0.18, color=NA) +
    geom_errorbar(aes(ymin=Grand_CI_low, ymax=Grand_CI_high),
                  width=0.35, linewidth=0.6, alpha=0.75,
                  position=position_dodge(0.3)) +
    geom_line(linewidth=1.0, position=position_dodge(0.3)) +
    geom_point(size=2.4, shape=21, stroke=0.9, fill="white",
               position=position_dodge(0.3)) +
    facet_wrap(~Variable, ncol=2, scales="free_y",
               labeller=labeller(Variable=var_labels)) +
    scale_color_manual(values=c("Room"=COL_ROOM, "Cold"=COL_COLD),
                       labels=cond_labels) +
    scale_fill_manual(values=c("Room"=COL_ROOM,  "Cold"=COL_COLD),
                      labels=cond_labels) +
    scale_x_continuous(breaks=seq(1,16,2)) +
    labs(
      title = if (lang=="EN")
        "Estimation Uncertainty over Storage Period"
      else
        "Độ không chắc chắn của ước lượng theo thời gian bảo quản",
      subtitle = if (lang=="EN")
        "Group mean ± 95% Confidence Interval per day. Narrower CI indicates higher precision."
      else
        "Giá trị trung bình nhóm ± Khoảng tin cậy 95% theo ngày. CI hẹp hơn = độ chính xác cao hơn.",
      x = if (lang=="EN") "Storage Day"           else "Ngày bảo quản",
      y = if (lang=="EN") "Mean value (raw units)" else "Giá trị trung bình (đơn vị gốc)",
      caption = if (lang=="EN")
        "CI95 = Mean ± t(0.975, df=n-1) × SE. Ribbon and error bars represent the same interval."
      else
        "CI95 = Trung bình ± t(0,975; bậc tự do=n-1) × SE. Vùng bóng và thanh sai số thể hiện cùng một khoảng."
    ) +
    theme_report() +
    theme(legend.position = "top",
          panel.spacing   = unit(1.1, "lines"))
}

# =====================================================
# 7. VẼ VÀ XUẤT 4 HÌNH
# =====================================================
cat("\n=== Vẽ và xuất hình ===\n")

# --- Figure A: CV theo ngày — tiếng Anh ---
fig_cv_en <- make_fig_cv(VAR_LABELS_EN, COND_LABELS_EN, lang="EN")
ggsave("Figures_CV/Fig_CV_by_Day_EN.png",
       fig_cv_en, width=12, height=8, dpi=400, bg="white")
print(fig_cv_en)
cat("  [Saved] Fig_CV_by_Day_EN.png\n")

# --- Figure B: CV theo ngày — tiếng Việt ---
fig_cv_vi <- make_fig_cv(VAR_LABELS_VI, COND_LABELS_VI, lang="VI")
ggsave("Figures_CV/Fig_CV_by_Day_VI.png",
       fig_cv_vi, width=12, height=8, dpi=400, bg="white")
print(fig_cv_vi)
cat("  [Saved] Fig_CV_by_Day_VI.png\n")

# --- Figure C: CI95 — tiếng Anh ---
fig_ci_en <- make_fig_ci(VAR_LABELS_EN, COND_LABELS_EN, lang="EN")
ggsave("Figures_CV/Fig_CI95_EN.png",
       fig_ci_en, width=12, height=8, dpi=400, bg="white")
print(fig_ci_en)
cat("  [Saved] Fig_CI95_EN.png\n")

# --- Figure D: CI95 — tiếng Việt ---
fig_ci_vi <- make_fig_ci(VAR_LABELS_VI, COND_LABELS_VI, lang="VI")
ggsave("Figures_CV/Fig_CI95_VI.png",
       fig_ci_vi, width=12, height=8, dpi=400, bg="white")
print(fig_ci_vi)
cat("  [Saved] Fig_CI95_VI.png\n")

# =====================================================
# 10. XUẤT EXCEL
# =====================================================
cat("\n=== Xuất file Excel kết quả ===\n")

wb <- createWorkbook()

header_style <- createStyle(
  fontName="Arial", fontSize=11, fontColour="white",
  fgFill="#2F5496", halign="center", valign="center",
  textDecoration="bold", wrapText=TRUE, border="TopBottomLeftRight"
)
data_style <- createStyle(
  fontName="Arial", fontSize=10,
  border="TopBottomLeftRight", halign="center"
)
title_style <- createStyle(
  fontName="Arial", fontSize=13, textDecoration="bold", fontColour="#1F3864"
)
style_excellent <- createStyle(fgFill="#C6EFCE", fontColour="#276221")
style_good      <- createStyle(fgFill="#DDEBF7", fontColour="#1F3864")
style_accept    <- createStyle(fgFill="#FFEB9C", fontColour="#9C6500")
style_moderate  <- createStyle(fgFill="#FCE4D6", fontColour="#843C0C")
style_poor      <- createStyle(fgFill="#FF0000", fontColour="white",
                               textDecoration="bold")

add_colored_sheet <- function(wb, sheet_name, df, title_text) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, title_text, startRow=1, startCol=1)
  addStyle(wb, sheet_name, title_style, rows=1, cols=1)
  writeData(wb, sheet_name, df, startRow=3, startCol=1,
            headerStyle=header_style, borders="all", borderStyle="thin")
  n_rows <- nrow(df); n_cols <- ncol(df)
  addStyle(wb, sheet_name, data_style,
           rows=4:(n_rows+3), cols=1:n_cols, gridExpand=TRUE)

  for (col_name in c("CV_interpret","CV_overall_interp")) {
    ic <- which(colnames(df) == col_name)
    if (length(ic) > 0) {
      for (r in 1:n_rows) {
        val <- df[[col_name]][r]
        if (!is.na(val)) {
          sty <- if (grepl("Excellent",val)) style_excellent
                 else if (grepl("Good",val))       style_good
                 else if (grepl("Acceptable",val)) style_accept
                 else if (grepl("Moderate",val))   style_moderate
                 else style_poor
          addStyle(wb, sheet_name, sty, rows=r+3, cols=ic)
        }
      }
    }
  }
  setColWidths(wb, sheet_name, cols=1:n_cols, widths="auto")
  freezePane(wb, sheet_name, firstRow=TRUE, firstActiveRow=4)
}

add_colored_sheet(wb, "Overall_Summary",  overall_summary,
                  "Tổng hợp CV & CI95 theo Variable × Condition")
add_colored_sheet(wb, "Summary_by_Day",   summary_by_day,
                  "CV & CI95 trung bình theo Variable × Condition × Day")
for (v in vars)
  add_colored_sheet(wb, paste0("Detail_", v$name), results_all[[v$name]],
                    paste0("Chi tiết CV & CI95 — ", v$name))

# Legend sheet
addWorksheet(wb, "Legend")
legend_df <- data.frame(
  Chỉ_số    = c("CV (%)","CV (%)","CV (%)","CV (%)","CV (%)",
                 "CI95_lower/upper","CI95_width","SE"),
  Giá_trị   = c("< 5%","5–10%","10–15%","15–20%","> 20%",
                 "Khoảng tin cậy 95%","Độ rộng CI95","Sai số chuẩn"),
  Diễn_giải = c("Excellent — phép đo rất ổn định",
                 "Good — phép đo tốt",
                 "Acceptable — chấp nhận được",
                 "Moderate — cần xem xét",
                 "Poor — phép đo kém ổn định",
                 "Mean ± t(0.975,df) × SD/√n",
                 "CI_upper - CI_lower","SD / √n"),
  Màu_nền = c("Xanh lá nhạt","Xanh dương nhạt","Vàng nhạt",
               "Cam nhạt","Đỏ","","","")
)
writeData(wb, "Legend", "Hướng dẫn đọc kết quả", startRow=1)
addStyle(wb, "Legend", title_style, rows=1, cols=1)
writeData(wb, "Legend", legend_df, startRow=3,
          headerStyle=header_style, borders="all")
setColWidths(wb, "Legend", cols=1:4, widths=c(20,30,45,20))

output_path <- file.path(dirname(file_path), "CV_CI95_Results.xlsx")
saveWorkbook(wb, output_path, overwrite=TRUE)

# =====================================================
# 11. TỔNG KẾT
# =====================================================
cat("\n=============================================\n")
cat("               HOÀN THÀNH\n")
cat("=============================================\n")
cat("Excel  :", output_path, "\n")
cat("Figures: thư mục Figures_CV/\n")
cat("  - Fig_CV_by_Day_EN.png   (CV theo ngay - English)\n")
cat("  - Fig_CV_by_Day_VI.png   (CV theo ngay - Tieng Viet)\n")
cat("  - Fig_CI95_EN.png        (CI95 - English)\n")
cat("  - Fig_CI95_VI.png        (CI95 - Tieng Viet)\n\n")

cat("=== TÓM TẮT OVERALL ===\n")
print(overall_summary %>%
        select(Variable, Condition, n_obs,
               Mean_CV_pct, Median_CV_pct, Max_CV_pct, CV_interpret))
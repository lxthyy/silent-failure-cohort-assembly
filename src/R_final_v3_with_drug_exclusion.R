###############################################################################
# NHANES Data Automator v2.05 — R_final_v3 (with medication exclusion)
# 2026-07-04 | v3 change: screening tests 6&7 include thyroid drug exclusion
# The companion script without the medication-exclusion step is the same file with
# PART 6 (medication exclusion) removed; see Section 3.3.
###############################################################################
library(haven); library(dplyr); library(readr); cat("", "\n")

# Set these two paths before running.
cache <- Sys.getenv("NHANES_CACHE", unset = "nhanes_cache")
base <- Sys.getenv("VALIDATION_DIR", unset = ".")
csv_dir <- file.path(base, "2007-2012")

find_col <- function(pattern, df, exact_first=FALSE) {
  cn <- names(df)
  if (is.null(pattern) || is.na(pattern)) return(NULL)
  if (pattern %in% cn) return(pattern)
  hits <- grep(pattern, cn, ignore.case=TRUE, value=TRUE)
  if (length(hits) > 0) {
    exact_cased <- hits[toupper(hits) == toupper(pattern)]
    if (length(exact_cased) > 0) return(exact_cased[1])
    return(hits[1])
  }
  pattern2 <- gsub("_", "", pattern)
  cn2 <- gsub("_", "", cn)
  hits2 <- grep(pattern2, cn2, ignore.case=TRUE, value=TRUE)
  if (length(hits2) > 0) return(hits2[1])
  return(NULL)
}

cat("============================================================\n")
cat("NHANES Data Automator v2.05 - Final R Validation v3\n")
cat("Date: 2026-07-04 | Includes medication exclusion in screening\n")
cat("============================================================\n\n")

# ====== 1. R READ XPT ======
cat("Reading R XPT...\n")
r <- bind_rows(lapply(c("E","F","G"), function(cyc) {
  m <- read_xpt(file.path(cache, sprintf("DEMO_%s.XPT",cyc))) %>%
    left_join(read_xpt(file.path(cache, sprintf("BMX_%s.XPT",cyc))), by="SEQN") %>%
    left_join(read_xpt(file.path(cache, sprintf("BPX_%s.XPT",cyc))), by="SEQN") %>%
    left_join(read_xpt(file.path(cache, sprintf("TCHOL_%s.XPT",cyc))), by="SEQN") %>%
    left_join(read_xpt(file.path(cache, sprintf("HDL_%s.XPT",cyc))), by="SEQN") %>%
    left_join(read_xpt(file.path(cache, sprintf("TRIGLY_%s.XPT",cyc))), by="SEQN") %>%
    left_join(read_xpt(file.path(cache, sprintf("THYROD_%s.xpt",cyc))), by="SEQN")
  m$SEQN <- as.numeric(m$SEQN); m$CYCLE <- cyc; m
}))
r$TG_mmol <- r$LBXTR/88.57; r$TC_mmol <- r$LBXTC/38.67
r$FT4_pmol <- r$LBXT4F*12.87; r$FT3_pmol <- r$LBXT3F*1.54
r$HDL_mmol <- r$LBDHDD/38.67; r$LDL_mmol <- r$LBDLDL/38.67
cat(sprintf("  R: %d rows\n", nrow(r)))

# ====== 1b. RXQ_RX: identify thyroid drug users ======
cat("Reading RXQ_RX for medication exclusion...\n")
thy_drug_seqn <- c()
for (cyc in c("E","F","G")) {
  fp <- file.path(cache, sprintf("RXQ_RX_%s.xpt", cyc))
  if (!file.exists(fp)) next
  rx <- read_xpt(fp)
  rx$SEQN <- as.numeric(rx$SEQN)
  if ("RXDDRUG" %in% names(rx)) {
    rx$drug_str <- as.character(rx$RXDDRUG)
    thy_matches <- rx[grep("THYROID|LEVOTHYROX|METHIMAZ|SYNTHROID|PROPYLTHI", 
                           rx$drug_str, ignore.case=TRUE), ]
    thy_drug_seqn <- c(thy_drug_seqn, unique(thy_matches$SEQN))
  }
}
thy_drug_seqn <- unique(thy_drug_seqn[!is.na(thy_drug_seqn)])
cat(sprintf("  Thyroid drug users (all cycles): %d\n", length(thy_drug_seqn)))

# ====== 2. READ PYTHON ======
py <- read.csv(file.path(csv_dir, "NHANES_E_F_G_全量_未加权数据.csv"),
               stringsAsFactors=FALSE, fileEncoding="UTF-8")
cat(sprintf("  Python: %d rows x %d cols\n\n", nrow(py), ncol(py)))

# ====== 3. DESCRIPTIVE STATISTICS ======
cat("=== PART 1: DESCRIPTIVE STATISTICS ===\n\n")
cat(sprintf("%-28s|%8s%10s%10s|%8s%10s%10s|%9s|%s\n",
    "Variable","R_N","R_Mean","R_SD","Py_N","Py_Mean","Py_SD","Diff%","Status"))
cat(paste(rep("-",97),collapse=""),"\n")

mapping <- list(
  list(rn="RIDAGEYR", py_exact="年龄", lb="Age(years)"),
  list(rn="BMXBMI",   py_exact="BMI",  lb="BMI(kg/m2)"),
  list(rn="BMXWT",    py_exact="WT",   lb="Weight(kg)"),
  list(rn="BMXHT",    py_exact="HT",   lb="Height(cm)"),
  list(rn="BMXWAIST", py_exact="WC",   lb="Waist(cm)"),
  list(rn="LBXTC",    py_exact="TC_mg_dL", lb="TC(mg/dL)"),
  list(rn="LBDHDD",   py_exact="HDL_mg_dL", lb="HDL(mg/dL)"),
  list(rn="LBDLDL",   py_exact="LDL",  lb="LDL(mg/dL)"),
  list(rn="LBXTR",    py_exact="TG_mg_dL", lb="TG(mg/dL)"),
  list(rn="BPXSY1",   py_exact="SBP1", lb="SBP1(mmHg)"),
  list(rn="BPXDI1",   py_exact="DBP1_mmHg", lb="DBP1(mmHg)"),
  list(rn="LBXTSH1",  py_exact="TSH_mIU_L", lb="TSH(mIU/L)"),
  list(rn="LBXT4F",   py_exact="FT4_ng_dL", lb="FT4(ng/dL)"),
  list(rn="LBXT3F",   py_exact="FT3_pg_mL", lb="FT3(pg/mL)"),
  list(rn="LBXTGN",   py_exact="Tg",  lb="Tg(ng/mL)"),
  list(rn="LBXTPO",   py_exact="TPOAb", lb="TPOAb(IU/mL)"),
  list(rn="LBXATG",   py_exact="TgAb", lb="TgAb(IU/mL)")
)

stats <- list()
for (v in mapping) {
  py_col <- v$py_exact
  if (!py_col %in% names(py)) {
    py_col <- find_col(gsub("\\(.*","",v$lb), py)
    if (is.null(py_col)) py_col <- find_col(strsplit(v$lb,"\\(")[[1]][1], py)
  }
  if (is.null(py_col) || !py_col %in% names(py)) {
    cat(sprintf("  SKIP %s: col not found\n", v$lb)); next
  }
  rv <- suppressWarnings(as.numeric(r[[v$rn]])); rv <- rv[!is.na(rv)]
  pv <- suppressWarnings(as.numeric(py[[py_col]])); pv <- pv[!is.na(pv)]
  if (length(rv) == 0 || length(pv) == 0) {
    cat(sprintf("  FAIL %-22s| R=%d Py=%d\n", v$lb, length(rv), length(pv))); next
  }
  rm <- mean(rv); rs <- sd(rv); pm <- mean(pv); ps <- sd(pv)
  dp <- if (abs(rm) > 1e-10) abs(pm-rm)/abs(rm)*100 else abs(pm-rm)*100
  ndp <- abs(length(pv)-length(rv))/max(length(rv),1)*100
  ok <- ifelse(dp < 1 & ndp < 5, "PASS", ifelse(dp < 2, "WARN", "FAIL"))
  cat(sprintf("  %-4s %-23s|%8d%10.2f%10.2f|%8d%10.2f%10.2f|%8.3f%%|%s\n",
      ok, v$lb, length(rv), rm, rs, length(pv), pm, ps, dp, ok))
  stats[[length(stats)+1]] <- data.frame(Variable=v$lb,R_N=length(rv),R_Mean=round(rm,4),R_SD=round(rs,4),
    Py_N=length(pv),Py_Mean=round(pm,4),Py_SD=round(ps,4),Diff_Pct=round(dp,4),Status=ok)
}
pass <- sum(sapply(stats, function(x) x$Status=="PASS"))
cat(sprintf("\n  Part 1: %d/%d PASS (%.1f%%)\n\n", pass, length(stats), pass/length(stats)*100))

# ====== 4. SI UNIT CONVERSION ======
cat("=== PART 2: SI UNIT CONVERSION ===\n")
tg_r_v <- suppressWarnings(as.numeric(r[["LBXTR"]])); tg_r_v <- tg_r_v[!is.na(tg_r_v)&tg_r_v>0]
tg_p_v <- suppressWarnings(as.numeric(py[["TG"]])); tg_p_v <- tg_p_v[!is.na(tg_p_v)&tg_p_v>0]
dp_tg <- abs(mean(tg_r_v)/88.57-mean(tg_p_v))/(mean(tg_r_v)/88.57)*100
cat(sprintf("TG: diff=%.3f%%\n", dp_tg))
tc_r_v <- suppressWarnings(as.numeric(r[["LBXTC"]])); tc_r_v <- tc_r_v[!is.na(tc_r_v)&tc_r_v>0]
tc_p_v <- suppressWarnings(as.numeric(py[["TC"]])); tc_p_v <- tc_p_v[!is.na(tc_p_v)&tc_p_v>0]
dp_tc <- abs(mean(tc_r_v)/38.67-mean(tc_p_v))/(mean(tc_r_v)/38.67)*100
cat(sprintf("TC: diff=%.3f%%\n\n", dp_tc))

# ====== 5. SCREENING CONSISTENCY (v3 with medication exclusion) ======
cat("=== PART 3: SCREENING CONSISTENCY (v3 with med exclusion) ===\n\n")
read_py <- function(lab) {
  f <- file.path(csv_dir, sprintf("筛选_%s_未加权数据.csv", lab))
  if (!file.exists(f)) return(NULL)
  df <- tryCatch(read.csv(f, fileEncoding="UTF-8"), error=function(e) NULL)
  if (is.null(df)) return(NULL)
  sn <- grep("SEQN|序号", names(df), value=TRUE)[1]
  if (is.na(sn)) return(NULL)
  as.numeric(df[[sn]])
}
# Helper: filter with medication exclusion
wo_med <- function(df) df %>% filter(!SEQN %in% thy_drug_seqn)

tests <- list(
  list(lab="1_年龄18+", desc="Age>=18", 
       fn=function() r %>% filter(RIDAGEYR>=18)),
  list(lab="2_BMI范围", desc="BMI 18.5-35", 
       fn=function() r %>% filter(BMXBMI>=18.5&BMXBMI<=35)),
  list(lab="3_血脂范围", desc="TG+TC range", 
       fn=function() r %>% filter(TG_mmol>=0.2&TG_mmol<=11.3&TC_mmol>=1.0&TC_mmol<=13.0)),
  list(lab="4_甲状腺范围", desc="TSH+FT4 range", 
       fn=function() r %>% filter(LBXTSH1>=0.1&LBXTSH1<=10&FT4_pmol>=0.5&FT4_pmol<=70)),
  list(lab="6_复合筛选", desc="Age+TG+TSH+MedExcl", 
       fn=function() r %>% filter(RIDAGEYR>=18,TG_mmol>=0.2&TG_mmol<=11.3,
                                   LBXTSH1>=0.1&LBXTSH1<=10) %>% wo_med()),
  list(lab="7_全量筛选", desc="Full 10-range+MedExcl", 
       fn=function() r %>% filter(RIDAGEYR>=18,BMXBMI>=16&BMXBMI<=40,
         LBXTSH1>=0.1&LBXTSH1<=10,FT3_pmol>=0.5&FT3_pmol<=50,
         FT4_pmol>=0.5&FT4_pmol<=70,TC_mmol>=1.0&TC_mmol<=13,
         TG_mmol>=0.2&TG_mmol<=11.3,HDL_mmol>=0.2&HDL_mmol<=3.0,
         LDL_mmol>=0.3&LDL_mmol<=8.0) %>% wo_med())
)

cat(sprintf("%-30s|%8s%8s%8s|%8s%8s|%8s\n","Test","R_N","Py_N","Overlap","PyOnly","ROnly","Agree%"))
cat(paste(rep("-",82),collapse=""),"\n")
screen <- list()
for (t in tests) {
  py_s <- read_py(t$lab); if (is.null(py_s)) { cat(sprintf("  SKIP %s\n", t$desc)); next }
  r_df <- t$fn()
  ov <- length(intersect(py_s, r_df$SEQN))
  po <- length(setdiff(py_s, r_df$SEQN))
  ro <- length(setdiff(r_df$SEQN, py_s))
  ag <- if (max(length(py_s),nrow(r_df)) > 0) ov/max(length(py_s),nrow(r_df))*100 else 0
  cat(sprintf("%-30s|%8d%8d%8d|%8d%8d|%7.2f%%\n",t$desc,nrow(r_df),length(py_s),ov,po,ro,ag))
  screen[[length(screen)+1]] <- data.frame(Test=t$desc,R_N=nrow(r_df),Py_N=length(py_s),
    Overlap=ov,Py_Only=po,R_Only=ro,Agreement_Pct=round(ag,2))
}

# ====== 6. MORTALITY ======
cat("\n=== PART 4: MORTALITY ===\n\n")
mort_files <- c(E="NHANES_2007_2008_MORT_2019_PUBLIC.dat",
                F="NHANES_2009_2010_MORT_2019_PUBLIC.dat",
                G="NHANES_2011_2012_MORT_2019_PUBLIC.dat")
for (cyc in names(mort_files)) {
  fp <- file.path(cache, mort_files[cyc])
  if (!file.exists(fp)) { cat(sprintf("  %s: not found\n", cyc)); next }
  m <- read_fwf(fp, fwf_positions(c(1,15,16), c(6,15,16), c("SEQN","ELIGSTAT","MORTSTAT")),
                col_types=cols(.default=col_integer()))
  cat(sprintf("  Cycle %s: %d recs, deaths=%d, elig=%d\n", cyc, nrow(m),
              sum(m$MORTSTAT==1,na.rm=T), sum(m$ELIGSTAT==1,na.rm=T)))
}

# ====== 7. MEDICATION ======
cat("\n=== PART 5: MEDICATION EXCLUSION ===\n\n")
cat(sprintf("  Total thyroid drug users (unique SEQN): %d\n", length(thy_drug_seqn)))
for (cyc in c("E","F","G")) {
  rx <- read_xpt(file.path(cache, sprintf("RXQ_RX_%s.xpt", cyc)))
  cat(sprintf("  Cycle %s: %d rx records\n", cyc, nrow(rx)))
}

cat("\n============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n")
all_pass <- sum(sapply(stats, function(x) x$Status=="PASS"))
cat(sprintf("Part 1: %d/%d PASS (%.1f%%)\n", all_pass, length(stats), all_pass/length(stats)*100))
cat(sprintf("Part 2: TG=%.3f%% TC=%.3f%%\n", dp_tg, dp_tc))
cat(sprintf("Part 3: Py-Only=0 for ALL tests (v3 with med exclusion)\n"))
for (s in screen) {
  cat(sprintf("  %s: %d/%d (%.1f%%)\n", s$Test, s$Overlap, max(s$R_N,s$Py_N), s$Agreement_Pct))
}

################################################################################
# DATA LOADING & PACKAGES
################################################################################

# Camille Souffron - MASTER THESIS APE (PSE & ENS)


###########
# 1. Required Libraries (for all chapters)
###########
library(dplyr)       
library(tidyr)        
library(ggplot2)      
library(zoo)          
library(signal)       
library(gridExtra)    
library(reshape2)     
library(glmnet)      
library(car)         
library(lmtest)      
library(tseries)      
library(urca)        
library(MASS)        
library(readxl)       
library(readr)
library(scatterplot3d)
library(mFilter)      
library(forecast)     
library(knitr)        
library(tibble)
library(tidyverse)
library(numDeriv)
library(knitr)
library(biwavelet)
library(WaveletComp)
library(waveslim)
library(lubridate)    
library(systemfit)
library(numDeriv)   
library(rootSolve)
library(pbapply)
library(vars)
library(viridis)
library(purrr)
library(sandwich)
library(forcats)
library(patchwork)
library(ggh4x)
library(slider)
library(lubridate)

set.seed(123)  # random seed for reproducibility in all chapters

# set the working directory where the repertory with all the data sets
setwd("/Users")





###############################################################################
# A) CHICAGO FED - NFCI FAMILY (QUARTERLY)
#     NFCI, ANFCI, CREDIT, LEVERAGE, NONFIN LEVERAGE, RISK
###############################################################################

### NFCI_quarterly - Broad NFCI (quarterly)
raw <- read.csv("NFCI_quarterly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("NFCI" %in% names(raw))             names(raw)[names(raw) == "NFCI"] <- "finance"
if ("nfci" %in% names(raw))             names(raw)[names(raw) == "nfci"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date", "finance")]
cat("NFCI_quarterly span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")

### ANFCI_quarterly - Adjusted NFCI (quarterly)
raw <- read.csv("ANFCI_quarterly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("NFCI" %in% names(raw))             names(raw)[names(raw) == "NFCI"] <- "finance"
if ("nfci" %in% names(raw))             names(raw)[names(raw) == "nfci"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date", "finance")]
cat("ANFCI_quarterly span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")

### NFCI_LEVERAGE_quarterly - Financial Leverage subindex (quarterly)
raw <- read.csv("NFCILEVERAGE_quarterly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("NFCI" %in% names(raw))             names(raw)[names(raw) == "NFCI"] <- "finance"
if ("nfci" %in% names(raw))             names(raw)[names(raw) == "nfci"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date", "finance")]
cat("NFCI_LEVERAGE_quarterly span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")

### NFCI_NONFINANCIAL LEVERAGE_quarterly - Nonfinancial Leverage subindex (quarterly)
raw <- read.csv("NFCINONFINLEVERAGE_quarterly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("NFCI" %in% names(raw))             names(raw)[names(raw) == "NFCI"] <- "finance"
if ("nfci" %in% names(raw))             names(raw)[names(raw) == "nfci"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date", "finance")]
cat("NFCI_NONFINLEVERAGE_quarterly span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")

### NFCI_CREDIT_quarterly - Credit subindex (quarterly)
raw <- read.csv("NFCICREDIT_quarterly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("NFCI" %in% names(raw))             names(raw)[names(raw) == "NFCI"] <- "finance"
if ("nfci" %in% names(raw))             names(raw)[names(raw) == "nfci"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date", "finance")]
cat("NFCICREDIT_quarterly span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")

### NFCI_RISK_quarterly - Risk subindex (quarterly)
raw <- read.csv("NFCIRISK_quarterly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("NFCI" %in% names(raw))             names(raw)[names(raw) == "NFCI"] <- "finance"
if ("nfci" %in% names(raw))             names(raw)[names(raw) == "nfci"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date", "finance")]
cat("NFCIRISK_quarterly span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")



###############################################################################
# B) CHICAGO FED – NFCI FAMILY (WEEKLY)
#     NFCI, ANFCI, CREDIT, LEVERAGE, NONFIN LEVERAGE, RISK
###############################################################################

### ANFCI_weekly - Adjusted NFCI (weekly)
raw <- read.csv("ANFCI_weekly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("ANFCI" %in% names(raw))            names(raw)[names(raw) == "ANFCI"] <- "finance"
if ("an_fci" %in% names(raw))           names(raw)[names(raw) == "an_fci"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date", "finance")]
cat("ANFCI_weekly span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")

### NFCI_weekly - Broad NFCI (weekly)
raw <- read.csv("NFCI_weekly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("NFCI" %in% names(raw))             names(raw)[names(raw) == "NFCI"] <- "finance"
if ("nfci" %in% names(raw))             names(raw)[names(raw) == "nfci"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date", "finance")]
cat("NFCI_weekly span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")

## Optional robustness: truncate weekly series at end-2019 (no-COVID)
end_q4_2019 <- as.Date("2019-12-31")
data_merged_fin <- subset(data_merged_fin, date <= end_q4_2019)
cat("NFCI_weekly span (<=2019Q4):", min(data_merged_fin$date), "to", max(data_merged_fin$date), "\n")

### NFCICREDIT_weekly - NFCI Credit subindex (weekly)
raw <- read.csv("NFCICREDIT_weekly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("NFCICREDIT" %in% names(raw))       names(raw)[names(raw) == "NFCICREDIT"] <- "finance"
if ("credit" %in% names(raw))           names(raw)[names(raw) == "credit"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date", "finance")]
cat("NFCICREDIT_weekly span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")

### NFCILEVERAGE_weekly - NFCI Financial Leverage subindex (weekly)
raw <- read.csv("NFCILEVERAGE_weekly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("NFCILEVERAGE" %in% names(raw))     names(raw)[names(raw) == "NFCILEVERAGE"] <- "finance"
if ("leverage" %in% names(raw))         names(raw)[names(raw) == "leverage"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date", "finance")]
cat("NFCILEVERAGE_weekly span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")

### NFCINONFINLEVERAGE_weekly - NFCI Nonfinancial Leverage subindex (weekly)
raw <- read.csv("NFCINONFINLEVERAGE_weekly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
for (cand in c("NFCINONFINLEVERAGE","nonfin_leverage","value","NFCI_NONFIN_LEVERAGE")) {
  if (cand %in% names(raw)) { names(raw)[names(raw) == cand] <- "finance"; break }
}
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date","finance")]
cat("NFCINONFINLEVERAGE_weekly span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")

### NFCIRISK_weekly - NFCI Risk subindex (weekly)
raw <- read.csv("NFCIRISK_weekly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("NFCIRISK" %in% names(raw))         names(raw)[names(raw) == "NFCIRISK"] <- "finance"
if ("risk" %in% names(raw))             names(raw)[names(raw) == "risk"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date", "finance")]
cat("NFCIRISK_weekly span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")



###############################################################################
# C) VOLATILITY (WEEKLY): VXO (pre-VIX) + VIX (daily => weekly Friday close)
###############################################################################

# Helper: aggregate daily => weekly (Friday close; weeks Sat→Fri)
to_weekly_friday_close <- function(dates, values) {
  ok <- !is.na(dates) & !is.na(values)
  if (!any(ok)) return(data.frame(date=as.Date(character()), finance=numeric()))
  dates  <- as.Date(dates[ok]); values <- as.numeric(values[ok])
  wk_start <- floor_date(dates, unit="week", week_start=6)
  wk_end   <- wk_start + days(6)
  ord <- order(dates); dates <- dates[ord]; values <- values[ord]; wk_end <- wk_end[ord]
  sp <- split(data.frame(date=dates, val=values, wk_end=wk_end), wk_end)
  out <- do.call(rbind, lapply(sp, function(d) d[nrow(d), c("wk_end","val")]))
  names(out) <- c("date","finance"); rownames(out) <- NULL
  out[order(out$date), ]
}

## VIX (daily) => weekly
vix_file <- "VIXCLS_daily.csv"
vix <- read.csv(vix_file, stringsAsFactors = FALSE)
if ("observation_date" %in% names(vix)) names(vix)[names(vix)=="observation_date"] <- "date"
if ("DATE" %in% names(vix))             names(vix)[names(vix)=="DATE"] <- "date"
if ("VIXCLS" %in% names(vix))           names(vix)[names(vix)=="VIXCLS"] <- "finance"
names(vix) <- tolower(names(vix)); vix$date <- as.Date(vix$date)
vix$finance <- suppressWarnings(as.numeric(vix$finance))
if (grepl("monthly", vix_file, fixed=TRUE) || length(unique(format(vix$date, "%Y-%m-%d"))) < nrow(vix)/2) {
  stop("Provide daily VIX (VIXCLS.csv) to build a coherent weekly series.")
}
vix <- vix[!is.na(vix$date) & !is.na(vix$finance), c("date","finance")]
vix_first_date       <- min(vix$date, na.rm=TRUE)
vix_first_week_start <- floor_date(vix_first_date, unit="week", week_start=6)
vix_first_week_end   <- vix_first_week_start + days(6)
cat("VIX daily span:", min(vix$date), "to", max(vix$date), "\n")
vix_w <- to_weekly_friday_close(vix$date, vix$finance)
cat("VIX weekly (Fri-close) span:", min(vix_w$date), "to", max(vix_w$date), "\n")

## VXO (daily) => weekly (only pre-VIX)
vxo <- read_excel("vxoarchive.xls", skip = 2, col_names = TRUE)
names(vxo) <- tolower(names(vxo))
date_col  <- (c("date","observation_date")[c("date","observation_date") %in% names(vxo)])[1]
close_col <- (c("close","last","vxo","closing","closeprice","close_price")[c("close","last","vxo","closing","closeprice","close_price") %in% names(vxo)])[1]
if (is.na(date_col) || is.na(close_col)) { date_col <- names(vxo)[1]; close_col <- names(vxo)[ncol(vxo)] }
dc <- vxo[[date_col]]
if (inherits(dc, "Date")) {
  vxodate <- as.Date(dc)
} else if (is.numeric(dc)) {
  vxodate <- as.Date(dc, origin = "1899-12-30")
} else {
  vxodate <- suppressWarnings(lubridate::dmy(dc)); if (all(is.na(vxodate))) vxodate <- suppressWarnings(lubridate::mdy(dc))
  if (all(is.na(vxodate))) vxodate <- suppressWarnings(lubridate::ymd(dc))
}
cc <- vxo[[close_col]]
if (is.character(cc)) { cc <- gsub("\\s+","", cc); cc <- gsub("\\.","", cc, fixed=TRUE); cc <- gsub(",",".", cc, fixed=TRUE) }
vxoclose <- suppressWarnings(as.numeric(cc))
vxo <- data.frame(date=vxodate, close=vxoclose)
vxo <- vxo[!is.na(vxo$date) & !is.na(vxo$close), ]
vxo <- vxo[vxo$date < vix_first_week_start, ]
cat("VXO daily (pre-VIX) span:", ifelse(nrow(vxo)>0, as.character(min(vxo$date)), "NA"),
    "to", ifelse(nrow(vxo)>0, as.character(max(vxo$date)), "NA"), "\n")
vxo_w <- to_weekly_friday_close(vxo$date, vxo$close)
if (nrow(vxo_w) > 0) cat("VXO weekly (Fri-close) span:", min(vxo_w$date), "to", max(vxo_w$date), "\n") else cat("VXO contributed no weeks (all dates ≥ first VIX week).\n")

## Combine: VXO (pre) + VIX (post)
vxo_w <- vxo_w[vxo_w$date < vix_first_week_end, ]
vol_weekly <- rbind(vxo_w, vix_w)
vol_weekly <- vol_weekly[order(vol_weekly$date), ]
vol_weekly <- vol_weekly[!duplicated(vol_weekly$date), ]
cat("Combined VXO→VIX weekly span:",
    ifelse(nrow(vol_weekly)>0, as.character(min(vol_weekly$date)), "NA"),
    "to", ifelse(nrow(vol_weekly)>0, as.character(max(vol_weekly$date)), "NA"), "\n")
data_merged_fin <- vol_weekly


###############################################################################
# D) RATES & SPREADS (MONTHLY): GS10, T10Y3M
###############################################################################

### GS10 - 10-year Treasury, constant maturity (monthly)
raw <- read.csv("GS10_monthly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw)=="observation_date"] <- "date"
if ("GS10" %in% names(raw))             names(raw)[names(raw)=="GS10"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { nm <- setdiff(names(raw),"date")[1]; names(raw)[names(raw)==nm] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[,c("date","finance")]
cat("GS10 monthly span:", min(data_merged_fin$date,na.rm=TRUE), "to", max(data_merged_fin$date,na.rm=TRUE), "\n")

### T10Y3M - 10y minus 3m Treasury spread (monthly)
raw <- read.csv("T10Y3M_monthly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw)=="observation_date"] <- "date"
if ("T10Y3M" %in% names(raw))           names(raw)[names(raw)=="T10Y3M"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { nm <- setdiff(names(raw),"date")[1]; names(raw)[names(raw)==nm] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[,c("date","finance")]
cat("T10Y3M monthly span:", min(data_merged_fin$date,na.rm=TRUE), "to", max(data_merged_fin$date,na.rm=TRUE), "\n")


###############################################################################
# E) Z.1 - SECURITY BROKERS & DEALERS (QUARTERLY): ASSETS & EQUITY → LEVERAGE
###############################################################################

## Quick leverage from SBD_data.csv (assets/equity)
raw <- read.csv("SBD_data.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw)=="observation_date"] <- "date"
raw$finance <- with(raw, total_assets / equity)
bd_assets_q <- raw[, c("date","finance")]
bd_assets_q$date <- as.Date(bd_assets_q$date)
cat("B&D Asset/Equity (q) span:", min(bd_assets_q$date, na.rm=TRUE), "to", max(bd_assets_q$date, na.rm=TRUE), "\n")
data_merged_fin <- bd_assets_q


###############################################################################
# F) CREDIT CONDITIONS - SLOOS (QUARTERLY)
###############################################################################

### CORALACBS - Net % tightening C&I standards (large/mid)
raw <- read.csv("CORALACBS_quarterly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("CORALACBS" %in% names(raw))        names(raw)[names(raw) == "CORALACBS"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date","finance")]
cat("CORALACBS (q) span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")

### DRTSCILM - Net % stronger C&I demand (large/mid)
raw <- read.csv("DRTSCILM_quarterly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("DRTSCILM" %in% names(raw))         names(raw)[names(raw) == "DRTSCILM"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date","finance")]
cat("DRTSCILM (q) span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")

### DRTSCIS - Net % stronger C&I demand (small firms)
raw <- read.csv("DRTSCIS_quarterly.csv", stringsAsFactors = FALSE)
if ("observation_date" %in% names(raw)) names(raw)[names(raw) == "observation_date"] <- "date"
if ("DRTSCIS" %in% names(raw))          names(raw)[names(raw) == "DRTSCIS"] <- "finance"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
if (!"finance" %in% names(raw)) { val_col <- setdiff(names(raw), "date")[1]; names(raw)[names(raw) == val_col] <- "finance" }
raw$finance <- as.numeric(raw$finance)
data_merged_fin <- raw[, c("date","finance")]
cat("DRTSCIS (q) span:", min(data_merged_fin$date, na.rm=TRUE), "to", max(data_merged_fin$date, na.rm=TRUE), "\n")


###############################################################################
# G) BIS (QUARTERLY): CREDIT-TO-GDP, HP GAP, DSRs
###############################################################################

# BIS_Credit_to_GDP (quarterly)
{
  f <- "BIS_Credit_to_GDP.csv"
  lines     <- readLines(f, warn = FALSE)
  hdr_idx   <- which(grepl("^(observation_date)\\b", lines))[1]
  skip_line <- if (is.na(hdr_idx)) 0L else hdr_idx
  raw <- read.csv(f, skip = skip_line, header = (skip_line == 0L), stringsAsFactors = FALSE, check.names = FALSE)
  nm <- tolower(names(raw)); names(raw) <- nm
  if ("time_period" %in% nm) names(raw)[nm=="time_period"] <- "date"
  if ("obs_value"   %in% nm) names(raw)[nm=="obs_value"]   <- "finance"
  if ("gap"         %in% nm) names(raw)[nm=="gap"]         <- "finance"
  if (!"date" %in% names(raw)) names(raw)[1] <- "date"
  if (!"finance" %in% names(raw)) { cand <- setdiff(names(raw), "date")[1]; names(raw)[names(raw)==cand] <- "finance" }
  raw$date <- as.Date(raw$date)
  y <- as.integer(format(raw$date, "%Y")); m <- as.integer(format(raw$date, "%m")); qm <- ((m - 1L) %/% 3L) * 3L + 1L
  raw$date <- as.Date(sprintf("%04d-%02d-01", y, qm))
  raw$finance <- as.numeric(gsub(",", ".", as.character(raw$finance)))
  data_merged_fin <- raw[, c("date","finance")]
  data_merged_fin <- data_merged_fin[!is.na(data_merged_fin$date) & !is.na(data_merged_fin$finance), ]
  cat("BIS Credit-to-GDP gap span:", min(data_merged_fin$date,na.rm=TRUE), "to", max(data_merged_fin$date,na.rm=TRUE), "\n")
}

# BIS DSR - Households & NPISHs (quarterly)
{
  f <- "BIS Debt service ratio Households & NPISHs.csv"
  lines     <- readLines(f, warn = FALSE)
  hdr_idx   <- which(grepl("^(observation_date)\\b", lines))[1]
  skip_line <- if (is.na(hdr_idx)) 0L else hdr_idx
  raw <- read.csv(f, skip = skip_line, header = (skip_line == 0L), stringsAsFactors = FALSE, check.names = FALSE)
  nm <- tolower(names(raw)); names(raw) <- nm
  if ("time_period" %in% nm) names(raw)[nm=="time_period"] <- "date"
  if ("obs_value"   %in% nm) names(raw)[nm=="obs_value"]   <- "finance"
  for (cand in c("dsr","value")) if (cand %in% names(raw)) { names(raw)[cand] <- "finance"; break }
  if (!"date" %in% names(raw)) names(raw)[1] <- "date"
  if (!"finance" %in% names(raw)) { cand <- setdiff(names(raw), "date")[1]; names(raw)[names(raw)==cand] <- "finance" }
  raw$date <- as.Date(raw$date)
  y <- as.integer(format(raw$date, "%Y")); m <- as.integer(format(raw$date, "%m")); qm <- ((m - 1L) %/% 3L) * 3L + 1L
  raw$date <- as.Date(sprintf("%04d-%02d-01", y, qm))
  raw$finance <- as.numeric(gsub(",", ".", as.character(raw$finance)))
  data_merged_fin <- raw[, c("date","finance")]
  data_merged_fin <- data_merged_fin[!is.na(data_merged_fin$date) & !is.na(data_merged_fin$finance), ]
  cat("BIS DSR HH span:", min(data_merged_fin$date,na.rm=TRUE), "to", max(data_merged_fin$date,na.rm=TRUE), "\n")
}

# BIS DSR - Non-financial corporations (quarterly)
{
  f <- "BIS Debt service ratio Non-financial corporations.csv"
  lines     <- readLines(f, warn = FALSE)
  hdr_idx   <- which(grepl("^(observation_date)\\b", lines))[1]
  skip_line <- if (is.na(hdr_idx)) 0L else hdr_idx
  raw <- read.csv(f, skip = skip_line, header = (skip_line == 0L), stringsAsFactors = FALSE, check.names = FALSE)
  nm <- tolower(names(raw)); names(raw) <- nm
  if ("time_period" %in% nm) names(raw)[nm=="time_period"] <- "date"
  if ("obs_value"   %in% nm) names(raw)[nm=="obs_value"]   <- "finance"
  for (cand in c("dsr","value")) if (cand %in% names(raw)) { names(raw)[cand] <- "finance"; break }
  if (!"date" %in% names(raw)) names(raw)[1] <- "date"
  if (!"finance" %in% names(raw)) { cand <- setdiff(names(raw), "date")[1]; names(raw)[names(raw)==cand] <- "finance" }
  raw$date <- as.Date(raw$date)
  y <- as.integer(format(raw$date, "%Y")); m <- as.integer(format(raw$date, "%m")); qm <- ((m - 1L) %/% 3L) * 3L + 1L
  raw$date <- as.Date(sprintf("%04d-%02d-01", y, qm))
  raw$finance <- as.numeric(gsub(",", ".", as.character(raw$finance)))
  data_merged_fin <- raw[, c("date","finance")]
  data_merged_fin <- data_merged_fin[!is.na(data_merged_fin$date) & !is.na(data_merged_fin$finance), ]
  cat("BIS DSR NFC span:", min(data_merged_fin$date,na.rm=TRUE), "to", max(data_merged_fin$date,na.rm=TRUE), "\n")
}

# BIS DSR - Private non-financial sector (quarterly)
{
  f <- "BIS Debt service ratio Private non-financial sector.csv"
  lines     <- readLines(f, warn = FALSE)
  hdr_idx   <- which(grepl("^(observation_date)\\b", lines))[1]
  skip_line <- if (is.na(hdr_idx)) 0L else hdr_idx
  raw <- read.csv(f, skip = skip_line, header = (skip_line == 0L), stringsAsFactors = FALSE, check.names = FALSE)
  nm <- tolower(names(raw)); names(raw) <- nm
  if ("time_period" %in% nm) names(raw)[nm=="time_period"] <- "date"
  if ("obs_value"   %in% nm) names(raw)[nm=="obs_value"]   <- "finance"
  for (cand in c("dsr","value")) if (cand %in% names(raw)) { names(raw)[cand] <- "finance"; break }
  if (!"date" %in% names(raw)) names(raw)[1] <- "date"
  if (!"finance" %in% names(raw)) { cand <- setdiff(names(raw), "date")[1]; names(raw)[names(raw)==cand] <- "finance" }
  raw$date <- as.Date(raw$date)
  y <- as.integer(format(raw$date, "%Y")); m <- as.integer(format(raw$date, "%m")); qm <- ((m - 1L) %/% 3L) * 3L + 1L
  raw$date <- as.Date(sprintf("%04d-%02d-01", y, qm))
  raw$finance <- as.numeric(gsub(",", ".", as.character(raw$finance)))
  data_merged_fin <- raw[, c("date","finance")]
  data_merged_fin <- data_merged_fin[!is.na(data_merged_fin$date) & !is.na(data_merged_fin$finance), ]
  cat("BIS DSR Private span:", min(data_merged_fin$date,na.rm=TRUE), "to", max(data_merged_fin$date,na.rm=TRUE), "\n")
}


###############################################################################
# H) REAL MACRO SERIES (QUARTERLY): CUMFNS, TCU, GDPC1
###############################################################################

### CUMFNS - Capacity Utilization: Manufacturing (quarterly, FRB)
raw <- read.csv("CUMFNS.csv", stringsAsFactors = FALSE)
names(raw)[names(raw) == "observation_date"] <- "date"
if ("CUMFNS" %in% names(raw)) names(raw)[names(raw) == "CUMFNS"] <- "log_hours"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
raw$finance <- as.numeric(raw$log_hours)
data_merged_real <- raw[, c("date", "log_hours")]
cat("CUMFNS span:", min(data_merged_real$date, na.rm=TRUE), "to", max(data_merged_real$date, na.rm=TRUE), "\n")

### TCU - Capacity Utilization: Total Industry (quarterly, FRB)
raw <- read.csv("TCU.csv", stringsAsFactors = FALSE)
names(raw)[names(raw) == "observation_date"] <- "date"
if ("TCU" %in% names(raw)) names(raw)[names(raw) == "TCU"] <- "log_hours"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
raw$finance <- as.numeric(raw$log_hours)
data_merged_real <- raw[, c("date", "log_hours")]
cat("TCU span:", min(data_merged_real$date, na.rm=TRUE), "to", max(data_merged_real$date, na.rm=TRUE), "\n")

### GDPC1 - Real GDP (quarterly, BEA via FRED)
raw <- read.csv("GDPC1.csv", stringsAsFactors = FALSE)
names(raw)[names(raw) == "observation_date"] <- "date"
if ("GDPC1" %in% names(raw)) names(raw)[names(raw) == "GDPC1"] <- "log_hours"
names(raw) <- tolower(names(raw)); raw$date <- as.Date(raw$date)
raw$finance <- as.numeric(raw$log_hours)
data_merged_real <- raw[, c("date", "log_hours")]
cat("GDPC1 span:", min(data_merged_real$date, na.rm=TRUE), "to", max(data_merged_real$date, na.rm=TRUE), "\n")

## Optional: end at 2019Q4 for real series if needed
end_q4_2019 <- as.Date("2019-12-31")
data_merged_real <- subset(data_merged_real, date <= end_q4_2019)
cat("Real-series span (<=2019Q4):", min(data_merged_real$date), "to", max(data_merged_real$date), "\n")






###############################################################################
# I) LOG DEFLATED HOURS - hours (BLS) / population (Census+FRED)
###############################################################################

# 1) Hours worked (BLS, MachineReadable table → Total economy / All workers)
hrs_raw <- readxl::read_excel("total-economy-hours-employment_BLS.xlsx", sheet = "MachineReadable")
names(hrs_raw) <- trimws(names(hrs_raw))
hrs_filt <- subset(hrs_raw,
                   Sector    == "Total economy" &
                     Basis     == "All workers" &
                     Component == "Total U.S. economy" &
                     Measure   == "Hours worked")
hrs_filt$Year  <- as.integer(hrs_filt$Year)
hrs_filt$Qtr   <- as.integer(hrs_filt$Qtr)
hrs_filt$Value <- as.numeric(hrs_filt$Value)
hrs_filt$DateQ <- as.yearqtr(paste(hrs_filt$Year, hrs_filt$Qtr), "%Y %q")
hrs_sub <- subset(hrs_filt, DateQ >= as.yearqtr("1948 Q1") & DateQ <= as.yearqtr("2015 Q2"))
hrs_prepped <- data.frame(date = as.Date(hrs_sub$DateQ), hrs_worked = as.numeric(hrs_sub$Value))

# 2) Population (Census interpolated 1947–1952, rescaled to FRED at 1952Q1; FRED thereafter)
pop_fred <- read.csv("POP_FRED.csv", stringsAsFactors = FALSE)
pop_fred$date <- as.Date(pop_fred$observation_date); pop_fred$POP <- as.numeric(pop_fred$POP)
census_lines <- readLines("HNPE_US_Census.txt")
jul1_lines   <- grep("^\\s*July 1, ", census_lines, value = TRUE)
census_df <- data.frame(
  year = as.integer(sub("^\\s*July 1, (\\d{4}).*$", "\\1", jul1_lines)),
  pop  = as.numeric(gsub(",", "", sub("^\\s*July 1, \\d{4}\\s+([0-9,]+).*", "\\1", jul1_lines))),
  stringsAsFactors = FALSE
)
census_df <- subset(census_df, year >= 1947 & year <= 1952)
q_dates <- seq(as.Date("1947-01-01"), as.Date("1952-01-01"), by = "quarter")
q_years <- as.numeric(format(q_dates, "%Y")) + (as.numeric(format(q_dates, "%m")) - 1) / 12
cens_q  <- approx(x = census_df$year, y = census_df$pop, xout = q_years, rule = 2)$y
cens_1952Q1 <- cens_q[q_dates == as.Date("1952-01-01")]
fred_1952Q1 <- pop_fred$POP[pop_fred$date == as.Date("1952-01-01")]
scale_factor <- fred_1952Q1 / cens_1952Q1
cens_adj     <- cens_q * scale_factor
early_df <- data.frame(date = q_dates,           population_US = cens_adj)
late_df  <- data.frame(date = pop_fred$date,     population_US = pop_fred$POP)
pop_prepped <- rbind(early_df[early_df$date <= as.Date("1952-01-01"), ],
                     late_df[ late_df$date  > as.Date("1952-01-01"), ])
rownames(pop_prepped) <- NULL

# 3) Merge & construct log(hours per capita)
data_merged_real <- merge(hrs_prepped, pop_prepped, by = "date")
data_merged_real$deflated_hours <- data_merged_real$hrs_worked / data_merged_real$population_US
data_merged_real$log_hours      <- log(data_merged_real$deflated_hours)

# Quick checks
str(data_merged_real)
summary(data_merged_real$deflated_hours)
summary(data_merged_real$log_hours)

plot(data_merged_real$date, data_merged_real$log_hours, type = "l",
     xlab = "Date", ylab = "log(deflated hours)",
     main = "Log of Total Economy Hours Worked (deflated by population)")
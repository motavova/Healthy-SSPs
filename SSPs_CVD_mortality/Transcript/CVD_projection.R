# =============================================================================
# CVD_PROJECTIONS
# Updated: 14/08/2026
# Author: Martina Otavova
# =============================================================================

.libPaths(c("/home/otavova/R/library"))  # point R at the project-specific package library

required_packages <- c("dplyr", "tidyr", "zoo", "lme4", "lmerTest", "countrycode", "ggplot2")  # splines ships with base R, no install needed
missing_packages <- setdiff(required_packages, rownames(installed.packages()))       # figure out which of those aren't installed yet
if (length(missing_packages) > 0) install.packages(missing_packages)                   # download and install anything missing

library(dplyr)       # data wrangling (mutate, filter, joins, group_by, etc.)
library(tidyr)        # tidying helpers (not directly used below but loaded for the pipeline)
library(zoo)           # na.approx() for linear interpolation of missing population values
library(lme4)            # lmer() mixed-effects model fitting
library(lmerTest)         # adds p-values/df to lmer models (loaded for downstream inspection)
library(countrycode)       # ISO3 <-> UN region/subregion lookups
library(splines)             # ns() natural cubic spline basis for the year trend
library(ggplot2)               # country-level plots for the group3 chapter (chapter 10)

# -----------------------------------------------------------------------------
# 1) PATHWAYS AND CONSTANTS
# -----------------------------------------------------------------------------

input_dir  <- "/home/otavova/Healthy-SSPs/SSPs_CVD_mortality/Data"                # where all source .rds/.rda/.csv files live
output_dir <- "/home/otavova/Healthy-SSPs/SSPs_CVD_mortality/Results/Projections" # where the final projection gets saved

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)  # create output_dir if it doesn't exist yet (silently)

eps <- 1e-12                                    # small constant used to keep ratios/rates away from exactly 0
covid_years <- c(2020, 2021, 2022)               # years excluded from model fitting (COVID-distorted mortality)
spline_df <- 3                                    # degrees of freedom for the natural spline on year
spline_term <- paste0("ns(year_c, df = ", spline_df, ")")  # build the spline formula term as a string, e.g. "ns(year_c, df = 3)"

risk_factor_cols <- c("bmi", "smoking", "fpg", "sbp", "ldl", "sodium", "albuminuria", "ckd5")  # candidate risk-factor covariates, availability varies by age group

base_required_cols <- c(              # columns that must be non-NA for a row to be usable when FITTING the model
  "log_ratio_cardio", "death", "allcause_death",
  "gdp", "edu", "urban", "tfr",
  "sex", "iso3", "region", "pop"
)

proj_base_required_cols <- c(         # columns that must be non-NA for a row to be usable when PROJECTING
  "pop", "gdp", "edu", "urban", "tfr",
  "mx_allcause", "region"
)

standardize_value <- function(x, mean_x, sd_x) {   # z-score a vector using externally supplied mean/sd (so train and projection use the same scale)
  if (is.na(sd_x) || sd_x == 0) return(rep(NA_real_, length(x)))  # guard against degenerate (missing or zero-variance) sd
  (x - mean_x) / sd_x                                              # standard z-score formula
}

make_un_subregion <- function(iso3) {                                                          # map an ISO3 country code to its UN subregion name
  region <- countrycode(iso3, origin = "iso3c", destination = "un.regionsub.name", warn = FALSE)  # lookup via countrycode package
  if_else(iso3 == "TWN", "Eastern Asia", region)                                                    # countrycode has no UN region for Taiwan, so hardcode it
}

# -----------------------------------------------------------------------------
# 2) UPLOAD DATA
# -----------------------------------------------------------------------------
cardio <- readRDS(file.path(input_dir, "cardio.rds"))  # observed CVD & all-cause deaths by country/age/sex/year (WPP-based)

pop <- readRDS(file.path(input_dir, "pop.rds")) %>%   # population by scenario/country/age/sex/year (already interpolated across years upstream)
  group_by(scenario, iso3, age, sex) %>%                      # interpolate any remaining gaps within each country/age/sex/scenario series
  arrange(year, .by_group = TRUE) %>%                          # make sure years are in order before interpolating
  mutate(
    pop = ifelse(pop == 0, NA_real_, pop),    # treat a literal 0 population as missing, not a real value
    pop = na.approx(pop, na.rm = FALSE),        # linearly interpolate across the missing (NA) years
    pop = ifelse(is.na(pop), 0.001, pop)          # anything still NA (e.g. leading/trailing gaps) gets a tiny floor instead of NA
  ) %>%
  ungroup()                                          # drop the grouping, back to a plain data frame

covariates <- readRDS(file.path(input_dir, "edu_gdp_urban.rds"))  # education/urbanization/GDP by scenario/country/age/sex/year


tfr <- readRDS(file.path(input_dir, "tfr.rds")) %>% select(-name)  # total fertility rate by scenario/country/year; drop the country-name column, iso3 is enough

mx_proj <- readRDS(file.path(input_dir, "ltable_mx.rds")) %>%  # projected all-cause mortality rates (Samir source), only available every 5 years
  mutate(mx_allcause = mx)                                                 # rename to the column name used everywhere else in this script


risk_factors <- readRDS(file.path(input_dir, "ihme_risk_factors.rds")) %>%  # smoothed risk-factor trajectories 1990-2100, all 5 SSPs
  select(iso3, age, sex, year, scenario, all_of(risk_factor_cols))            # keep only the identifying columns + the 8 candidate risk factors


rf_hist <- risk_factors %>%          # the historical (observed-period) slice of risk factors, used when fitting the model
  filter(scenario == "SSP1") %>%      # all SSPs share identical values before 2024, so SSP1 is an arbitrary-but-valid pick for the past
  select(-scenario)                     # scenario is no longer needed once we've filtered to one

# -----------------------------------------------------------------------------
# 3) MODIFY DATA FOR FITTING THE MODEL
# -----------------------------------------------------------------------------

obs_base <- cardio %>%                                                     # start from observed CVD/all-cause deaths
  mutate(scenario = "SSP1") %>%                                              # tag historical rows as SSP1 so they can be joined against scenario-keyed tables
  left_join(covariates, by = c("scenario", "iso3", "age", "sex", "year")) %>%  # attach education/urbanization/GDP
  left_join(tfr, by = c("scenario", "iso3", "year")) %>%                       # attach total fertility rate
  left_join(rf_hist, by = c("iso3", "age", "sex", "year")) %>%                  # attach historical risk factors (bmi, smoking, etc.)
  mutate(
    pop = pmax(pop, 0.001),                                    # floor population at a tiny positive value to avoid division/log issues
    gdp = pmax(gdp, 1),                                          # floor GDP at 1 so log(gdp) below is always defined
    year_c = year - 1990,                                          # center year at 1990 for numerical stability in the model
    region = make_un_subregion(iso3),                                # derive UN subregion from ISO3 for the random-effects grouping
    ratio_cardio = death / allcause_death,                             # share of all-cause deaths that are cardiovascular
    ratio_cardio = pmin(pmax(ratio_cardio, eps), 1 - eps),               # clamp the ratio strictly inside (0, 1) so log() is defined
    log_ratio_cardio = log(ratio_cardio)                                   # the actual outcome variable the model predicts
  ) %>%
  filter(!year %in% covid_years)     # drop COVID years, their mortality patterns would distort the fitted trend

obs_complete <- obs_base %>%                                          # restrict to rows with everything the model needs
  filter(complete.cases(across(all_of(base_required_cols)))) %>%       # drop any row missing a required covariate/outcome
  mutate(
    sex = factor(sex),          # convert to factor for use as a categorical fixed effect
    iso3 = factor(iso3),          # convert to factor; its levels define which countries the model "knows"
    region = factor(region),        # convert to factor; used for the nested random effect
    log_gdp = log(gdp)                # log-transform GDP (already floored at 1 above, so this is always finite)
  )

std_full <- obs_complete %>%                                                                            # compute mean/sd of every covariate across the FULL fitting sample
  summarise(                                                                                              # these standardization parameters get reused for projection too, so train/projection share one scale
    edu_mean = mean(edu, na.rm = TRUE), edu_sd = sd(edu, na.rm = TRUE),
    urban_mean = mean(urban, na.rm = TRUE), urban_sd = sd(urban, na.rm = TRUE),
    log_gdp_mean = mean(log_gdp, na.rm = TRUE), log_gdp_sd = sd(log_gdp, na.rm = TRUE),
    tfr_mean = mean(tfr, na.rm = TRUE), tfr_sd = sd(tfr, na.rm = TRUE),
    bmi_mean = mean(bmi, na.rm = TRUE), bmi_sd = sd(bmi, na.rm = TRUE),
    smoking_mean = mean(smoking, na.rm = TRUE), smoking_sd = sd(smoking, na.rm = TRUE),
    fpg_mean = mean(fpg, na.rm = TRUE), fpg_sd = sd(fpg, na.rm = TRUE),
    sbp_mean = mean(sbp, na.rm = TRUE), sbp_sd = sd(sbp, na.rm = TRUE),
    ldl_mean = mean(ldl, na.rm = TRUE), ldl_sd = sd(ldl, na.rm = TRUE),
    sodium_mean = mean(sodium, na.rm = TRUE), sodium_sd = sd(sodium, na.rm = TRUE),
    albuminuria_mean = mean(albuminuria, na.rm = TRUE), albuminuria_sd = sd(albuminuria, na.rm = TRUE),
    ckd5_mean = mean(ckd5, na.rm = TRUE), ckd5_sd = sd(ckd5, na.rm = TRUE)
  )

obs_model_full <- obs_complete %>%                                                                  # add the standardized (z-scored) version of every covariate
  mutate(
    edu_z = standardize_value(edu, std_full$edu_mean, std_full$edu_sd),
    urban_z = standardize_value(urban, std_full$urban_mean, std_full$urban_sd),
    gdp_z = standardize_value(log_gdp, std_full$log_gdp_mean, std_full$log_gdp_sd),
    tfr_z = standardize_value(tfr, std_full$tfr_mean, std_full$tfr_sd),
    bmi_z = standardize_value(bmi, std_full$bmi_mean, std_full$bmi_sd),
    smoking_z = standardize_value(smoking, std_full$smoking_mean, std_full$smoking_sd),
    fpg_z = standardize_value(fpg, std_full$fpg_mean, std_full$fpg_sd),
    sbp_z = standardize_value(sbp, std_full$sbp_mean, std_full$sbp_sd),
    ldl_z = standardize_value(ldl, std_full$ldl_mean, std_full$ldl_sd),
    sodium_z = standardize_value(sodium, std_full$sodium_mean, std_full$sodium_sd),
    albuminuria_z = standardize_value(albuminuria, std_full$albuminuria_mean, std_full$albuminuria_sd),
    ckd5_z = standardize_value(ckd5, std_full$ckd5_mean, std_full$ckd5_sd)
  )

# -----------------------------------------------------------------------------
# 4) FIT THE MODEL
# -----------------------------------------------------------------------------
models_by_age <- obs_model_full %>%     # fit one mixed model per age group, since available risk factors differ by age
  group_by(age) %>%                       # split the data by age group
  group_modify(~ {                          # for each age group's slice (.x) and key (.y), return a one-row summary tibble
    age_label <- unique(.y$age)                # the age-group label for this iteration (e.g. "0-4")
    age_rfs <- risk_factor_cols[vapply(           # figure out which of the 8 candidate risk factors actually have data for this age
      .x %>% select(all_of(risk_factor_cols)),
      function(column) any(!is.na(column)),
      logical(1)
    )]

    required_cols <- c(base_required_cols, age_rfs)                    # this age's full set of required (non-NA) columns
    obs_age <- .x %>% filter(complete.cases(across(all_of(required_cols))))  # rows usable for fitting this age's model

    if (nrow(obs_age) == 0) {                    # if nothing survives the completeness filter, skip fitting entirely
      return(tibble(
        model = list(NULL),      # no model was fit
        age_rfs = list(character()),  # no risk factors were used
        resid = list(tibble())          # no residuals to carry forward
      ))
    }

    age_rfs_z <- paste0(age_rfs, "_z")            # z-scored column names for this age's available risk factors
    formula_terms <- c(                             # assemble the right-hand side of the model formula
      "edu_z", "urban_z", "gdp_z", "tfr_z",
      age_rfs_z,
      "urban_z:gdp_z", "sex", spline_term
    )
    model_formula <- as.formula(                     # build the full lmer formula, fixed + random effects
      paste(
        "log_ratio_cardio ~",
        paste(formula_terms, collapse = " + "),
        "+ (1 | region / iso3) + (0 + year_c | iso3)"  # country nested in region intercept, plus per-country year slope
      )
    )

    fitted_model <- lmer(                    # fit the mixed-effects model for this age group
      model_formula,
      data = obs_age,
      REML = FALSE,                            # use ML (not REML) so nested models remain comparable
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))  # more robust optimizer settings for convergence
    )

    residuals_by_group <- obs_age %>%                       # compute each country/sex's most recent residual, to carry forward into projections
      mutate(
        fitted_lmer = predict(                                # in-sample fitted values including random effects up to country level
          fitted_model,
          newdata = obs_age,
          re.form = ~ (1 | region / iso3),
          allow.new.levels = TRUE
        ),
        resid_lmer = log_ratio_cardio - fitted_lmer             # residual = observed minus fitted, on the log-ratio scale
      ) %>%
      group_by(iso3, sex) %>%                                     # one residual trajectory per country x sex
      arrange(year, .by_group = TRUE) %>%                           # order by year so "last" below is the most recent
      summarise(
        last_resid = last(resid_lmer),        # the most recent residual, used to anchor projections at the last observed value
        last_train_year = max(year),            # the year that residual corresponds to
        .groups = "drop"
      ) %>%
      mutate(age = age_label)      # tag with the age group (this is a nested tibble, not the top-level group_modify output)

    tibble(                              # the actual row returned for this age group by group_modify
      model = list(fitted_model),           # the fitted lmer model object
      age_rfs = list(age_rfs),                # which risk factors were used for this age
      resid = list(residuals_by_group)          # per-country/sex last residual + last observed year
    )
  }) %>%
  ungroup()   # drop the grouping, models_by_age is now one row per age group

# -----------------------------------------------------------------------------
# 5) RUN PROJECTIONS
# -----------------------------------------------------------------------------

proj_iso3 <- unique(obs_complete$iso3)   # only project for countries that had usable data when fitting

proj_frame_base <- pop %>%                                    # start from interpolated population (annual, not just every 5 years)
  filter(year > 2023, iso3 %in% proj_iso3) %>%                   # future years only, restricted to countries in the fitting sample
  mutate(
    pop = pop * 1000,             # pop.rds stores population in thousands; convert to actual head-count units
    pop = pmax(pop, 0.001)          # floor again after rescaling, just in case
  ) %>%
  left_join(covariates, by = c("scenario", "iso3", "age", "sex", "year")) %>%   # attach projected education/urbanization/GDP
  left_join(tfr, by = c("scenario", "iso3", "year")) %>%                        # attach projected fertility rate
  left_join(mx_proj, by = c("scenario", "iso3", "age", "sex", "year")) %>%       # attach projected all-cause mortality rate (only present every 5 years)
  left_join(risk_factors, by = c("scenario", "iso3", "age", "sex", "year")) %>%   # attach projected risk factors, all 5 SSPs
  mutate(
    gdp = pmax(gdp, 1),                                          # floor GDP again so log_gdp below is always defined
    mx_allcause = ifelse(mx_allcause == 0, NA_real_, mx_allcause),  # treat a literal 0 mortality rate as missing, not real
    mx_allcause = pmax(mx_allcause, eps),                            # floor to avoid a zero death-rate downstream (NA stays NA)
    year_c = year - 1990,                                              # same centering used when fitting
    region = make_un_subregion(iso3),                                    # same region derivation used when fitting
    log_gdp = log(gdp),                                                    # same log-GDP transform used when fitting
    edu_z = standardize_value(edu, std_full$edu_mean, std_full$edu_sd),          # standardize using the SAME mean/sd as the fitting sample
    urban_z = standardize_value(urban, std_full$urban_mean, std_full$urban_sd),
    gdp_z = standardize_value(log_gdp, std_full$log_gdp_mean, std_full$log_gdp_sd),
    tfr_z = standardize_value(tfr, std_full$tfr_mean, std_full$tfr_sd),
    bmi_z = standardize_value(bmi, std_full$bmi_mean, std_full$bmi_sd),
    smoking_z = standardize_value(smoking, std_full$smoking_mean, std_full$smoking_sd),
    fpg_z = standardize_value(fpg, std_full$fpg_mean, std_full$fpg_sd),
    sbp_z = standardize_value(sbp, std_full$sbp_mean, std_full$sbp_sd),
    ldl_z = standardize_value(ldl, std_full$ldl_mean, std_full$ldl_sd),
    sodium_z = standardize_value(sodium, std_full$sodium_mean, std_full$sodium_sd),
    albuminuria_z = standardize_value(albuminuria, std_full$albuminuria_mean, std_full$albuminuria_sd),
    ckd5_z = standardize_value(ckd5, std_full$ckd5_mean, std_full$ckd5_sd)
  )

proj_spline <- proj_frame_base %>%      # generate the actual projected cardiovascular death counts
  group_by(age) %>%                       # process one age group at a time, matching that age's fitted model
  group_modify(~ {                          # .x = this age's projection rows, .y = the age-group key
    age_label <- unique(.y$age)                # this iteration's age label
    model_row <- models_by_age %>% filter(age == age_label)   # look up the model fitted for this age group

    if (nrow(model_row) == 0L || is.null(model_row$model[[1]])) {  # if no model was fit for this age (e.g. no data), skip it
      return(tibble())
    }

    age_rfs <- model_row$age_rfs[[1]]                                          # risk factors this age's model actually uses
    required_cols <- setdiff(c(proj_base_required_cols, age_rfs), "mx_allcause")  # completeness check excludes mx_allcause on purpose, so years without 5-year mortality data aren't dropped entirely

    proj_age <- .x %>%
      filter(complete.cases(across(all_of(required_cols)))) %>%   # keep only rows with all OTHER required covariates present
      mutate(
        sex = factor(sex, levels = levels(obs_complete$sex)),        # align factor levels with the ones the model was trained on
        iso3 = factor(iso3, levels = levels(obs_complete$iso3)),        # unseen countries become NA rather than new factor levels
        region = factor(region, levels = levels(obs_complete$region))
      ) %>%
      filter(!is.na(sex), !is.na(iso3), !is.na(region))    # drop anything that fell outside the trained factor levels

    if (nrow(proj_age) == 0L) {    # nothing left to project for this age group
      return(tibble())
    }

    fitted_model <- model_row$model[[1]]     # the lmer model fitted for this age group
    proj_age <- proj_age %>%
      mutate(
        pred_log_m_cardio_base = predict(       # model's baseline prediction (fixed effects + region/country intercepts)
          fitted_model,
          newdata = proj_age,
          re.form = ~ (1 | region / iso3),          # deliberately excludes the per-country year-slope random effect here
          allow.new.levels = TRUE                     # allow country/region levels the model hasn't seen without erroring
        )
      )

    residuals_by_group <- model_row$resid[[1]]     # per-country/sex last observed residual + last observed year, from step 4

    proj_age %>%
      left_join(residuals_by_group, by = c("iso3", "sex")) %>%    # attach each country/sex's carried-forward residual
      mutate(
        last_resid = coalesce(last_resid, 0),                 # countries with no fitted residual (new/unseen) get zero offset
        last_train_year = coalesce(last_train_year, 2023L),      # and default to anchoring at the last historical year
        years_ahead = year - last_train_year,                      # how far this projection year is from the anchor year
        resid_damped = last_resid * exp(-0.05 * years_ahead),         # exponentially decay the carried-forward residual over time
        pred_log_ratio = pred_log_m_cardio_base + resid_damped,         # final predicted log cardio-share, baseline + decayed residual
        pred_ratio_cardio = pmin(pmax(exp(pred_log_ratio), eps), 0.85),   # back-transform and clamp the share into a plausible (eps, 0.85) range
        pred_allcause_deaths = mx_allcause * pop,                           # NA here whenever mx_allcause is NA (i.e. off the 5-year mortality grid)
        pred_cardio_deaths = pred_ratio_cardio * pred_allcause_deaths         # inherits that same NA when all-cause deaths are unknown
      ) %>%
      select(scenario, iso3, sex, year, pop, pred_ratio_cardio, pred_cardio_deaths, pred_allcause_deaths)  # keep only the columns needed downstream
  }) %>%
  ungroup()   # drop the grouping, proj_spline is now one row per scenario/country/sex/age/year


proj_spline <- proj_spline[complete.cases(proj_spline),]   # drop rows still missing anything (e.g. years off the mx_allcause 5-year grid)

head(proj_spline)   # quick sanity check of the final projection frame

# -----------------------------------------------------------------------------
# 6) SAVE OUTCOME PROJECTION IN RESULTS
# -----------------------------------------------------------------------------
saveRDS(proj_spline, file.path(output_dir, "proj_spline.rds"))   # persist the final projection to disk

# -----------------------------------------------------------------------------
# 7) SAVE COMBINED OBSERVED + PROJECTED SERIES
# -----------------------------------------------------------------------------

observed_series <- cardio %>%                                       # raw historical data (1990-2023), independent of the model's completeness filters
  mutate(
    scenario = "historical",             # historical years are pre-SSP-divergence, so tag them as their own pseudo-scenario
    iso3 = as.character(iso3),           # match the type used in the projected block below
    sex = as.character(sex),               # match the type used in the projected block below
    cardio_deaths = death,                   # rename to the shared output column name
    allcause_deaths = allcause_death,          # rename to the shared output column name
    cardio_share = share                         # cardio.rds already carries the observed cardio/all-cause share
  ) %>%
  select(scenario, iso3, age, sex, year, pop, cardio_deaths, allcause_deaths, cardio_share)

projected_series <- proj_spline %>%                                  # model-based projections (2024-2100)
  mutate(
    iso3 = as.character(iso3),         # proj_spline's iso3 is a factor; match observed_series's character type
    sex = as.character(sex),             # proj_spline's sex is a factor; match observed_series's character type
    cardio_deaths = pred_cardio_deaths,      # rename to the shared output column name
    allcause_deaths = pred_allcause_deaths,    # rename to the shared output column name
    cardio_share = pred_ratio_cardio             # rename to the shared output column name
  ) %>%
  select(scenario, iso3, age, sex, year, pop, cardio_deaths, allcause_deaths, cardio_share)

obs_and_proj <- bind_rows(observed_series, projected_series)   # stack historical rows on top of projected rows into one long table

saveRDS(obs_and_proj, file.path(output_dir, "obs_and_proj.rds"))   # persist the combined observed + projected series to disk

# -----------------------------------------------------------------------------
# 8) GROUP3 COUNTRIES (13): SEPARATE MODEL, NO TFR
# -----------------------------------------------------------------------------
# ATG, BRB, BRN, DJI, ERI, GRD, GUM, LBY, MRT, PNG, SYC, UZB, VIR never make it
# into the model above: tfr.rds only has 2025-2100 fertility rates for them (no
# 1990-2024 history), so every historical row fails base_required_cols's tfr
# check; and ihme_risk_factors.rds (the smoothed risk-factor file) has no rows
# for them at all. This chapter fits a second, separate pooled model for just
# these 13 countries -- same model family (natural-spline year term, region/
# country + per-country year-slope random effects) but with tfr dropped from
# the formula/required columns, and risk factors taken from their own raw
# (unsmoothed) file. All-cause mortality still comes from the same ltable_mx.rds
# projection (mx_proj) used above -- that file already covers all 13.

group3_iso3 <- c("ATG", "BRB", "BRN", "DJI", "ERI", "GRD", "GUM",
                 "LBY", "MRT", "PNG", "SYC", "UZB", "VIR")  # the 13 countries this chapter covers

base_required_cols_g3      <- setdiff(base_required_cols, "tfr")       # same completeness check as the main model, minus tfr
proj_base_required_cols_g3 <- setdiff(proj_base_required_cols, "tfr")   # same, for the projection frame

risk_factors_g3 <- readRDS(file.path(input_dir, "ihme_risk_factors_13countries.rds")) %>%  # raw/unsmoothed risk factors, group3 only
  select(iso3, age, sex, year, scenario, all_of(risk_factor_cols))

rf_hist_g3 <- risk_factors_g3 %>%    # historical (SSP1) slice, same convention as rf_hist above
  filter(scenario == "SSP1") %>%
  select(-scenario)

# --- 8a) Build fitting data for the 13 countries (no tfr join, no tfr requirement) ---

obs_base_g3 <- cardio %>%                                                     # observed CVD/all-cause deaths, restricted to group3
  filter(iso3 %in% group3_iso3) %>%
  mutate(scenario = "SSP1") %>%
  left_join(covariates, by = c("scenario", "iso3", "age", "sex", "year")) %>%   # education/urbanization/GDP (already covers group3)
  left_join(rf_hist_g3, by = c("iso3", "age", "sex", "year")) %>%                # group3's own historical risk factors
  mutate(
    pop = pmax(pop, 0.001),
    gdp = pmax(gdp, 1),
    year_c = year - 1990,
    region = make_un_subregion(iso3),
    ratio_cardio = death / allcause_death,
    ratio_cardio = pmin(pmax(ratio_cardio, eps), 1 - eps),
    log_ratio_cardio = log(ratio_cardio)
  ) %>%
  filter(!year %in% covid_years)

obs_complete_g3 <- obs_base_g3 %>%                                            # rows usable for fitting (tfr not required here)
  filter(complete.cases(across(all_of(base_required_cols_g3)))) %>%
  mutate(
    sex = factor(sex),
    iso3 = factor(iso3),
    region = factor(region),
    log_gdp = log(gdp)
  )

std_g3 <- obs_complete_g3 %>%                                                 # standardize using THIS 13-country sample's own mean/sd --
  summarise(                                                                   # a separate model needs its own scale, not the main model's
    edu_mean = mean(edu, na.rm = TRUE), edu_sd = sd(edu, na.rm = TRUE),
    urban_mean = mean(urban, na.rm = TRUE), urban_sd = sd(urban, na.rm = TRUE),
    log_gdp_mean = mean(log_gdp, na.rm = TRUE), log_gdp_sd = sd(log_gdp, na.rm = TRUE),
    bmi_mean = mean(bmi, na.rm = TRUE), bmi_sd = sd(bmi, na.rm = TRUE),
    smoking_mean = mean(smoking, na.rm = TRUE), smoking_sd = sd(smoking, na.rm = TRUE),
    fpg_mean = mean(fpg, na.rm = TRUE), fpg_sd = sd(fpg, na.rm = TRUE),
    sbp_mean = mean(sbp, na.rm = TRUE), sbp_sd = sd(sbp, na.rm = TRUE),
    ldl_mean = mean(ldl, na.rm = TRUE), ldl_sd = sd(ldl, na.rm = TRUE),
    sodium_mean = mean(sodium, na.rm = TRUE), sodium_sd = sd(sodium, na.rm = TRUE),
    albuminuria_mean = mean(albuminuria, na.rm = TRUE), albuminuria_sd = sd(albuminuria, na.rm = TRUE),
    ckd5_mean = mean(ckd5, na.rm = TRUE), ckd5_sd = sd(ckd5, na.rm = TRUE)
  )

obs_model_full_g3 <- obs_complete_g3 %>%                                     # add the standardized (z-scored) version of every covariate, no tfr_z
  mutate(
    edu_z = standardize_value(edu, std_g3$edu_mean, std_g3$edu_sd),
    urban_z = standardize_value(urban, std_g3$urban_mean, std_g3$urban_sd),
    gdp_z = standardize_value(log_gdp, std_g3$log_gdp_mean, std_g3$log_gdp_sd),
    bmi_z = standardize_value(bmi, std_g3$bmi_mean, std_g3$bmi_sd),
    smoking_z = standardize_value(smoking, std_g3$smoking_mean, std_g3$smoking_sd),
    fpg_z = standardize_value(fpg, std_g3$fpg_mean, std_g3$fpg_sd),
    sbp_z = standardize_value(sbp, std_g3$sbp_mean, std_g3$sbp_sd),
    ldl_z = standardize_value(ldl, std_g3$ldl_mean, std_g3$ldl_sd),
    sodium_z = standardize_value(sodium, std_g3$sodium_mean, std_g3$sodium_sd),
    albuminuria_z = standardize_value(albuminuria, std_g3$albuminuria_mean, std_g3$albuminuria_sd),
    ckd5_z = standardize_value(ckd5, std_g3$ckd5_mean, std_g3$ckd5_sd)
  )

# --- 8b) Fit one model per age group, formula identical to the main model minus tfr_z ---
# 13 countries is a much thinner sample than the main 175 -- some age groups throw
# rank-deficiency/singular-fit warnings (fine, just noisier variance estimates from
# regions with only 1 country in them) but a tryCatch guards against a hard failure
# stopping the whole script if one age group's model can't be fit at all.

fit_lmer_safe_g3 <- function(formula, data, label) {
  tryCatch(
    lmer(formula, data = data, REML = FALSE,
        control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
    error = function(e) NULL
  )
}

models_by_age_g3 <- obs_model_full_g3 %>%
  group_by(age) %>%
  group_modify(~ {
    age_label <- unique(.y$age)
    age_rfs <- risk_factor_cols[vapply(
      .x %>% select(all_of(risk_factor_cols)),
      function(column) any(!is.na(column)),
      logical(1)
    )]

    required_cols <- c(base_required_cols_g3, age_rfs)
    obs_age <- .x %>% filter(complete.cases(across(all_of(required_cols))))

    if (nrow(obs_age) == 0) {
      return(tibble(model = list(NULL), age_rfs = list(character()), resid = list(tibble())))
    }

    age_rfs_z <- paste0(age_rfs, "_z")
    formula_terms <- c(                             # same terms as the main model, minus tfr_z
      "edu_z", "urban_z", "gdp_z",
      age_rfs_z,
      "urban_z:gdp_z", "sex", spline_term
    )
    model_formula <- as.formula(
      paste(
        "log_ratio_cardio ~",
        paste(formula_terms, collapse = " + "),
        "+ (1 | region / iso3) + (0 + year_c | iso3)"
      )
    )

    fitted_model <- fit_lmer_safe_g3(model_formula, obs_age, paste("Group3 age", age_label))
    if (is.null(fitted_model)) {
      return(tibble(model = list(NULL), age_rfs = list(age_rfs), resid = list(tibble())))
    }

    residuals_by_group <- obs_age %>%
      mutate(
        fitted_lmer = predict(
          fitted_model,
          newdata = obs_age,
          re.form = ~ (1 | region / iso3),
          allow.new.levels = TRUE
        ),
        resid_lmer = log_ratio_cardio - fitted_lmer
      ) %>%
      group_by(iso3, sex) %>%
      arrange(year, .by_group = TRUE) %>%
      summarise(
        last_resid = last(resid_lmer),
        last_train_year = max(year),
        .groups = "drop"
      ) %>%
      mutate(age = age_label)

    tibble(
      model = list(fitted_model),
      age_rfs = list(age_rfs),
      resid = list(residuals_by_group)
    )
  }) %>%
  ungroup()

# --- 8c) Project the 13 countries, 2024-2100, all 5 SSPs (mx_allcause from the same ltable_mx.rds/mx_proj used above) ---

proj_frame_base_g3 <- pop %>%
  filter(year > 2023, iso3 %in% group3_iso3) %>%
  mutate(
    pop = pop * 1000,
    pop = pmax(pop, 0.001)
  ) %>%
  left_join(covariates, by = c("scenario", "iso3", "age", "sex", "year")) %>%
  left_join(mx_proj, by = c("scenario", "iso3", "age", "sex", "year")) %>%       # same all-cause mortality source as the main pipeline
  left_join(risk_factors_g3, by = c("scenario", "iso3", "age", "sex", "year")) %>%  # projected risk factors, all 5 SSPs, group3's own file
  mutate(
    gdp = pmax(gdp, 1),
    mx_allcause = ifelse(mx_allcause == 0, NA_real_, mx_allcause),
    mx_allcause = pmax(mx_allcause, eps),
    year_c = year - 1990,
    region = make_un_subregion(iso3),
    log_gdp = log(gdp),
    edu_z = standardize_value(edu, std_g3$edu_mean, std_g3$edu_sd),          # standardize using the group3 sample's own mean/sd
    urban_z = standardize_value(urban, std_g3$urban_mean, std_g3$urban_sd),
    gdp_z = standardize_value(log_gdp, std_g3$log_gdp_mean, std_g3$log_gdp_sd),
    bmi_z = standardize_value(bmi, std_g3$bmi_mean, std_g3$bmi_sd),
    smoking_z = standardize_value(smoking, std_g3$smoking_mean, std_g3$smoking_sd),
    fpg_z = standardize_value(fpg, std_g3$fpg_mean, std_g3$fpg_sd),
    sbp_z = standardize_value(sbp, std_g3$sbp_mean, std_g3$sbp_sd),
    ldl_z = standardize_value(ldl, std_g3$ldl_mean, std_g3$ldl_sd),
    sodium_z = standardize_value(sodium, std_g3$sodium_mean, std_g3$sodium_sd),
    albuminuria_z = standardize_value(albuminuria, std_g3$albuminuria_mean, std_g3$albuminuria_sd),
    ckd5_z = standardize_value(ckd5, std_g3$ckd5_mean, std_g3$ckd5_sd)
  )

proj_spline_g3 <- proj_frame_base_g3 %>%
  group_by(age) %>%
  group_modify(~ {
    age_label <- unique(.y$age)
    model_row <- models_by_age_g3 %>% filter(age == age_label)

    if (nrow(model_row) == 0L || is.null(model_row$model[[1]])) {
      return(tibble())
    }

    age_rfs <- model_row$age_rfs[[1]]
    required_cols <- setdiff(c(proj_base_required_cols_g3, age_rfs), "mx_allcause")

    proj_age <- .x %>%
      filter(complete.cases(across(all_of(required_cols)))) %>%
      mutate(
        sex = factor(sex, levels = levels(obs_complete_g3$sex)),
        iso3 = factor(iso3, levels = levels(obs_complete_g3$iso3)),
        region = factor(region, levels = levels(obs_complete_g3$region))
      ) %>%
      filter(!is.na(sex), !is.na(iso3), !is.na(region))

    if (nrow(proj_age) == 0L) {
      return(tibble())
    }

    fitted_model <- model_row$model[[1]]
    proj_age <- proj_age %>%
      mutate(
        pred_log_m_cardio_base = predict(
          fitted_model,
          newdata = proj_age,
          re.form = ~ (1 | region / iso3),
          allow.new.levels = TRUE
        )
      )

    residuals_by_group <- model_row$resid[[1]]

    proj_age %>%
      left_join(residuals_by_group, by = c("iso3", "sex")) %>%
      mutate(
        last_resid = coalesce(last_resid, 0),
        last_train_year = coalesce(last_train_year, 2023L),
        years_ahead = year - last_train_year,
        resid_damped = last_resid * exp(-0.05 * years_ahead),
        pred_log_ratio = pred_log_m_cardio_base + resid_damped,
        pred_ratio_cardio = pmin(pmax(exp(pred_log_ratio), eps), 0.85),
        pred_allcause_deaths = mx_allcause * pop,
        pred_cardio_deaths = pred_ratio_cardio * pred_allcause_deaths
      ) %>%
      select(scenario, iso3, sex, year, pop, pred_ratio_cardio, pred_cardio_deaths, pred_allcause_deaths)
  }) %>%
  ungroup()

proj_spline_g3 <- proj_spline_g3[complete.cases(proj_spline_g3), ]

# -----------------------------------------------------------------------------
# 9) COMBINE MAIN (175) + GROUP3 (13) PROJECTIONS INTO ONE FINAL FILE
# -----------------------------------------------------------------------------

proj_spline <- bind_rows(                                       # overwrite proj_spline so it now covers all 188 countries
  proj_spline %>% mutate(iso3 = as.character(iso3), sex = as.character(sex)),
  proj_spline_g3 %>% mutate(iso3 = as.character(iso3), sex = as.character(sex))
)

saveRDS(proj_spline, file.path(output_dir, "proj_spline.rds"))   # re-save proj_spline.rds, now the one final 188-country projection file

projected_series <- proj_spline %>%                              # rebuild the projected half of the combined series (same recipe as chapter 7)
  mutate(
    iso3 = as.character(iso3),
    sex = as.character(sex),
    cardio_deaths = pred_cardio_deaths,
    allcause_deaths = pred_allcause_deaths,
    cardio_share = pred_ratio_cardio
  ) %>%
  select(scenario, iso3, age, sex, year, pop, cardio_deaths, allcause_deaths, cardio_share)

obs_and_proj <- bind_rows(observed_series, projected_series)     # observed_series (chapter 7) already covers all 204 historical countries

saveRDS(obs_and_proj, file.path(output_dir, "obs_and_proj.rds"))   # re-save obs_and_proj.rds, now the one final combined observed + projected file

# -----------------------------------------------------------------------------
# 10) PLOTS FOR THE 13 GROUP3 COUNTRIES
# -----------------------------------------------------------------------------
# Same plotting logic/output layout as Transcript/CountryPlots.R (country-level
# deaths + share, and country x age x sex deaths + share), restricted to the 13
# countries added in this chapter. The other 175 countries already have plots on
# disk from a prior CountryPlots.R run, so it isn't re-run here for all of them.

plots_out_dir           <- "/home/otavova/Healthy-SSPs/SSPs_CVD_mortality/Results/Plots/CountryPlots"  # same root as CountryPlots.R, so plots land side by side
country_plots_dir       <- file.path(plots_out_dir, "country_plots")
country_plots_share_dir <- file.path(plots_out_dir, "country_plots_share")
age_sex_dir             <- file.path(plots_out_dir, "country_age_sex_plots")
age_sex_share_dir       <- file.path(plots_out_dir, "country_age_sex_plots_share")
for (d in c(country_plots_dir, country_plots_share_dir, age_sex_dir, age_sex_share_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

age_order <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34",         # canonical age-band ordering for facet/filename sorting
              "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69",
              "70-74", "75-79", "80-84", "85-89", "90-94", "95+")

age_to_token <- function(age) {                                                    # "95+" -> "95_", "0-4" -> "0_4", for safe filenames
  if (grepl("\\+$", age)) sub("\\+$", "_", age) else gsub("-", "_", age)
}

ssp_colors_obs <- c(                                                                # fixed color mapping, matches CountryPlots.R exactly
  "Observed" = "black",
  "SSP1" = "#1B9E77", "SSP2" = "#D95F02", "SSP3" = "#7570B3",
  "SSP4" = "#E7298A", "SSP5" = "#66A61E"
)
legend_order <- c("Observed", "SSP1", "SSP2", "SSP3", "SSP4", "SSP5")

sanitize_filename <- function(iso) {                                                # ISO3 -> a safe filename stem based on the full country name
  name <- countrycode(iso, "iso3c", "country.name", warn = FALSE)
  if (is.na(name) || iso == "TWN") name <- if (iso == "TWN") "Taiwan" else iso
  gsub("[^A-Za-z0-9]+", "_", name)
}

obs_raw_g3 <- obs_and_proj %>%
  filter(scenario == "historical", iso3 %in% group3_iso3) %>%
  select(iso3, age, sex, year, death = cardio_deaths, allcause_death = allcause_deaths)

proj_raw_g3 <- obs_and_proj %>%
  filter(scenario != "historical", iso3 %in% group3_iso3) %>%
  select(scenario, iso3, age, sex, year, pop,
        pred_cardio_deaths = cardio_deaths, pred_allcause_deaths = allcause_deaths)

countries_plot_g3 <- sort(unique(proj_raw_g3$iso3))

last_obs_year_g3   <- max(obs_raw_g3$year)
first_proj_year_g3 <- min(proj_raw_g3$year)

# --- 10a) Country-level (both sexes combined): deaths + share ---

obs_country_g3 <- obs_raw_g3 %>%
  group_by(iso3, year) %>%
  summarise(cardio_deaths = sum(death, na.rm = TRUE),
           allcause_deaths = sum(allcause_death, na.rm = TRUE), .groups = "drop") %>%
  mutate(cardio_share_pct = 100 * cardio_deaths / allcause_deaths)

proj_country_g3 <- proj_raw_g3 %>%
  group_by(iso3, scenario, year) %>%
  summarise(cardio_deaths = sum(pred_cardio_deaths, na.rm = TRUE),
           allcause_deaths = sum(pred_allcause_deaths, na.rm = TRUE), .groups = "drop") %>%
  mutate(cardio_share_pct = 100 * cardio_deaths / allcause_deaths)

make_country_plot <- function(iso, value_col, y_label, out_path) {
  obs_i  <- obs_country_g3  %>% filter(iso3 == iso) %>% arrange(year)
  proj_i <- proj_country_g3 %>% filter(iso3 == iso) %>% arrange(scenario, year)
  if (nrow(obs_i) == 0 || nrow(proj_i) == 0) return(invisible(NULL))

  anchor_obs  <- obs_i[[value_col]][obs_i$year == last_obs_year_g3]
  anchor_proj <- mean(proj_i[[value_col]][proj_i$year == first_proj_year_g3], na.rm = TRUE)
  offset <- anchor_proj - anchor_obs

  obs_plot  <- obs_i  %>% transmute(year, series = "Observed", value = .data[[value_col]] + offset)
  proj_plot <- proj_i %>% transmute(year, series = scenario,   value = .data[[value_col]])

  plot_df <- bind_rows(obs_plot, proj_plot) %>% mutate(series = factor(series, levels = legend_order))

  country_name <- countrycode(iso, "iso3c", "country.name", warn = FALSE)
  if (is.na(country_name)) country_name <- iso

  proj_plot <- proj_plot %>% mutate(series = factor(series, levels = legend_order))

  p <- ggplot(plot_df, aes(x = year, y = value, color = series)) +
    geom_line(linewidth = 0.9) +
    geom_point(data = proj_plot, size = 1.3) +
    scale_color_manual(values = ssp_colors_obs, breaks = legend_order) +
    labs(title = country_name, x = "Year", y = y_label, color = NULL) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", plot.title = element_text(face = "plain"))

  ggsave(out_path, p, width = 7, height = 5, dpi = 150)
}

for (iso in countries_plot_g3) {
  make_country_plot(iso, "cardio_deaths", "CVD deaths",
                    file.path(country_plots_dir, paste0(sanitize_filename(iso), ".png")))
  make_country_plot(iso, "cardio_share_pct", "CVD share of all-cause deaths (%)",
                    file.path(country_plots_share_dir, paste0(sanitize_filename(iso), ".png")))
}

# --- 10b) Country x age x sex: deaths + share ---

group_keys_as <- c("iso3", "age", "sex")

obs_as_g3 <- obs_raw_g3 %>%
  select(iso3, age, sex, year, death) %>%
  rename(cardio_deaths = death)

proj_as_g3 <- proj_raw_g3 %>%
  select(scenario, iso3, age, sex, year, pred_cardio_deaths) %>%
  rename(cardio_deaths = pred_cardio_deaths)

anchor_obs_as <- obs_as_g3 %>%
  filter(year == last_obs_year_g3) %>%
  select(all_of(group_keys_as), anchor_obs = cardio_deaths)

anchor_proj_as <- proj_as_g3 %>%
  filter(year == first_proj_year_g3) %>%
  group_by(across(all_of(group_keys_as))) %>%
  summarise(anchor_proj = mean(cardio_deaths, na.rm = TRUE), .groups = "drop")

offsets_as <- anchor_obs_as %>%
  inner_join(anchor_proj_as, by = group_keys_as) %>%
  mutate(offset = anchor_proj - anchor_obs) %>%
  select(all_of(group_keys_as), offset)

obs_shifted_as <- obs_as_g3 %>%
  inner_join(offsets_as, by = group_keys_as) %>%
  mutate(series = "Observed", value = cardio_deaths + offset) %>%
  select(all_of(group_keys_as), year, series, value)

proj_series_as <- proj_as_g3 %>%
  mutate(series = scenario, value = cardio_deaths) %>%
  select(all_of(group_keys_as), year, series, value)

plot_df_as_g3 <- bind_rows(obs_shifted_as, proj_series_as) %>%
  mutate(series = factor(series, levels = legend_order))

obs_as_share_g3 <- obs_raw_g3 %>%
  filter(allcause_death > 0) %>%
  mutate(cardio_share_pct = 100 * death / allcause_death) %>%
  select(iso3, age, sex, year, cardio_share_pct)

proj_as_share_g3 <- proj_raw_g3 %>%
  filter(pred_allcause_deaths > 0) %>%
  mutate(cardio_share_pct = 100 * pred_cardio_deaths / pred_allcause_deaths) %>%
  select(scenario, iso3, age, sex, year, cardio_share_pct)

anchor_obs_share <- obs_as_share_g3 %>%
  filter(year == last_obs_year_g3) %>%
  select(all_of(group_keys_as), anchor_obs = cardio_share_pct)

anchor_proj_share <- proj_as_share_g3 %>%
  filter(year == first_proj_year_g3) %>%
  group_by(across(all_of(group_keys_as))) %>%
  summarise(anchor_proj = mean(cardio_share_pct, na.rm = TRUE), .groups = "drop")

offsets_share <- anchor_obs_share %>%
  inner_join(anchor_proj_share, by = group_keys_as) %>%
  mutate(offset = anchor_proj - anchor_obs) %>%
  select(all_of(group_keys_as), offset)

obs_shifted_share <- obs_as_share_g3 %>%
  inner_join(offsets_share, by = group_keys_as) %>%
  mutate(series = "Observed", value = cardio_share_pct + offset) %>%
  select(all_of(group_keys_as), year, series, value)

proj_series_share <- proj_as_share_g3 %>%
  mutate(series = scenario, value = cardio_share_pct) %>%
  select(all_of(group_keys_as), year, series, value)

plot_df_share_g3 <- bind_rows(obs_shifted_share, proj_series_share) %>%
  mutate(series = factor(series, levels = legend_order))

make_as_plot <- function(df_cell, title, y_label, out_path) {
  if (nrow(df_cell) == 0) return(invisible(NULL))
  proj_pts <- df_cell %>% filter(series != "Observed")
  p <- ggplot(df_cell, aes(x = year, y = value, color = series)) +
    geom_line(linewidth = 0.8) +
    geom_point(data = proj_pts, size = 1) +
    scale_color_manual(values = ssp_colors_obs, breaks = legend_order) +
    labs(title = title, x = "Year", y = y_label, color = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(size = 11))
  ggsave(out_path, p, width = 6, height = 4, dpi = 120)
}

n_saved_g3 <- 0
t0_g3 <- Sys.time()

for (cty in countries_plot_g3) {
  cty_dir_deaths <- file.path(age_sex_dir, sanitize_filename(cty))
  cty_dir_share  <- file.path(age_sex_share_dir, sanitize_filename(cty))
  dir.create(cty_dir_deaths, showWarnings = FALSE)
  dir.create(cty_dir_share,  showWarnings = FALSE)

  df_cty_deaths <- plot_df_as_g3    %>% filter(iso3 == cty)
  df_cty_share  <- plot_df_share_g3 %>% filter(iso3 == cty)
  ages_present  <- age_order[age_order %in% unique(as.character(df_cty_deaths$age))]

  for (ag in ages_present) {
    for (sx in c("Male", "Female")) {
      make_as_plot(df_cty_deaths %>% filter(age == ag, sex == sx),
                  paste0(sanitize_filename(cty), ": age ", ag, ", ", sx), "CVD deaths",
                  file.path(cty_dir_deaths, paste0(age_to_token(ag), "_", sx, ".png")))
      make_as_plot(df_cty_share %>% filter(age == ag, sex == sx),
                  paste0(sanitize_filename(cty), ": age ", ag, ", ", sx), "CVD share of all-cause deaths (%)",
                  file.path(cty_dir_share, paste0(age_to_token(ag), "_", sx, ".png")))
      n_saved_g3 <- n_saved_g3 + 2
    }
  }
}

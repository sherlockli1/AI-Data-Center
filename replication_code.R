# =====================================================================================
#  AI Data Center, Neighborhood Characteristics, and Drinking Water Contamination
#  Replication code (analyses run in R 4.5.2)
# =====================================================================================
#  Set the working directory to a folder containing these public data:
#    ICIS-AIR_downloads/            EPA ICIS-Air facility and program files   (echo.epa.gov/tools/data-downloads)
#    ECHO_EXPORTER.csv              EPA ECHO Exporter                         (same page)
#    SDWA_latest_downloads/         EPA Safe Drinking Water Information System (same page)
#    im3_open_source_data_center_atlas*.csv     IM3 Open Source Data Center Atlas (doi:10.57931/17054299)
#    data/registry/suedatacenters_20261110.csv  Compute Atlas registry (suedatacenters.org), 10 Nov 2026
#    data/transmission/*.shp        HIFLD transmission lines (doi:10.3886/E240591V1)
#    data/land_value/.../places_fmv_all.tif     Nolte (2020) land values (doi:10.1073/pnas.2012865117)
#  ACS (tidycensus; requires a Census API key) and NLCD (FedData) are downloaded and cached in data/.
#
#  Output: output/main (Tables 1-2, Figures 1-2), output/supplement (Tables S1-S3),
#  and the other numbers reported in the text, printed at the end.
# =====================================================================================

suppressPackageStartupMessages({
  library(terra); library(exactextractr); library(FedData); library(sf)
  library(tidyverse); library(lubridate); library(did); library(fixest); library(tidycensus)
})
options(tigris_use_cache = TRUE, timeout = 3600, readr.show_col_types = FALSE)

YEARS         <- 2000:2025   # drinking water panel
FIRST_ARRIVAL <- 2016        # EPA loaded its backlog into ICIS-Air in 2014-2015; earlier first records are not arrivals
AI_START      <- 2022
LAST_YEAR     <- 2026
STATES        <- c(setdiff(state.abb, c("AK", "HI")), "DC")
REG_FILE      <- "data/registry/suedatacenters_20261110.csv"
OUT_MAIN      <- "output/main"
OUT_SUPP      <- "output/supplement"
for (d in c("data", OUT_MAIN, OUT_SUPP)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

zscore <- function(x) as.numeric(scale(x))
num    <- function(x) prettyNum(x, big.mark = ",")
pct    <- function(x) sprintf("%.0f%%", 100 * mean(x))
ci     <- function(est, se) sprintf("%.2f (%.2f–%.2f)", exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se))
norm_county <- function(x) str_to_lower(str_squish(str_remove(x, regex(" (county|parish)$", ignore_case = TRUE))))

logit <- function(rhs, data, y, fe)
  feglm(as.formula(paste(y, "~", paste(rhs, collapse = " + "), "|", fe)), data,
        family = binomial, cluster = ~county)
or_table <- function(m, model)
  tibble(model = model, term = names(coef(m)), est = unname(coef(m)), se = unname(se(m)), n = nobs(m))
lincom <- function(m, a, b, model) {
  V <- vcov(m)
  tibble(model = model, term = a, est = coef(m)[[a]] + coef(m)[[b]],
         se = sqrt(V[a, a] + V[b, b] + 2 * V[a, b]), n = nobs(m))
}

LABELS <- c(z_share_hisp   = "% Hispanic",
            z_share_black  = "% Black",
            z_share_renter = "% renter-occupied",
            z_share_lep    = "% limited-English households",
            z_share_noncit = "% non-citizen",
            z_log_inc      = "Log median household income",
            z_share_ba     = "% with bachelor’s degree")

# Odds ratios (95% CI) per SD for the neighborhood characteristics, one column per model
or_wide <- function(long) {
  long <- filter(long, term %in% names(LABELS))
  body <- long |>
    mutate(cell = ci(est, se), term = factor(term, names(LABELS))) |>
    select(term, model, cell) |>
    pivot_wider(names_from = model, values_from = cell, values_fill = "") |>
    arrange(term) |>
    mutate(term = unname(LABELS[as.character(term)])) |>
    rename(`Characteristic (per SD)` = term)
  obs_row <- long |>
    distinct(model, n) |>
    mutate(n = num(n)) |>
    pivot_wider(names_from = model, values_from = n) |>
    mutate(`Characteristic (per SD)` = "Observations", .before = 1)
  bind_rows(body, obs_row) |> select(`Characteristic (per SD)`, any_of(unique(long$model)))
}

# NLCD land cover shares (per 10 percentage points), downloaded and cached by county
nlcd_shares <- function(polys, id_col, year, prefix, dir) {
  dir.create(file.path(dir, "by_county"), recursive = TRUE, showWarnings = FALSE)
  one <- function(cty) {
    f <- file.path(dir, "by_county", paste0(prefix, "_", cty, ".rds"))
    for (i in 1:3) if (!file.exists(f)) try({
      pp <- polys[polys$county == cty, ]
      r  <- get_nlcd(template = pp, label = paste0(prefix, "_", cty), year = year, dataset = "landcover",
                     landmass = "L48", extraction.dir = dir, force.redo = i > 1)
      levels(r) <- NULL
      x <- exact_extract(r, st_transform(pp, terra::crs(r)), "frac", progress = FALSE)
      x$id <- pp[[id_col]]
      saveRDS(x, f)
    })
    if (file.exists(f)) readRDS(f)
  }
  lc  <- map(unique(polys$county), one) |> compact() |> bind_rows() |>
    mutate(across(starts_with("frac_"), \(x) replace_na(x, 0)))
  cls <- function(codes) rowSums(lc[, intersect(paste0("frac_", codes), names(lc)), drop = FALSE])
  tibble(id = lc$id, dev_high = 10 * cls(24), dev_any = 10 * cls(21:24),
         open_land = 10 * cls(c(31, 52, 71, 81, 82)))
}

theme_paper <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
save_table <- function(df, file, dir = OUT_MAIN) write_csv(df, file.path(dir, file), na = "")
save_fig   <- function(p, file, w, h, dir = OUT_MAIN)
  ggsave(file.path(dir, file), p, width = w, height = h, dpi = 300, bg = "white")


# ===== 1. Data =======================================================================

# Data centers in EPA ICIS-Air: NAICS 518210, SIC 7374, or a data center name or operator.
# Arrival = year of first ICIS-Air record.
site_file <- "data/data_centers_epa_icis.csv"
if (!file.exists(site_file)) {
  read_icis <- function(file, ...) read_csv(file.path("ICIS-AIR_downloads", file), ...,
                                            col_types = cols(.default = "c"), locale = locale(encoding = "latin1"))
  name_pat <- regex(paste("DATA ?CENTER|DATA CTR|EQUINIX|\\bQTS\\b|DIGITAL REALTY|CYRUSONE|VANTAGE DATA",
                          "ALIGNED DATA|STACK INFRASTRUCTURE|SABEY|CLOUDHQ|CORESITE|DUPONT FABROS",
                          "RAGINGWIRE|NTT GLOBAL DATA|COMPASS DATA|EDGECORE|IRON MOUNTAIN DATA",
                          "AMAZON DATA|VADATA|AMAZON WEB SERVICES|COREWEAVE", sep = "|"), ignore_case = TRUE)
  fac <- read_icis("ICIS-AIR_FACILITIES.csv") |>
    filter(str_detect(coalesce(NAICS_CODES, ""), "\\b518210\\b") |
             str_detect(coalesce(SIC_CODES, ""), "\\b7374\\b") |
             str_detect(coalesce(FACILITY_NAME, ""), name_pat))
  first_rec <- read_icis("ICIS-AIR_PROGRAMS.csv", col_select = c(PGM_SYS_ID, BEGIN_DATE)) |>
    filter(PGM_SYS_ID %in% fac$PGM_SYS_ID) |>
    mutate(year = as.integer(year(mdy(BEGIN_DATE)))) |>
    filter(!is.na(year)) |>
    group_by(PGM_SYS_ID) |>
    summarise(first_year = min(year), .groups = "drop")
  echo_file <- list.files(".", "^ECHO_EXPORTER\\.csv$", recursive = TRUE, ignore.case = TRUE, full.names = TRUE)[1]
  coords <- read_csv(echo_file, col_select = c(REGISTRY_ID, FAC_LAT, FAC_LONG), col_types = cols(.default = "c")) |>
    distinct(REGISTRY_ID, .keep_all = TRUE)
  fac |>
    inner_join(first_rec, by = "PGM_SYS_ID") |>
    filter(first_year <= LAST_YEAR) |>                                   # later dates are typos
    left_join(coords, by = "REGISTRY_ID", na_matches = "never") |>
    mutate(site_id = coalesce(REGISTRY_ID, PGM_SYS_ID)) |>
    group_by(site_id) |>
    summarise(name = first(FACILITY_NAME), state = first(STATE), county = first(COUNTY_NAME),
              first_year = min(first_year), lat = as.numeric(first(FAC_LAT)), lon = as.numeric(first(FAC_LONG)),
              .groups = "drop") |>
    mutate(lon = if_else(lon > 0, -lon, lon)) |>
    write_csv(site_file)
}
dc_sites <- read_csv(site_file, col_types = cols(first_year = "i", lat = "d", lon = "d", .default = "c")) |>
  mutate(county_name = norm_county(county), legacy = first_year < FIRST_ARRIVAL)   # legacy: arrival year unknown
dc_pts <- dc_sites |>
  filter(between(lat, 24, 50), between(lon, -125, -66)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(5070)

# Coverage check: IM3 atlas data centers within 1 km of an EPA site
im3 <- list.files(".", "^im3_open_source_data_center_atlas.*\\.csv$", recursive = TRUE, full.names = TRUE)[1] |>
  read_csv(col_types = cols(lon = "d", lat = "d", .default = "c")) |>
  filter(between(lat, 24, 50), between(lon, -125, -66)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(5070)
im3$recorded <- lengths(st_is_within_distance(im3, dc_pts, dist = 1000)) > 0

# SDWIS: community water systems and health-based violations
read_sdwa <- function(file, vars) read_csv(file.path("SDWA_latest_downloads", file), col_select = all_of(vars),
                                           col_types = cols(.default = "c"))
systems <- read_sdwa("SDWA_PUB_WATER_SYSTEMS.csv", c("PWSID", "PWS_TYPE_CODE", "STATE_CODE")) |>
  filter(PWS_TYPE_CODE == "CWS") |>
  distinct(PWSID, .keep_all = TRUE) |>
  transmute(PWSID, state = STATE_CODE) |>
  left_join(read_sdwa("SDWA_GEOGRAPHIC_AREAS.csv", c("PWSID", "AREA_TYPE_CODE", "COUNTY_SERVED")) |>
              filter(AREA_TYPE_CODE == "CN", !is.na(COUNTY_SERVED)) |>
              transmute(PWSID, county = norm_county(COUNTY_SERVED)) |>
              distinct(PWSID, .keep_all = TRUE),                        # multi-county systems: first county
            by = "PWSID")

viol_file <- "data/sdwa_health_based_violations.rds"
if (!file.exists(viol_file)) {
  read_sdwa("SDWA_VIOLATIONS_ENFORCEMENT.csv",
            c("PWSID", "VIOLATION_ID", "NON_COMPL_PER_BEGIN_DATE", "IS_HEALTH_BASED_IND")) |>
    filter(IS_HEALTH_BASED_IND == "Y") |>
    distinct(PWSID, VIOLATION_ID, .keep_all = TRUE) |>                 # one record per violation
    mutate(year = as.integer(year(mdy(NON_COMPL_PER_BEGIN_DATE)))) |>
    filter(year %in% YEARS) |>
    distinct(PWSID, year) |>
    saveRDS(viol_file)
}
viol <- readRDS(viol_file)

# Census tracts (2010 boundaries): ACS 2008-2012 and 2015-2019
acs12_file <- "data/acs2012_tracts.rds"
if (!file.exists(acs12_file)) {
  v12 <- c(pop = "B01003_001", med_inc = "B19013_001", race_d = "B03002_001", blk = "B03002_004",
           hisp = "B03002_012", edu_d = "B15003_001", ba = "B15003_022", ma = "B15003_023",
           prof = "B15003_024", phd = "B15003_025")
  saveRDS(map(STATES, \(s) get_acs(geography = "tract", variables = v12, year = 2012, state = s,
                                    geometry = TRUE, output = "wide")) |> bind_rows(), acs12_file)
}
tracts_geo <- readRDS(acs12_file) |> st_transform(5070)
tracts12 <- tracts_geo |>
  st_drop_geometry() |>
  transmute(GEOID, county = substr(GEOID, 1, 5), pop = popE,
            log_inc = log(med_incE), share_black = blkE / race_dE, share_hisp = hispE / race_dE,
            share_ba = (baE + maE + profE + phdE) / edu_dE)
tracts_geo <- tracts_geo |> select(GEOID) |> mutate(county = substr(GEOID, 1, 5))

voice_file <- "data/acs2012_voice_all_states.rds"
if (!file.exists(voice_file)) {
  vv <- c(cit_d = "B05001_001", noncit = "B05001_006", hh_lang = "B16002_001", lep1 = "B16002_004",
          lep2 = "B16002_007", lep3 = "B16002_010", lep4 = "B16002_013", ten_d = "B25003_001",
          renter = "B25003_003")
  saveRDS(map_dfr(STATES, \(s) get_acs(geography = "tract", variables = vv, year = 2012, state = s,
                                        output = "wide")), voice_file)
}
voice12 <- readRDS(voice_file) |>
  transmute(GEOID, share_noncit = noncitE / cit_dE, share_renter = renterE / ten_dE,
            share_lep = (lep1E + lep2E + lep3E + lep4E) / hh_langE)

acs19_file <- "data/acs2019_all_states.rds"
if (!file.exists(acs19_file)) {
  va <- c(pop = "B01003_001", race_d = "B03002_001", blk = "B03002_004", hisp = "B03002_012",
          ten_d = "B25003_001", renter = "B25003_003", med_inc = "B19013_001", edu_d = "B15003_001",
          ba = "B15003_022", ma = "B15003_023", prof = "B15003_024", phd = "B15003_025",
          hh_lang = "C16002_001", lep1 = "C16002_004", lep2 = "C16002_007", lep3 = "C16002_010",
          lep4 = "C16002_013")
  vb <- c(units = "B25024_001", sfd = "B25024_002", u5_9 = "B25024_006", u10_19 = "B25024_007",
          u20_49 = "B25024_008", u50 = "B25024_009", mobile = "B25024_010", hu = "B25002_001",
          vacant = "B25002_003", yr_built = "B25035_001")
  saveRDS(map_dfr(STATES, \(s)
    get_acs(geography = "tract", variables = va, year = 2019, state = s, output = "wide") |>
      left_join(get_acs(geography = "tract", variables = vb, year = 2019, state = s, output = "wide") |>
                  select(-NAME), by = "GEOID")), acs19_file)
}
acs19 <- readRDS(acs19_file) |>
  transmute(GEOID, pop = popE, race_d = race_dE, blk = blkE, hisp = hispE, ten_d = ten_dE,
            renter = renterE, med_inc = med_incE, hh_lang = hh_langE, lep = lep1E + lep2E + lep3E + lep4E,
            share_ba = (baE + maE + profE + phdE) / edu_dE,
            share_sfd = sfdE / unitsE, share_multi5 = (u5_9E + u10_19E + u20_49E + u50E) / unitsE,
            share_mobile = mobileE / unitsE, vacancy = vacantE / huE,
            yr_built = if_else(yr_builtE >= 1900, yr_builtE, NA_real_))

# Siting needs: transmission lines (HIFLD) and land value (Nolte 2020, ln $/ha)
tx <- st_read(list.files("data/transmission", "\\.shp$", full.names = TRUE, recursive = TRUE)[1], quiet = TRUE) |>
  st_transform(5070) |>
  mutate(kv = as.numeric(VOLTAGE))                                     # unknown voltage is coded -999999
tx_counts <- c(nrow(tx), sum(tx$kv < 0, na.rm = TRUE))
tx115 <- filter(tx, kv >= 115); tx345 <- filter(tx, kv >= 345); rm(tx)
km_to <- function(pts, lines)
  as.numeric(st_distance(pts, lines[st_nearest_feature(pts, lines), ], by_element = TRUE)) / 1000
lv_file <- list.files("data/land_value", "^places_fmv_all\\.tif$", full.names = TRUE, recursive = TRUE)
lv_r <- rast(lv_file[!str_detect(lv_file, "__MACOSX")][1])

ctrl_file <- "data/tract_siting_controls.rds"
if (!file.exists(ctrl_file)) {
  pts <- st_point_on_surface(st_geometry(tracts_geo))
  saveRDS(tibble(GEOID = tracts_geo$GEOID,
                 area_km2 = as.numeric(st_area(tracts_geo)) / 1e6,
                 log_dist_tx115 = log1p(km_to(pts, tx115)),
                 log_dist_tx345 = log1p(km_to(pts, tx345)),
                 log_land_value = exact_extract(lv_r, st_transform(tracts_geo, terra::crs(lv_r)), "mean",
                                                progress = FALSE)), ctrl_file)
}
tract_ctrl <- readRDS(ctrl_file)

# Compute Atlas registry of proposed data centers (snapshot used in the paper)
stopifnot(tools::md5sum(REG_FILE) == "90cf9204980c119989bccb5424e3abda")
reg_all <- read_csv(REG_FILE, guess_max = 5000)
reg <- reg_all |>
  filter(facility_type == "data_center") |>
  mutate(ann_year_num = as.integer(substr(announced_date, 1, 4)),
         latitude = as.numeric(latitude), longitude = as.numeric(longitude),
         is_proposal = coalesce(status %in% c("proposed", "permitted", "under_construction", "cancelled") |
                                  (status == "operational" & coalesce(ann_year_num, 0L) >= 2020), FALSE),
         is_legacy = coalesce(status == "operational", FALSE) & !is_proposal)


# ===== 2. Drinking water: staggered difference-in-differences =========================
# Treated from the first data center arrival in the county; counties with a pre-existing site are dropped
g_cty <- dc_sites |>
  group_by(state, county_name) |>
  summarise(g = min(first_year), legacy = any(legacy), .groups = "drop")
water <- systems |>
  filter(!is.na(county)) |>
  left_join(g_cty, by = c("state", "county" = "county_name")) |>
  filter(!coalesce(legacy, FALSE)) |>
  mutate(g = replace_na(g, 0L), g = if_else(g > max(YEARS), 0L, g)) |>   # arrivals after 2025 serve as controls
  filter(g == 0 | g > min(YEARS)) |>
  select(PWSID, state, county, g) |>
  crossing(year = YEARS) |>
  left_join(mutate(viol, y = 1L), by = c("PWSID", "year")) |>
  mutate(y = replace_na(y, 0L),
         id = as.integer(factor(PWSID)), cty = as.integer(factor(paste(state, county))))

set.seed(2026)                                                         # bootstrap confidence bands
att <- att_gt(yname = "y", tname = "year", idname = "id", gname = "g", data = water,
              control_group = "notyettreated", clustervars = c("id", "cty"), est_method = "dr",
              print_details = FALSE)
es  <- aggte(att, type = "dynamic", min_e = -5, max_e = 5, na.rm = TRUE)
ov  <- aggte(att, type = "simple", na.rm = TRUE)

fig1 <- with(es, tibble(e = egt, att = att.egt, lo = att.egt - crit.val.egt * se.egt,
                        hi = att.egt + crit.val.egt * se.egt)) |>
  mutate(when = if_else(e < 0, "Before", "After")) |>
  ggplot(aes(e, 100 * att, colour = when)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = 100 * lo, ymax = 100 * hi), size = 0.3) +
  scale_colour_manual(values = c(Before = "grey50", After = "#2C6FBB"), guide = "none") +
  scale_x_continuous(breaks = -5:5) +
  labs(x = "Years since data center arrival", y = "Effect (percentage points)") +
  theme_paper


# ===== 3. Siting of built data centers ================================================
dc_tract <- dc_pts |>
  st_join(select(tracts_geo, GEOID), join = st_within) |>
  st_drop_geometry() |>
  filter(!is.na(GEOID)) |>
  group_by(GEOID) |>
  summarise(first_dc = min(first_year), .groups = "drop")
risk <- tracts12 |>
  filter(pop > 0) |>
  left_join(dc_tract, by = "GEOID") |>
  filter(is.na(first_dc) | first_dc >= FIRST_ARRIVAL)                 # drop tracts with a pre-existing site
panel_b <- bind_rows(
  risk |> mutate(ai = 0L, host = as.integer(!is.na(first_dc) & first_dc < AI_START)),
  risk |> filter(is.na(first_dc) | first_dc >= AI_START) |> mutate(ai = 1L, host = as.integer(!is.na(first_dc)))
) |>
  group_by(county, ai) |>
  filter(any(host == 1)) |>
  ungroup()

land11 <- nlcd_shares(filter(tracts_geo, county %in% panel_b$county), "GEOID", 2011, "cty", "data/nlcd")
built <- panel_b |>
  left_join(tract_ctrl, by = "GEOID") |>
  left_join(land11, by = c("GEOID" = "id")) |>
  left_join(voice12, by = "GEOID") |>
  mutate(log_density = log1p(pop / area_km2)) |>
  filter(if_all(c(log_inc, share_hisp, share_black, share_ba, share_renter, share_noncit, share_lep,
                  log_density, log_dist_tx115, log_dist_tx345, dev_high, dev_any, open_land,
                  log_land_value), is.finite)) |>
  mutate(across(c(log_inc, share_hisp, share_black, share_ba, share_renter, share_noncit, share_lep,
                  log_land_value), zscore, .names = "z_{.col}"))

COMP  <- c("z_share_hisp", "z_share_black", "z_log_inc", "z_share_ba")
VOICE <- c("z_share_renter", "z_share_noncit", "z_share_lep")
NEEDS <- c("log_density", "log_dist_tx115", "log_dist_tx345", "dev_high", "dev_any", "open_land", "z_log_land_value")
rhs_b <- list("A. Neighborhood characteristics and density" = c(COMP, "log_density"),
              "B. + Power, land cover, land value"          = c(COMP, NEEDS),
              "C. + Renter, non-citizen, limited English"   = c(COMP, VOICE, NEEDS))
fits_b <- map(rhs_b, logit, data = built, y = "host", fe = "county^ai")

# Model C with each neighborhood characteristic interacted with the AI era
lab_all <- paste0("Built, ", FIRST_ARRIVAL, "–", LAST_YEAR)
lab_pre <- paste0("Built, ", FIRST_ARRIVAL, "–", AI_START - 1)
lab_ai  <- paste0("Built, ", AI_START, "–", LAST_YEAR)
m_per  <- logit(c(sprintf("(%s) * ai", paste(c(COMP, VOICE), collapse = " + ")), NEEDS), built, "host", "county^ai")
or_per <- bind_rows(or_table(m_per, lab_pre) |> filter(term %in% c(COMP, VOICE)),
                    map_dfr(c(COMP, VOICE), \(v) lincom(m_per, v, paste0(v, ":ai"), lab_ai)))


# ===== 4. Siting of proposed data centers =============================================
reg_tr <- reg |>
  filter(!is.na(latitude), !is.na(longitude)) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(5070) |>
  st_join(select(tracts_geo, GEOID), join = st_within) |>
  st_drop_geometry() |>
  filter(!is.na(GEOID)) |>
  group_by(GEOID) |>
  summarise(n_prop = sum(is_proposal), legacy_dc = as.integer(any(is_legacy)), .groups = "drop")
prop_counties <- unique(substr(reg_tr$GEOID[reg_tr$n_prop > 0], 1, 5))

land19 <- nlcd_shares(filter(tracts_geo, county %in% prop_counties), "GEOID", 2019, "s5", "data/nlcd2019")
prop <- acs19 |>
  mutate(county = substr(GEOID, 1, 5)) |>
  filter(county %in% prop_counties, pop >= 100) |>
  left_join(tract_ctrl, by = "GEOID") |>
  left_join(land19, by = c("GEOID" = "id")) |>
  left_join(reg_tr, by = "GEOID") |>
  mutate(across(c(n_prop, legacy_dc), \(x) replace_na(x, 0L)),
         proposal = as.integer(n_prop > 0),
         share_hisp = hisp / race_d, share_black = blk / race_d, share_renter = renter / ten_d,
         share_lep = lep / hh_lang, log_inc = log(med_inc), log_density = log1p(pop / area_km2)) |>
  filter(if_all(c(share_hisp, share_black, share_renter, share_lep, log_inc, share_ba, log_density,
                  log_dist_tx115, log_dist_tx345, log_land_value, dev_high, open_land), is.finite)) |>
  mutate(across(c(share_hisp, share_black, share_renter, share_lep, log_inc, share_ba, log_land_value,
                  share_sfd, share_multi5, share_mobile, vacancy, yr_built), zscore, .names = "z_{.col}"))

NEEDS_P <- c("log_density", "log_dist_tx115", "log_dist_tx345", "dev_high", "open_land", "z_log_land_value", "legacy_dc")
HOUSING <- c("z_share_sfd", "z_share_multi5", "z_share_mobile", "z_vacancy", "z_yr_built")
rhs_p <- list("A. Neighborhood characteristics and density"            = c(COMP, "log_density"),
              "B. + Power, land, land value, pre-existing data center" = c(COMP, NEEDS_P),
              "C. + Renter, limited English"                           = c(COMP, "z_share_renter", "z_share_lep", NEEDS_P),
              "D. + Housing type, vacancy, age"                        = c(COMP, "z_share_renter", "z_share_lep", NEEDS_P, HOUSING))
fits_p <- map(rhs_p, logit, data = prop, y = "proposal", fe = "county")


# ===== 5. Cancellation of proposed data centers =======================================
proj <- reg |>
  filter(is_proposal, !is.na(latitude), !is.na(longitude)) |>
  mutate(canceled   = as.integer(status == "cancelled"),
         ann_year   = if_else(is.na(ann_year_num), "unknown", as.character(ann_year_num)),
         mw         = as.numeric(capacity_mw_planned),
         acres      = as.numeric(land_acres),
         size_known = case_when(!is.na(mw) ~ "mw", !is.na(acres) ~ "acres", TRUE ~ "neither"),
         log_mw     = if_else(is.na(mw) | mw <= 0, 0, log(mw)),
         log_acres  = if_else(is.na(acres) | acres <= 0, 0, log(acres))) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) |>
  st_transform(5070) |>
  st_join(select(tracts_geo, GEOID), join = st_within) |>
  filter(!is.na(GEOID)) |>
  mutate(county = substr(GEOID, 1, 5), state_fips = substr(GEOID, 1, 2))

# Neighborhood characteristics within 3 km: tracts weighted by their share of area inside the buffer
buf3 <- st_buffer(select(proj, id, county), 3000)
near_tracts <- tracts_geo[lengths(st_intersects(tracts_geo, buf3)) > 0, ] |>
  select(GEOID) |>
  inner_join(acs19, by = "GEOID") |>
  mutate(t_area = as.numeric(st_area(geometry)))
comp3 <- suppressWarnings(st_intersection(select(buf3, id), near_tracts)) |>
  mutate(w = as.numeric(st_area(geometry)) / t_area) |>
  st_drop_geometry() |>
  group_by(id) |>
  summarise(pop_3km      = sum(w * pop),
            share_hisp   = sum(w * hisp) / sum(w * race_d),
            share_black  = sum(w * blk) / sum(w * race_d),
            share_renter = sum(w * renter) / sum(w * ten_d),
            share_lep    = sum(w * lep) / sum(w * hh_lang),
            log_inc      = weighted.mean(log(med_inc), w * pop, na.rm = TRUE), .groups = "drop") |>
  mutate(log_density_3km = log1p(pop_3km / (pi * 3^2)))

proj$log_dist_tx115 <- log1p(km_to(st_geometry(proj), tx115))
proj$log_dist_tx345 <- log1p(km_to(st_geometry(proj), tx345))
proj$log_land_value <- exact_extract(lv_r, st_transform(buf3, terra::crs(lv_r)), "mean", progress = FALSE)
land_buf <- nlcd_shares(buf3, "id", 2019, "reg19", "data/nlcd2019")

canc <- proj |>
  st_drop_geometry() |>
  left_join(comp3, by = "id") |>
  left_join(land_buf, by = "id") |>
  filter(pop_3km >= 100,
         if_all(c(share_hisp, share_renter, log_inc, share_black, share_lep, log_density_3km,
                  log_dist_tx115, log_dist_tx345, dev_high, open_land, log_land_value), is.finite)) |>
  mutate(across(c(share_hisp, share_renter, log_inc, share_black, share_lep), zscore, .names = "z_{.col}"))

RHS_C <- c("z_share_hisp", "z_share_renter", "z_log_inc", "z_share_black", "z_share_lep",
           "log_density_3km", "log_dist_tx115", "log_dist_tx345", "dev_high", "open_land",
           "log_land_value", "i(size_known)", "log_mw", "log_acres")
FE_C   <- "state_fips + ann_year"
m_canc <- logit(RHS_C, canc, "canceled", FE_C)
m_h4   <- logit(c(RHS_C, "z_share_hisp:z_share_lep"), canc, "canceled", FE_C)

# H1-H4 (one-sided p values; Holm adjustment across H1-H3)
tests <- tibble(Hypothesis = c("H1", "H2", "H3", "H4"),
                term = c("z_share_hisp", "z_share_renter", "z_log_inc", "z_share_hisp:z_share_lep"),
                `Predictor (per SD)` = c("% Hispanic", "% renter-occupied", "Log median household income",
                                         "% Hispanic × % limited-English households"),
                fewer = c(TRUE, TRUE, FALSE, TRUE),
                fit = list(m_canc, m_canc, m_canc, m_h4)) |>
  mutate(est = map2_dbl(fit, term, \(m, t) coef(m)[[t]]),
         se  = map2_dbl(fit, term, \(m, t) se(m)[[t]]),
         p   = if_else(fewer, pnorm(est / se), 1 - pnorm(est / se)))
tests$p_holm <- c(p.adjust(tests$p[1:3], method = "holm"), NA)


# ===== 6. Tables and figures ==========================================================
# Estimation samples (fixed-effect groups without variation in the outcome drop out)
built_m <- built[fixest::obs(fits_b[[3]]), ]
prop_m  <- prop[fixest::obs(fits_p[[3]]), ]
canc_m  <- canc[fixest::obs(m_canc), ]

describe <- function(d, flag, yes, no) d |>
  mutate(group = if_else(.data[[flag]] == 1, yes, no)) |>
  group_by(group) |>
  summarise(`Tracts, n` = num(n()),
            `Median household income, $` = num(round(median(exp(log_inc)), -2)),
            `Hispanic, %` = sprintf("%.1f", 100 * mean(share_hisp)),
            `Black, %` = sprintf("%.1f", 100 * mean(share_black)),
            `Renter-occupied housing, %` = sprintf("%.1f", 100 * mean(share_renter)),
            `Limited-English households, %` = sprintf("%.1f", 100 * mean(share_lep)),
            `Bachelor’s degree or higher, %` = sprintf("%.1f", 100 * mean(share_ba)),
            `Population density, persons/km² (median)` = num(round(median(expm1(log_density)))),
            `Distance to ≥115 kV line, km (median)` = sprintf("%.1f", median(expm1(log_dist_tx115))),
            `High-intensity developed land, %` = sprintf("%.1f", 10 * mean(dev_high))) |>
  pivot_longer(-group, names_to = "Characteristic") |>
  pivot_wider(names_from = group, values_from = value) |>
  select(Characteristic, all_of(c(yes, no)))
table1 <- left_join(describe(built_m, "host", "Built: with data centers", "Built: without data centers"),
                    describe(prop_m, "proposal", "Proposed: with proposed data centers", "Proposed: without proposed data centers"),
                    by = "Characteristic")

t2_long <- bind_rows(or_table(fits_b[[3]], lab_all), or_per,
                     or_table(fits_p[[3]], "Proposed"), or_table(m_canc, "Cancelled"))
events <- c(sum(built_m$host), sum(built_m$host[built_m$ai == 0]), sum(built_m$host[built_m$ai == 1]),
            sum(prop_m$proposal), sum(canc_m$canceled))
names(events) <- c(lab_all, lab_pre, lab_ai, "Proposed", "Cancelled")
table2 <- bind_rows(or_wide(t2_long),
                    as_tibble(as.list(num(events))) |> mutate(`Characteristic (per SD)` = "Events", .before = 1))

tableS3 <- tests |>
  transmute(Hypothesis, `Predictor (per SD)`,
            `Predicted direction` = if_else(fewer, "Fewer cancellations", "More cancellations"),
            `Odds ratio (95% CI)` = ci(est, se),
            p                     = sprintf("%.3f", p),
            `adjusted p`          = if_else(is.na(p_holm), "", sprintf("%.3f", p_holm)))

fig2 <- t2_long |>
  filter(term %in% c("z_share_hisp", "z_share_renter")) |>
  mutate(stage = factor(model, rev(names(events))), variable = LABELS[term]) |>
  ggplot(aes(exp(est), stage)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = exp(est - 1.96 * se), xmax = exp(est + 1.96 * se)), colour = "#2C6FBB") +
  facet_wrap(~ variable) +
  scale_x_log10() +
  labs(x = "Odds ratio per 1 SD, within county (log scale)", y = NULL,
       caption = "Cancelled: an odds ratio below 1 means projects in such areas were cancelled less often.") +
  theme_paper

save_table(table1, "Table1_descriptives.csv")
save_table(table2, "Table2_siting_stages.csv")
save_fig(fig1, "Figure1_water_event_study.png", 4.6, 4)
save_fig(fig2, "Figure2_siting_stages.png", 9, 4)
save_table(or_wide(imap_dfr(fits_b, or_table)), "TableS1_built_models.csv", OUT_SUPP)
save_table(or_wide(imap_dfr(fits_p, or_table)), "TableS2_proposal_models.csv", OUT_SUPP)
save_table(tableS3, "TableS3_cancellation_tests.csv", OUT_SUPP)


# ===== 7. Other numbers reported in the text ==========================================
arr  <- filter(dc_sites, !legacy)
yrs  <- table(arr$first_year)
top3 <- head(sort(table(arr$state), decreasing = TRUE), 3)
intext <- c(
  "EPA data centers (with contiguous-US coordinates)" = sprintf("%d (%d)", nrow(dc_sites), nrow(dc_pts)),
  "  pre-existing / arriving 2016-2026"           = sprintf("%d / %d", sum(dc_sites$legacy), nrow(arr)),
  "  arriving 2016-2021 / 2022-2026"              = sprintf("%d / %d", sum(arr$first_year < AI_START), sum(arr$first_year >= AI_START)),
  "  arrivals by year"                            = paste(names(yrs), yrs, sep = ": ", collapse = ", "),
  "  top states"                                  = paste(names(top3), top3, sep = ": ", collapse = ", "),
  "IM3 data centers (within 1 km of an EPA site)" = sprintf("%d (%d, %s)", nrow(im3), sum(im3$recorded), pct(im3$recorded)),
  "Transmission segments (unknown voltage)"       = sprintf("%d (%d)", tx_counts[1], tx_counts[2]),
  "Registry: crypto mining excluded"              = sum(reg_all$facility_type == "crypto_mining"),
  "Proposed data centers"                         = sum(reg$is_proposal),
  "  with announcement date (all / cancelled)"    = with(reg, sprintf("%s / %s", pct(!is.na(ann_year_num[is_proposal])),
                                                                     pct(!is.na(ann_year_num[status %in% "cancelled"])))),
  "Cancellation sample (cancelled)"               = sprintf("%d (%d, %.1f%%)", nrow(canc), sum(canc$canceled), 100 * mean(canc$canceled)),
  "Water systems in treated counties"             = n_distinct(water$id[water$g > 0]),
  "Effect on any health-based violation, pp (95% CI)" =
    sprintf("%.2f (%.2f to %.2f)", 100 * ov$overall.att, 100 * (ov$overall.att - 1.96 * ov$overall.se),
            100 * (ov$overall.att + 1.96 * ov$overall.se)))
cat(sprintf("%-50s %s", names(intext), intext), sep = "\n")


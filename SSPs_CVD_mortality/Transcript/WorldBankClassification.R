# =============================================================================
# WorldBankClassification.R
# =============================================================================
# Purpose : Region lookup table -- six geographic World Bank regions + a
#           High-Income bucket pulled out of geography, covering all 204
#           countries modeled in this project (Data/cardio.rds).
#           Updated 08/09/2026 to the World Bank's CURRENT country
#           classification (https://api.worldbank.org/v2/country). Formerly
#           named Sellers_regions.R, based on Sellers et al. (2020) Table S1
#           (10.1007/s10584-020-02824-4), which used the World Bank's 2019
#           scheme -- renamed since the table no longer reflects that 2019
#           snapshot.
#
#           Two things changed relative to the original Sellers table:
#           1. High-Income membership is rebuilt from the WB's current
#              income classification (incomeLevel == "High income"), not
#              the 2019-vintage fixed list -- e.g. Bulgaria, Romania, Russia,
#              Costa Rica and Guyana are now High-Income; Argentina no longer
#              is. 20 countries missing from the original table (mostly
#              small territories) are now included.
#           2. The World Bank renamed "Middle East & North Africa" to
#              "Middle East, North Africa, Afghanistan & Pakistan" and moved
#              Afghanistan and Pakistan into it out of South Asia -- this
#              table follows that current membership and name.
#
#           4 project countries aren't classified by the World Bank at all
#           (Taiwan, Cook Islands, Niue, Tokelau) and are assigned manually
#           below (see comments at each iso3).
# Usage   : source("WorldBankClassification.R")  ->  adds `world_bank_classification` data frame
#           then: left_join(your_data, world_bank_classification, by = "iso3")
# =============================================================================

world_bank_classification <- data.frame(stringsAsFactors = FALSE,

  iso3 = c(
    # ── East Asia & Pacific ──────────────────────────────────────────────────
    "CHN","COK","FJI","FSM","IDN","KHM","KIR","LAO","MHL","MMR",
    "MNG","MYS","NIU","PHL","PNG","PRK","SLB","THA","TKL","TLS",
    "TON","TUV","VNM","VUT","WSM",
    # COK, NIU, TKL: small Pacific territories not classified by the World
    # Bank; grouped here geographically, not treated as High-Income.

    # ── Europe & Central Asia ────────────────────────────────────────────────
    "ALB","ARM","AZE","BIH","BLR","GEO","KAZ","KGZ","MDA","MKD",
    "MNE","SRB","TJK","TKM","TUR","UKR","UZB",

    # ── Latin America & Caribbean ────────────────────────────────────────────
    "ARG","BLZ","BOL","BRA","COL","CUB","DMA","DOM","ECU","GRD",
    "GTM","HND","HTI","JAM","LCA","MEX","NIC","PER","PRY","SLV",
    "SUR","VCT","VEN",

    # ── Middle East, North Africa, Afghanistan & Pakistan ────────────────────
    # (World Bank's current name for this region; AFG and PAK moved here out
    # of South Asia along with the rename)
    "AFG","DJI","DZA","EGY","IRN","IRQ","JOR","LBN","LBY","MAR",
    "PAK","PSE","SYR","TUN","YEM",

    # ── South Asia ───────────────────────────────────────────────────────────
    "BGD","BTN","IND","LKA","MDV","NPL",

    # ── Sub-Saharan Africa ───────────────────────────────────────────────────
    "AGO","BDI","BEN","BFA","BWA","CAF","CIV","CMR","COD","COG",
    "COM","CPV","ERI","ETH","GAB","GHA","GIN","GMB","GNB","GNQ",
    "KEN","LBR","LSO","MDG","MLI","MOZ","MRT","MUS","MWI","NAM",
    "NER","NGA","RWA","SDN","SEN","SLE","SOM","SSD","STP","SWZ",
    "TCD","TGO","TZA","UGA","ZAF","ZMB","ZWE",

    # ── High-Income ───────────────────────────────────────────────────────────
    "AND","ARE","ASM","ATG","AUS","AUT","BEL","BGR","BHR","BHS",
    "BMU","BRB","BRN","CAN","CHE","CHL","CRI","CYP","CZE","DEU",
    "DNK","ESP","EST","FIN","FRA","GBR","GRC","GRL","GUM","GUY",
    "HRV","HUN","IRL","ISL","ISR","ITA","JPN","KNA","KOR","KWT",
    "LTU","LUX","LVA","MCO","MLT","MNP","NLD","NOR","NRU","NZL",
    "OMN","PAN","PLW","POL","PRI","PRT","QAT","ROU","RUS","SAU",
    "SGP","SMR","SVK","SVN","SWE","SYC","TTO","TWN","URY","USA",
    "VIR"
    # TWN: not classified by the World Bank (not a member); treated as
    # High-Income to match its actual economic status.
  ),

  region = c(
    rep("East Asia & Pacific",                                25),
    rep("Europe & Central Asia",                               17),
    rep("Latin America & Caribbean",                           23),
    rep("Middle East, North Africa, Afghanistan & Pakistan",   15),
    rep("South Asia",                                           6),
    rep("Sub-Saharan Africa",                                  47),
    rep("High-Income",                                         71)
  )
)

print(table(world_bank_classification$region))

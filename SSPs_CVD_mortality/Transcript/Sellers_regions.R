# =============================================================================
# Sellers_regions.R
# =============================================================================
# Purpose : Region lookup table based on Sellers (2020) Table S1.
#           Seven regions: six geographic (World Bank 2019) + High-Income.
#           Source: 10.1007/s10584-020-02824-4, ESM Table S1.
# Usage   : source("Sellers_regions.R")  ->  adds `sellers_regions` data frame
#           then: left_join(your_data, sellers_regions, by = "iso3")
# =============================================================================

sellers_regions <- data.frame(stringsAsFactors = FALSE,

  iso3 = c(
    # ── East Asia & Pacific ──────────────────────────────────────────────────
    "KHM","MNG","THA","CHN","MMR","TLS","FJI","PRK","TON",
    "IDN","PNG","VUT","LAO","PHL","VNM","MYS","WSM","FSM","SLB",

    # ── Europe & Central Asia ────────────────────────────────────────────────
    "ALB","KAZ","RUS","ARM","XKX","SRB","AZE","KGZ","TJK",
    "BLR","MKD","TUR","BIH","MDA","TKM","BGR","MNE","UKR",
    "GEO","ROU","UZB",

    # ── Latin America & Caribbean ────────────────────────────────────────────
    "BLZ","SLV","NIC","BOL","GRD","PRY","BRA","GTM","PER",
    "COL","GUY","LCA","CRI","HTI","VCT","CUB","HND","SUR",
    "DOM","JAM","VEN","ECU","MEX",

    # ── Middle East & North Africa ───────────────────────────────────────────
    "DZA","JOR","PSE","DJI","LBN","SYR","EGY","LBY","TUN",
    "IRN","MAR","YEM","IRQ",

    # ── South Asia ───────────────────────────────────────────────────────────
    "AFG","IND","PAK","BGD","MDV","LKA","BTN","NPL",

    # ── Sub-Saharan Africa ───────────────────────────────────────────────────
    "AGO","GMB","RWA","BEN","GHA","STP","BWA","GIN","SEN",
    "BFA","GNB","SLE","BDI","KEN","SOM","CPV","LSO","ZAF",
    "CAF","LBR","SSD","TCD","MDG","SDN","COM","MWI","SWZ",
    "COG","MLI","TZA","CIV","MRT","TGO","COD","MUS","UGA",
    "GNQ","MOZ","ZMB","ERI","NAM","ZWE","ETH","NER","GAB","NGA",

    # ── High-Income Countries ────────────────────────────────────────────────
    "ARG","HKG","PRI","AUS","HUN","QAT","AUT","ISL","SAU",
    "BHS","IRL","SYC","BHR","ISR","SGP","BRB","ITA","SVK",
    "BEL","JPN","SVN","BRN","KWT","KOR","CAN","LVA","ESP",
    "CHL","LTU","SWE","HRV","LUX","CHE","CYP","MLT","TWN",
    "CZE","NLD","TTO","DNK","NZL","ARE","EST","NOR","GBR",
    "FIN","OMN","USA","FRA","PAN","URY","DEU","POL","GRC","PRT"
  ),

  region = c(
    rep("East Asia & Pacific",        19),
    rep("Europe & Central Asia",       21),
    rep("Latin America & Caribbean",   23),
    rep("Middle East & North Africa",  13),
    rep("South Asia",                   8),
    rep("Sub-Saharan Africa",          46),
    rep("High-Income",                 55)
  )
)

print(table(sellers_regions$region))

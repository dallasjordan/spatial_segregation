# Script for additional analyses, per Scott's suggestions on improvements to manuscript
  # KDE Area
  # Centroid
  # Boxplots of lat/lon by species and colony
  # Oceanographic analysis
# Dallas Jordan Aug 7 2021


# Setup -------------------------------------------------------------------

# Load in KDEs and lat/lon data

  # loading lat/lon files
    
    # for pc
    setwd("E:/project_data/spatial_segregation/data")
    
    load("LAALdata_midway_withTrackID.Rdata")
    LAALmid <- LAAL
    LAALmid$id <- paste0("lm",LAALmid$track)
    load("LAALdata_tern_withTrackID.Rdata")
    LAALtern <- LAAL
    LAALtern$id <- paste0("lt",LAALtern$track)
    
    LAAL <- rbind(LAALmid, LAALtern)
    
    load("BFALdata_midway_withTrackID.Rdata")
    BFALmid <- BFAL
    BFALmid$id <- paste0("bm",BFALmid$track)
    load("BFALdata_tern_withTrackID.Rdata")
    BFALtern <- BFAL
    BFALtern$id <- paste0("bt",BFALtern$track)
    
    BFAL <- rbind(BFALmid, BFALtern)
    
    # comparisons: allLAAL v allBFAL
    #              ternLAAL v midwayLAAL
    #              ternBFAL v midwayBFAL
    #              ternLAAL v ternBFAL
    #              midwayLAAL v midwayBFAL
    all_data <- rbind(LAAL,BFAL)
    lm <- all_data[grep("lm", all_data$id), ]
    lt <- all_data[grep("lt", all_data$id), ]
    bm <- all_data[grep("bm", all_data$id), ]
    bt <- all_data[grep("bt", all_data$id), ]
    
  # Load in KDE contours to calculate area
    
    # for pc - You can load in the UDs and make contours or load in the contours
    allLAALud <- load("E:/project_data/spatial_segregation/data/final_ud/allLAAL.Rdata")
    allBFALud <- load("E:/project_data/spatial_segregation/data/final_ud/allBFAL.Rdata")
    midLAALud <- load("E:/project_data/spatial_segregation/data/final_ud/midLAAL.Rdata")
    ternLAALud <- load("E:/project_data/spatial_segregation/data/final_ud/midLAAL.Rdata")
    
    # for pc - loading in the contours
    setwd("E:/project_data/spatial_segregation/figures/individual/midLAAL/master_script_contours/")
    load("vert95_midLAAL.Rdata")
    load("vert50_midLAAL.Rdata")
    load("vert10_midLAAL.Rdata")

    
    plot(getverticeshr(allBFAL,percent=50), add=T)
    vert50_allBFAL <- getverticeshr(allBFAL,percent=50)

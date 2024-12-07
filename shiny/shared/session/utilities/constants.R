refTypes <- c('genome', 'spike_in')
insertTypes <- c('subnucleosomal', 'nucleosomal')
centerTypes <- c("rpba", "subnucleosomal_fraction")
paColors <- list(
    # BROWN  = rgb(0.2, 0,   0),   # extreme gain    
    # RED    = rgb(0.9, 0,   0),   # 3, CN3, full gain
    # YELLOW = rgb(0.8, 0.8, 0),   # 2.5, mosaic state
    # GREEN  = rgb(0,   0.8, 0),   # 2, CN2, no CNV
    # CYAN   = rgb(0,   0.6, 1),   # 1.5, mosaic state
    # BLUE   = rgb(0,   0,   1),   # 1, CN1, full loss
    # PURPLE = rgb(0.4, 0,   0.4), # extreme loss
    # BLACK  = rgb(0,   0,   0)    # 0, nothing    
    RED     = rgb(0.9, 0,   0),
    GREY    = rgb(0.75, 0.75, 0.75),
    BLUE    = rgb(0,   0,   1)
)
stageColors <- c(
    early_round = CONSTANTS$plotlyColors$black,
    late_round  = CONSTANTS$plotlyColors$blue,
    early_elong = CONSTANTS$plotlyColors$orange,
    int_elong   = CONSTANTS$plotlyColors$green,
    late_elong  = CONSTANTS$plotlyColors$purple
)

# 1,24290X4,day35-wt-rs-rep1,early_round,Early Round Spermatids,1,#440154FF
# 2,24290X6,day35-wt-rs-rep2,early_round,Early Round Spermatids,2,#440154FF
# 3,24290X8,day38-wt-rs-rep1,late_round,Late Round Spermatids ,1,magenta
# 4,24290X10,day38-wt-rs-rep2,late_round,Late Round Spermatids ,2,magenta
# 5,24290X3,day32-wt-es,early_elong,Early Elongating Spermatids,1,#21908CFF
# 6,24290X2,day34-wt-es,early_elong,Early Elongating Spermatids,2,#21908CFF
# 7,24290X5,day35-wt-es-rep1,int_elong,Intermediate Elongating Spermatids,1,#5DC863FF
# 8,24290X7,day35-wt-es-rep2,int_elong,Intermediate Elongating Spermatids,2,#5DC863FF
# 9,24290X9,day38-wt-es-rep1,late_elong,Late Elongating Spermatids ,1,#FDE725FF
# 10,24290X11,day38-wt-es-rep2,late_elong,Late Elongating Spermatids ,2,#FDE725FF

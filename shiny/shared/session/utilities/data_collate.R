#----------------------------------------------------------------------
# data recovery from atac/collate data packages
#----------------------------------------------------------------------
paCollate_ttl <- CONSTANTS$ttl$month
paCollate_force <- FALSE
paCollate_create <- "asNeeded"
paCollate_loadPersistent <- function(..., sep = "\t", header = FALSE, force = NULL, spinnerMessage = NULL){
    if (!is.null(spinnerMessage)) startSpinner(session, message = spinnerMessage)
    if (is.null(force)) force <- paCollate_force
    filePath <- loadPersistentFile(..., sep = sep, header = header, force = force, ttl = paCollate_ttl)
    persistentCache[[filePath]]$data
}
paCollate_getCached <- function(..., create = NULL, spinnerMessage = NULL){
    if (!is.null(spinnerMessage)) startSpinner(session, message = spinnerMessage)
    if (is.null(create)) create <- paCollate_create
    paCollateCache$get(..., permanent = TRUE, create = create)$value
}

# break a data package into its component object types and return the requested type object
paCollate_break_file <- function(sourceId, type){
    expandSourceFilePath(sourceId, paste0("collate-", type, ".rds"))
}   
paCollate_get_component <- function(sourceId, type){
    typeFile <- paCollate_break_file(sourceId, type)
    if(file.exists(typeFile)){
        readRDS(typeFile)
    } else {
        startSpinner(session, message = "loading collate data")
        x <- readRDS(getSourceFilePath(sourceId, "collate"))
        for(type_ in names(x)){
            startSpinner(session, message = paste("breaking", type_))
            filePath <- paCollate_break_file(sourceId, type_)
            if(!file.exists(filePath)) saveRDS(x[[type_]], file = filePath)
        }
        x[[type]]
    }
}

# load smaller collate components into RAM
paCollate_load_ram <- function(sourceId, type){
    paCollate_getCached(
        type,
        key = sourceId,
        from = 'ram',
        createFn = function(...) {
            paCollate_get_component(sourceId, type)
        },
        spinnerMessage = paste("loading", type)
    )
}
paCollate_load_ram_reactive <- function(sourceId, type) reactive({
    sourceId <- sourceId()
    req(sourceId)
    paCollate_load_ram(sourceId, type)
})

#----------------------------------------------------------------------
# getters for collate components
#----------------------------------------------------------------------
# parsed environment variables
paCollate_env <- function(sourceId) paCollate_load_ram(sourceId, "env")

# samples sorted by staging_order
paCollate_samples <- function(sourceId) paCollate_load_ram(sourceId, "samples")[order(staging_order)]
paCollate_samples_reactive <- function(sourceId) reactive({
    sourceId <- sourceId()
    req(sourceId)
    paCollate_samples(sourceId)
})

# genome references
paCollate_refs <- function(sourceId) paCollate_load_ram(sourceId, "references")

# bins (without associate count data)
paCollate_bins <- function(sourceId) paCollate_load_ram(sourceId, "bins")

# gc bias models, i.e., negative binomial regresssion fits
paCollate_gc_bias_models <- function(sourceId) paCollate_loadPersistent(
    sourceId = sourceId,
    contentFileType = "gcBiasModels",
    spinnerMessage = "loading gc bias models"
)

# bin parameters by sample name
paCollate_gc_bin_smp <- function(sourceId, sample_name) paCollate_load_ram(sourceId, "gc_bin_smp")[, sample_name]
paCollate_n_ins_wgt_smp <- function(sourceId, sample_name) paCollate_load_ram(sourceId, "n_ins_wgt_smp")[[sample_name]]

#----------------------------------------------------------------------
# handle reads per bin (rpb) normalization patterns
#----------------------------------------------------------------------
paCollate_rpb_smp <- function(sourceId, normalizeTn5Site){
    paCollate_getCached(
        "rpb_smp",
        keyObject = list(sourceId, normalizeTn5Site),
        from = 'ram',
        # create = "once",
        createFn = function(...) {

            # collect the counts or weights per bin per sample, based on normalizeTn5Site
            n_obs_bin_smp <- paCollate_get_component(sourceId, "n_obs_bin_smp")
            rpb_smp <- if(normalizeTn5Site){
                paCollate_get_component(sourceId, "n_wgt_bin_smp")
            } else {
                n_obs_bin_smp
            }

            # increase the counts of poorly mappable bins
            mpp_bin_smp <- paCollate_get_component(sourceId, "mpp_bin_smp")
            rpb_smp <- ifelse(rpb_smp > 0 & mpp_bin_smp > 0, rpb_smp / mpp_bin_smp, 0)

            # rescale counts to sum to the number of observed reads
            # this is done to maintain the same statistical weight pre- and post-normalization
            env  <- paCollate_env(sourceId)
            bins <- paCollate_bins(sourceId)
            gc_bin <- rowMeans(paCollate_load_ram(sourceId, "gc_bin_smp"))
            I_ref  <- getGenomeBins(bins, env$PRIMARY_GENOME)
            I_norm <- getIncludedAutosomeBins(bins, gc_bin, env$PRIMARY_GENOME)
            for(j in 1:ncol(rpb_smp)){
                n_obs  <- sum(n_obs_bin_smp[I_norm, j])
                n_norm <- sum(rpb_smp[I_norm, j])
                rpb_smp[I_ref, j] <- rpb_smp[I_ref, j] / n_norm * n_obs
            }
            rpb_smp
        },
        spinnerMessage = "loading rpb_smp"
    )
}

# obj <- list(
#     env = env[c(
#         'PRIMARY_GENOME',
#         'SPIKE_IN_GENOME',
#         'GENOME',
#         'MAPPABILITY_SIZE_LEVELS',
#         'MIN_INSERT_SIZE',
#         'MAX_INSERT_SIZE',
#         'BIN_SIZE'
#     )],
#     samples    = samples,
#     references = refs,
#     bins = bins[, .(
#         chrom = chr, start0 = bs0, # end1 = be1, 
#         included = included #, nAlleles = nAlleles
#     )],
#     n_obs_bin_smp = n_obs_bin_smp,
#     n_wgt_bin_smp = n_wgt_bin_smp,
#     mpp_bin_smp   = mpp_bin_smp,
#     gc_bin_smp    = gc_bin_smp,
#     f_obs_ref_isl_smp   = f_obs_ref_isl_smp,
#     n_ref_wgt_is_gc_smp = n_ref_wgt_is_gc_smp
#     n_ins_wgt_smp = n_ins_wgt_smp
# )

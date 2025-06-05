# data recover from atac/site data packages
paTn5Site_ttl <- CONSTANTS$ttl$month
paTn5Site_force <- FALSE
paTn5Site_create <- "asNeeded"
paTn5Sites_loadPersistent <- function(..., sep = "\t", header = FALSE, force = NULL, spinnerMessage = NULL){
    # dmsg(spinnerMessage)
    if (!is.null(spinnerMessage)) startSpinner(session, message = spinnerMessage)
    if (is.null(force)) force <- paTn5Site_force
    filePath <- loadPersistentFile(..., sep = sep, header = header, force = force, ttl = paTn5Site_ttl)
    persistentCache[[filePath]]$data
}
paTn5Sites_getCached <- function(..., create = NULL, spinnerMessage = NULL){
    # dmsg(spinnerMessage)
    if (!is.null(spinnerMessage)) startSpinner(session, message = spinnerMessage)
    if (is.null(create)) create <- paTn5Site_create
    paTn5SitesCache$get(..., permanent = TRUE, create = create)$value
}

# helper list for paTn5Sites data access
paTn5Sites_cfg <- list(
    prm_obs = list(
        type = "primary-sites-observed"
    ),
    prm_exp = list(
        type = "primary-sites-expected"
    ),
    spk_obs = list(
        type = "spike-in-sites-observed"
    ),
    spk_exp = list(
        type = "spike-in-sites-expected"
    )
)

# samples sorted by staging_order
paTn5Sites_samples <- function(sourceId) paTn5Sites_loadPersistent(
    sourceId = sourceId,
    contentFileType = "samples-metadata-file",
    sep = ",",
    header = TRUE
)[order(staging_order)]
paTn5Sites_samples_reactive <- function(sourceId) reactive({
    sourceId <- sourceId()
    req(sourceId)
    paTn5Sites_samples(sourceId)
})

# all possible Tn5 kmers, with expansion to local sequence context
paTn5Sites_kmers <- function(sourceId) paTn5Sites_loadPersistent(
    sourceId = sourceId,
    contentFileType = "tn5-kmers-file",
    spinnerMessage = "loading Tn5 kmers"
)[[1]]
paTn5Sites_expand <- function(kmers) {
    unname(sapply(kmers, function(kmer) {
        kmer_ <- chartr("ACGT", "TGCA", kmer)
        kmer  <- strsplit(kmer,  "")[[1]]
        kmer_ <- strsplit(kmer_, "")[[1]]
        paste(
            paste(c(kmer [1:2], "|", kmer [3], "nn", kmer [4:6], "nn", kmer [7], "-", kmer [8:9]), collapse = ""),
            paste(c(kmer_[1:2], "-", kmer_[3], "nn", kmer_[4:6], "nn", kmer_[7], "|", kmer_[8:9]), collapse = ""),
            sep = "<br>"
        )
    }))
}
paTn5Sites_expand_simple <- function(kmers) {
    unname(sapply(kmers, function(kmer) {
        kmer <- strsplit(kmer, "")[[1]]
        paste(c(kmer [1:3], "nn", kmer [4:6], "nn", kmer [7:9]), collapse = "")
    }))
}
paTn5Sites_kmer_hasDinuc <- function(sourceId, dinuc = "CG", I = TRUE, n = 2) {
    paTn5Sites_getCached(
        "kmer_hasDinuc",
        keyObject = list(sourceId, dinuc, n),
        from = 'disk',
        createFn = function(...){
            kmers <- paTn5Sites_kmers(sourceId)
            stringr::str_count(paTn5Sites_expand_simple(kmers), dinuc) >= n
        },
        spinnerMessage = paste("finding dinuc", dinuc)
    )[I]
}

# count recovery, by kmer, insert size level(, and sample)
paTn5Sites_n_prm_obs_isl_smps <- function(sourceId) paTn5Sites_loadPersistent(
    sourceId = sourceId,
    contentFileType = paTn5Sites_cfg$prm_obs$type,
    spinnerMessage = "loading n_prm_smps"
)
paTn5Sites_n_prm_obs_isl_smp <- function(sourceId, sampleName) paTn5Sites_getCached(
    "n_prm_obs_isl_smp",
    keyObject = list(sourceId, sampleName),
    from = 'disk',
    createFn = function(...){
        paTn5Sites_n_prm_obs_isl_smps(sourceId)[,,sampleName]
    },
    spinnerMessage = "loading n_prm_smp"
)
paTn5Sites_n_prm_obs_isl_set <- function(sourceId, sampleIs, setId) paTn5Sites_getCached(
    "n_prm_obs_isl_set",
    keyObject = list(sourceId, sampleIs),
    from = 'disk',
    # create = "asNeeded",
    createFn = function(...){
        paTn5Sites_n_prm_obs_isl_smps(sourceId)[,,sampleIs] %>% 
        apply(c(1, 2), sum)
    },
    spinnerMessage = paste("loading n_prm_obs", setId)
)
paTn5Sites_n_prm_obs_isl_all <- function(sourceId) paTn5Sites_getCached(
    "n_prm_obs_isl_all",
    key = sourceId,
    from = 'disk',
    create = "asNeeded",
    createFn = function(...){
        # int[kmer, insert size level, sample] 1:262144, 1:20, 1:28
        paTn5Sites_n_prm_obs_isl_smps(sourceId) %>% 
        apply(c(1, 2), sum) # int[kmer, insert size level]
    },
    spinnerMessage = "loading n_prm_obs"
)
paTn5Sites_n_prm_exp_isl <- function(sourceId) paTn5Sites_loadPersistent(
    sourceId = sourceId,
    contentFileType = paTn5Sites_cfg$prm_exp$type,
    spinnerMessage = "loading n_prm_exp"
)

# weight insert size levels by summed counts across all samples
paTn5Sites_w_prm_obs_isl_smp <- function(sourceId, sampleName) paTn5Sites_getCached(
    "w_prm_obs_isl_smp",
    keyObject = list(sourceId, sampleName),
    createFn = function(...){
        paTn5Sites_n_prm_obs_isl_smp(sourceId, sampleName) %>% colSums
    },
    spinnerMessage = "loading weights"
)
paTn5Sites_w_prm_obs_isl_set <- function(sourceId, sampleIs, setId) paTn5Sites_getCached(
    "w_prm_obs_isl_set",
    keyObject = list(sourceId, sampleIs),
    createFn = function(...){
        paTn5Sites_n_prm_obs_isl_set(sourceId, sampleIs, setId) %>% colSums
    },
    spinnerMessage = paste("loading weights", setId)
)
paTn5Sites_w_prm_obs_isl_all <- function(sourceId) paTn5Sites_getCached(
    "w_prm_obs_isl_all",
    key = sourceId,
    createFn = function(...){
        paTn5Sites_n_prm_obs_isl_all(sourceId) %>% colSums
    },
    spinnerMessage = "loading weights"
)

# convert counts to frequencies by insert size level
paTn5Sites_n_to_f_isl <- function(n_isl){
    sapply(1:ncol(n_isl), function(islI) n_isl[, islI] / sum(n_isl[, islI]))
}
paTn5Sites_f_prm_obs_isl_smp <- function(sourceId, sampleName) paTn5Sites_getCached(
    "f_prm_obs_isl_smp",
    keyObject = list(sourceId, sampleName),
    from = 'disk',
    createFn = function(...){
        paTn5Sites_n_prm_obs_isl_smp(sourceId, sampleName) %>% paTn5Sites_n_to_f_isl()
    },
    spinnerMessage = "loading f_prm_obs"
)
paTn5Sites_f_prm_obs_isl_set <- function(sourceId, sampleIs, setId) paTn5Sites_getCached(
    "f_prm_obs_isl_set",
    keyObject = list(sourceId, sampleIs),
    from = 'disk',
    createFn = function(...){
        paTn5Sites_n_prm_obs_isl_set(sourceId, sampleIs, setId) %>% paTn5Sites_n_to_f_isl()
    },
    spinnerMessage = paste("loading f_prm_obs", setId)
)
paTn5Sites_f_prm_obs_isl_all <- function(sourceId) paTn5Sites_getCached(
    "f_prm_obs_isl_all",
    key = sourceId,
    from = 'disk',
    createFn = function(...){
        paTn5Sites_n_prm_obs_isl_all(sourceId) %>% paTn5Sites_n_to_f_isl()
    },
    spinnerMessage = "loading f_prm_obs"
)
paTn5Sites_f_prm_exp_isl <- function(sourceId) paTn5Sites_getCached(
    "f_prm_exp_isl",
    key = sourceId,
    from = 'disk',
    createFn = function(...){
        paTn5Sites_n_prm_exp_isl(sourceId) %>% paTn5Sites_n_to_f_isl()
    },
    spinnerMessage = "loading f_prm_exp"
)

# aggregate frequencies weighted across insert size levels
paTn5Sites_f_prm_obs_smp <- function(sourceId, sampleName) paTn5Sites_getCached(
    "f_prm_obs_smp",
    keyObject = list(sourceId, sampleName),
    from = 'disk',
    createFn = function(...){
        w <- paTn5Sites_w_prm_obs_isl_smp(sourceId, sampleName)
        f <- paTn5Sites_f_prm_obs_isl_smp(sourceId, sampleName) %>% 
             apply(1, weighted.mean, w)
        f / sum(f) # enforce normalization to 1 (in fact, it already is...)
    }, 
    spinnerMessage = "aggregating f_prm_obs"
)
paTn5Sites_f_prm_obs_set <- function(sourceId, sampleIs, setId) paTn5Sites_getCached(
    "f_prm_obs_set",
    keyObject = list(sourceId, sampleIs),
    from = 'disk',
    createFn = function(...){
        w <- paTn5Sites_w_prm_obs_isl_set(sourceId, sampleIs, setId)
        f <- paTn5Sites_f_prm_obs_isl_set(sourceId, sampleIs, setId) %>% 
             apply(1, weighted.mean, w)
        f / sum(f) # enforce normalization to 1 (in fact, it already is...)
    }, 
    spinnerMessage = paste("aggregating f_prm_obs", setId)
)
paTn5Sites_f_prm_obs_all <- function(sourceId) paTn5Sites_getCached(
    "f_prm_obs_all",
    key = sourceId,
    from = 'disk',
    createFn = function(...){
        w <- paTn5Sites_w_prm_obs_isl_all(sourceId)
        f <- paTn5Sites_f_prm_obs_isl_all(sourceId) %>% 
             apply(1, weighted.mean, w)
        f / sum(f) # enforce normalization to 1 (in fact, it already is...)
    }, 
    spinnerMessage = "aggregating f_prm_obs"
)
paTn5Sites_f_prm_exp <- function(sourceId) paTn5Sites_getCached(
    "f_prm_exp",
    key = sourceId,
    from = 'disk',
    createFn = function(...){
        w <- paTn5Sites_w_prm_obs_isl_all(sourceId)
        f <- paTn5Sites_f_prm_exp_isl(sourceId) %>% 
             apply(1, weighted.mean, w)
        f / sum(f) # enforce normalization to 1 (in fact, it already is...)
    }, 
    spinnerMessage = "aggregating f_prm_obs"
)

# aggregate frequencies to fold enrichment, weighted by insert size levels
paTn5Sites_e_prm_obs_set <- function(sourceId, sampleIs, setId) paTn5Sites_getCached(
    "e_prm_obs_set",
    keyObject = list(sourceId, sampleIs),
    from = 'disk',
    createFn = function(...){
        w <- paTn5Sites_w_prm_obs_isl_set(sourceId, sampleIs, setId)
        (
            paTn5Sites_f_prm_obs_isl_set(sourceId, sampleIs, setId) / 
            paTn5Sites_f_prm_exp_isl(    sourceId)
        ) %>% apply(1, weighted.mean, w)
    }, 
    spinnerMessage = paste("aggregating f_prm_set", setId)
)
paTn5Sites_e_prm_obs_all <- function(sourceId) paTn5Sites_getCached(
    "e_prm_obs_all",
    key = sourceId,
    from = 'disk',
    createFn = function(...){
        w <- paTn5Sites_w_prm_obs_isl_all(sourceId)
        (
            paTn5Sites_f_prm_obs_isl_all(sourceId) / 
            paTn5Sites_f_prm_exp_isl(    sourceId)
        ) %>% apply(1, weighted.mean, w)
    }, 
    spinnerMessage = "aggregating f_prm_exp"
)

# sort Tn5 kmers by either fold enrichment or frequency
paTn5Sites_order <- function(sourceId, by = "enrichment") {
    req(by %in% c("enrichment", "frequency"))
    if (by == "enrichment") {
        x <- paTn5Sites_e_prm_obs_all(sourceId)
    } else {
        x <- paTn5Sites_f_prm_obs_all(sourceId)
    }
    order(-x)
}

# tabulate all sites
paTn5Sites_table <- function(sourceId, by) paTn5Sites_getCached(
    "paTn5Sites_table",
    keyObject = list(sourceId, by),
    from = 'disk',
    createFn = function(...){
        order <- paTn5Sites_order(sourceId, by)
        kmers <- paTn5Sites_kmers(sourceId)[order]
        f_obs <- paTn5Sites_f_prm_obs_all(sourceId)[order]
        f_exp <- paTn5Sites_f_prm_exp(sourceId)[order]
        list(
            order = order,
            dt = data.table(
                rank    = 1:length(kmers),
                kmer    = kmers,
                e       = paTn5Sites_e_prm_obs_all(sourceId)[order],
                f_obs   = f_obs,
                cdf_obs = cumsum(f_obs),
                f_exp   = f_exp,
                cdf_exp = cumsum(f_exp)
            )
        )
    }
)

# emprically determined sample sets staging_order (thus, row I) to yield:
paTn5Sites_sampleIs_RS <- 3:11  # round spermatids with nucleosomal signal
paTn5Sites_sampleIs_ES <- 15:24 # elongating spermatids with high accessibility signal
paTn5Sites_tn5SetData <- function(sourceId) paTn5Sites_getCached(
    "tn5SetData",
    keyObject = list(sourceId, paTn5Sites_sampleIs_RS, paTn5Sites_sampleIs_ES),
    from = 'disk',
    # create = "once",
    createFn = function(...){
        f_rs <- paTn5Sites_f_prm_obs_set(sourceId, paTn5Sites_sampleIs_RS, "RS") %>% log10()
        f_es <- paTn5Sites_f_prm_obs_set(sourceId, paTn5Sites_sampleIs_ES, "ES") %>% log10()
        e_rs <- paTn5Sites_e_prm_obs_set(sourceId, paTn5Sites_sampleIs_RS, "RS") %>% log10()
        e_es <- paTn5Sites_e_prm_obs_set(sourceId, paTn5Sites_sampleIs_ES, "ES") %>% log10()
        has_CG <- paTn5Sites_kmer_hasDinuc(sourceId, "CG")
        has_GC <- paTn5Sites_kmer_hasDinuc(sourceId, "GC")
        list(
            sourceId = sourceId,
            sampleIs_RS = paTn5Sites_sampleIs_RS,
            sampleIs_ES = paTn5Sites_sampleIs_ES,
            f_rs = f_rs,
            f_es = f_es,
            e_rs = e_rs,
            e_es = e_es,
            has_dinuc = list(
                CG = has_CG,
                GC = has_GC,
                neither = !has_CG & !has_GC
            ),
            n_kmers = length(f_rs)
        ) 
    }, 
    spinnerMessage = "loading Tn5 set data"
)

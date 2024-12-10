# build and solve a Hidden Markov model

HMM_init <- function(transProb, emissProbs){ # expects logged probs
    T <- nrow(emissProbs) # i.e., number of observations
    N <- ncol(emissProbs) # i.e., number of hidden states
    transProbs <- matrix(transProb, nrow = N, ncol = N)
    diag(transProbs) <- log(1 - exp(transProb) * (N - 1))
    list(
        emissProbs = emissProbs,
        transProbs = transProbs
    )
}
HMM_viterbi <- function(hmm, observations = TRUE){
    ep <- hmm$emissProbs[observations,] # is a matrix, not a data.table
    tp <- hmm$transProbs
    ep[is.na(ep)] <- -Inf # block unusable paths as having zero probability
    ep[apply(ep, 1, max) == -Inf, ] <- log(1) # mask unusable bins, i.e., those with no usable paths

    # 1. initialization (observation t=1)
    T          <- nrow(ep) # length of the sequence of observations
    N          <- ncol(ep) # number of states
    delta      <- log(matrix(0, nrow = T, ncol = N))
    delta[1, ] <- sapply(1:N, function(i) log(1 / N) + ep[1, i])
    phi        <- matrix(NA, nrow = T, ncol = N)

    # 2. recursion;
    # NB: these 'for' loops are faster than apply methods with array as implemented and given recursion restrictions
    for (t in 2:T){
        pt <- t - 1
        for (j in 1:N){     # j = this hs
            ep_j <- ep[t, j]
            for (i in 1:N){ # i = prev hs
                delta_ <- delta[pt, i] + tp[i, j] + ep_j
                if(delta[t, j] < delta_){
                    delta[t, j] <- delta_
                    phi[pt, j]  <- i
                }
            }
        }
    }
    
    # 3. termination
    prob <- -Inf
    hsi  <- rep(1, T)
    for (j in 1:N){
        if(prob < delta[T, j]){
            prob <- delta[T, j]
            hsi[T] <- j
        }
    }
    
    # 4. reconstruction and return the hidden state indices
    for (t in (T - 1):1) hsi[t] <- phi[t, hsi[t + 1]]
    hsi
}

# likelihood of specific hidden state path using only ep, not tp
HMM_likelihood_obs <- function(HMM, t, n){
    sum(
        HMM$emissProbs[cbind(t,n)], # t and n are vectors of observations and states, respectively
        na.rm = TRUE # remove NA values because some probes in HMM may not be used
    )
}

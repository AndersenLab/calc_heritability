#!/usr/bin/env Rscript H2_script.R
library(boot)
library(lme4)
library(tidyverse)
library(futile.logger)
library(sommer)
library(data.table)

########################
### define functions ###
########################
# Heritability
# data is data frame that contains strain and Value column
# indicies are used by the boot function to sample from the 'data' data.frame
H2.test.boot <- function(data, indices){
    
    d <- data[indices,]
    
    pheno <- as.data.frame(dplyr::select(d, Value))[,1]
    
    strain <- as.factor(d$strain)
    
    reffMod <- lme4::lmer(pheno ~ 1 + (1|strain))
    
    Variances <- as.data.frame(lme4::VarCorr(reffMod, comp = "Variance"))
    
    Vg <- Variances$vcov[1]
    Ve <- Variances$vcov[2]
    H2 <- Vg/(Vg+Ve)
    
    # errors <- sqrt(diag(lme4::VarCorr(reffMod, comp = "Variance")$strain))
    
    return(H2)
}

# data is data frame that contains strain and Value column
H2.test <- function(data){
    
    pheno <- as.data.frame(dplyr::select(data, Value))[,1]
    strain <- as.factor(data$strain)
    
    reffMod <- lme4::lmer(pheno ~ 1 + (1|strain))
    
    Variances <- as.data.frame(lme4::VarCorr(reffMod, comp = "Variance"))
    
    Vg <- Variances$vcov[1]
    Ve <- Variances$vcov[2]
    H2 <- Vg/(Vg+Ve)
    
    # errors <- sqrt(diag(lme4::VarCorr(reffMod, comp = "Variance")$strain))
    
    return(H2)
}

# narrow sense heritability with sommer::mmer (with bootstrap)
narrowh2.boot <- function(data, indices){
    
    d <- data[indices,]
    
    h2_res <- sommer::mmes(Value ~ 1, random = ~sommer::vsm(sommer::ism(strain), Gu = A), data = d)
    
    h2 <- as.numeric(sommer::vpredict(h2_res, h2 ~ (V1) / (V1+V2))[[1]][1])
    
    return(h2)
    
}

# narrow sense heritability with sommer::mmer (no bootstrap)
narrowh2 <- function(df_h){
    
    h2_res <- sommer::mmes(Value ~ 1, random = ~sommer::vsm(sommer::ism(strain), Gu = A), data = df_h)
    
    h2 <- as.numeric(sommer::vpredict(h2_res, h2 ~ (V1) / (V1+V2))[[1]][1])
    
    return(h2)
    
}

# df is data frame that contains strain and Value column
H2.calc <- function(data, boot = TRUE, type = "broad", reps = 500){
    df <- dplyr::select(data, strain, Value)
    
    flog.info("Running bootstrapping")
    if(boot == TRUE){
        # bootstrapping with 10000 replications
        # can reduce value to save time (500 is reasonable most of the time).
        # if you Error in bca.ci(boot.out, conf, index[1L], L = L, t = t.o, t0 = t0.o,  : estimated adjustment 'a' is NA, then you need to increase R value.
        if(type == "broad") {
            results <- boot::boot(data = df, statistic = H2.test.boot, R = reps) 
        } else {
            results <- boot::boot(data = df, statistic = narrowh2.boot, R = reps) 
        }
        
        # get 95% confidence interval
        ci <- boot.ci(results, type="bca")
        
        H2_errors <- data.frame(H2 = ci$t0,
                                ci_l = ci$bca[4],
                                ci_r = ci$bca[5])
        
        return(H2_errors)
        
    } else {
        if(type == "broad") {
            H2 <- data.frame(H2 = H2.test(data = df), ci_l = NA, ci_r = NA)
        } else {
            H2 <- data.frame(H2 = narrowh2(df_h = df), ci_l = NA, ci_r = NA)
        }
        return(H2)
    }
    
}

# RUN
args = commandArgs(trailingOnly=TRUE)

# args <- NULL
# args[1] <- "ExampleTraitData.csv"
# args[2] <- "testResult.tsv"
# args[3] <- "hash.txt"
# args[4] <- "v2"
# args[5] <- "C_elegans_WI_strain_info_20210115"


# Read in data
processed_data <- data.table::fread(args[1])
colnames(processed_data) <- c("strain", "Value")

geno_matrix = data.table::fread(args[2])
# output_fname = args[2]
# heritability_version = args[4]
# hash <- readLines(args[3])

num_rep <- as.numeric(args[3])

# first - check if there is replicate data > 2
test <- processed_data %>%
    dplyr::group_by(strain) %>%
    dplyr::summarize(num = dplyr::n()) %>%
    dplyr::filter(num > 2)


# additive matrix - first filter by strain
A <- sommer::A.mat(t(geno_matrix %>% dplyr::select(dplyr::one_of(processed_data$strain))))

# Run H2 calculation
result_broad <- NULL

# if most of the strains have more than 2 reps, try to bootstrap
if(nrow(test) >= 0.8*length(unique(processed_data$strain))) {
    result_broad <- try(H2.calc(processed_data, boot = T, type = "broad", reps = num_rep) %>%
                            dplyr::mutate(type = "broad-sense"), silent = TRUE)
}


# if result doesn't converge, just give point estimate...
if(!is.data.frame(result_broad)) {
    # if there is no replicates, just return NA
    if(nrow(test) == 0) {
        result_broad <- data.frame(H2 = NA, ci_l = NA, ci_r = NA, type = "broad-sense")
    } else {
        result_broad <- H2.calc(processed_data, boot = F, type = "broad", reps = num_rep) %>%
            dplyr::mutate(type = "broad-sense") 
    }
}

# run h2 calculation
result_narrow <- NULL

if(nrow(test) >= 0.8*length(unique(processed_data$strain))) {
    result_narrow <- try(H2.calc(processed_data, boot = T, type = "narrow", reps = num_rep) %>%
                             dplyr::mutate(type = "narrow-sense"), silent = TRUE)
}

# if result doesn't converge, just give point estimate...
if(!is.data.frame(result_narrow) || nrow(result_narrow) == 0) {
    result_narrow <- H2.calc(processed_data, boot = F, type = "narrow", reps = num_rep) %>%
        dplyr::mutate(type = "narrow-sense")
}

# # combine broad and narrow into one dataframe
result <- rbind(result_broad, result_narrow) %>%
    dplyr::select(type, H2, ci_l, ci_r)

# Write the result
readr::write_tsv(result, "heritability_result.tsv")
# data.table::fwrite(result, output_fname, sep = '\t')
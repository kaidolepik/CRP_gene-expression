library(dplyr)
library(sem)
library(tibble)

# Functions for calculating normal and conditional normal
get_gaussian <- function(obs_values, mean_values, sd_value) {
    gaussian <- -1/2*log(2*pi*sd_value^2) - (obs_values - mean_values)^2/(2*sd_value^2)
    
    return(gaussian)
}

get_a_given_b <- function(a_values, b_values, a_mean, b_mean, a_sd, b_sd, rho) {
    mean_a_given_b <- a_mean + a_sd/b_sd*rho*(b_values - b_mean)
    sd_a_given_b <- (1-rho^2) * a_sd^2
    a_given_b <- get_gaussian(a_values, mean_a_given_b, sd_a_given_b)
    
    return(a_given_b)
}

get_a_given_xb <- function(a_values, b_values, a_means_by_snp, b_mean, a_sd, b_sd, rho) {
    mean_a_given_xb <- a_means_by_snp + a_sd/b_sd*rho*(b_values - b_mean)
    sd_a_given_xb <- (1-rho^2) * a_sd^2
    a_given_xb <- get_gaussian(a_values, mean_a_given_xb, sd_a_given_xb)
    
    return(a_given_xb)
}

# Objective function to optimize according to the specified model
model_func <- function(params, model) {
    mean_gene = params[length(params)-4]
    mean_trait = params[length(params)-3]
    sd_gene = params[length(params)-2]
    sd_trait = params[length(params)-1]
    rho = params[length(params)]
    
    x <- log(global_start_params_by_snp[as.character(snp), "snp_value_prob"])
    
    if (model == "snp->gene->trait") {
        gene_means_by_snp <- params[global_start_params_by_snp[as.character(snp), "index"]]
        
        y_given_x <- get_gaussian(gene, gene_means_by_snp, sd_gene)
        z_given_y <- get_a_given_b(trait, gene, mean_trait, mean_gene, sd_trait, sd_gene, rho)
        likelihoods <- x + y_given_x + z_given_y
    }
    else if (model == "snp->gene<-trait") {
        gene_means_by_snp <- params[global_start_params_by_snp[as.character(snp), "index"]]
        
        z <- get_gaussian(trait, mean_trait, sd_trait)
        y_given_xz <- get_a_given_xb(gene, trait, gene_means_by_snp, mean_trait, sd_gene, sd_trait, rho)
        likelihoods <- x + z + y_given_xz
    }
    else if (model == "snp->trait->gene") {
        trait_means_by_snp <- params[global_start_params_by_snp[as.character(snp), "index"]]
        
        z_given_x <- get_gaussian(trait, trait_means_by_snp, sd_trait)
        y_given_z <- get_a_given_b(gene, trait, mean_gene, mean_trait, sd_gene, sd_trait, rho)
        likelihoods <- x + z_given_x + y_given_z
    }
    else if (model == "snp->gene->trait_and_snp->trait") {
        gene_means_by_snp <- params[global_start_params_by_snp[as.character(snp), "index"]]
        trait_means_by_snp <- params[global_start_params_by_snp[as.character(snp), "index"] + max(global_start_params_by_snp$index)]
        
        y_given_x <- get_gaussian(gene, gene_means_by_snp, sd_gene)
        z_given_xy <- get_a_given_xb(trait, gene, trait_means_by_snp, mean_gene, sd_trait, sd_gene, rho)
        likelihoods <- x + y_given_x + z_given_xy
    }
    
    -sum(likelihoods)
}

# Calculates model AIC
model_AIC <- function(model) {
    start_params <- c()
    if (model == "snp->gene->trait" || model == "snp->gene<-trait" || model == "snp->gene->trait_and_snp->trait")
        start_params <- global_start_params_by_snp$gene_mean_by_snp
    if (model == "snp->trait->gene" || model == "snp->gene->trait_and_snp->trait")
        start_params <- c(start_params, global_start_params_by_snp$trait_mean_by_snp)
    lower <- c(rep(-Inf, length(start_params)), -Inf, -Inf, 1e-8, 1e-8, -1+1e-8)
    start_params <- c(start_params, global_start_params)
    upper <- rep(Inf, length(start_params))
    
    min_log_likelihood <- optim(par = start_params, fn = model_func, model = model, 
                                method = "L-BFGS-B", lower = lower, upper = upper)
    
    2*(min_log_likelihood$value + length(start_params))
}

# Finds AIC differences of the ML approach
ml_AIC_diffs <- function() {
    causal_AIC <- model_AIC("snp->gene->trait")
    colliding_AIC <- model_AIC("snp->gene<-trait")
    reactive_AIC <- model_AIC("snp->trait->gene")
    independent_AIC <- model_AIC("snp->gene->trait_and_snp->trait")
    AICs <- c(causal_AIC, colliding_AIC, reactive_AIC, independent_AIC)
    
    diffs <- outer(AICs, AICs, "-")
    colnames(diffs) <- rownames(diffs) <- c("causal", "colliding", "reactive", "independent")
    
    as.dist(diffs)
}

# AIC differences of the SEM approach
sem_AIC_diffs <- function() {
    causal_equations <- specifyEquations(text = "
        gene = snp_to_gene*snp
        trait = gene_to_trait*gene
        V(snp) = var_snp
        V(gene) = var_gene
        V(trait) = var_trait
    ")
    colliding_equations <- specifyEquations(text = "
        gene = snp_to_gene*snp + trait_to_gene*trait
        V(snp) = var_snp
        V(gene) = var_gene
        V(trait) = var_trait
    ")
    reactive_equations <- specifyEquations(text = "
        trait = snp_to_trait*snp
        gene = trait_to_gene*trait
        V(snp) = var_snp
        V(gene) = var_gene
        V(trait) = var_trait
    ")
    independent_equations <- specifyEquations(text = "
        gene = snp_to_gene*snp
        trait = snp_to_trait*snp + gene_to_trait*gene
        V(snp) = var_snp
        V(gene) = var_gene
        V(trait) = var_trait
    ")
    
    causal_AIC <- summary(sem(causal_equations, data = toy_data))$AIC
    colliding_AIC <- summary(sem(colliding_equations, data = toy_data))$AIC
    reactive_AIC <- summary(sem(reactive_equations, data = toy_data))$AIC
    independent_AIC <- summary(sem(independent_equations, data = toy_data))$AIC
    AICs <- c(causal_AIC, colliding_AIC, reactive_AIC, independent_AIC)
    
    diffs <- outer(AICs, AICs, "-")
    colnames(diffs) <- rownames(diffs) <- c("causal", "colliding", "reactive", "independent")
    
    as.dist(diffs)
}


# Generate toy data
set.seed(12345)
n <- 500
snp <- sample(c(0, 1, 2), n, replace = TRUE, prob = c(0.3^2, 2*0.3*0.7, 0.7^2))
trait <- rnorm(n)
gene <- 0.6*snp + rnorm(n)
trait <- 0.2*gene + rnorm(n)
toy_data <- data.frame(snp, gene, trait)

# Define starting parameter values
global_start_params <- c(mean_gene = mean(gene), 
                         mean_trait = mean(trait),
                         sd_gene = sd(gene), 
                         sd_trait = sd(trait), 
                         rho = cor(gene, trait))
global_start_params_by_snp <- toy_data %>%
    group_by(snp) %>%
    summarize(snp_value_prob = n()/nrow(toy_data),
              gene_mean_by_snp = mean(gene), 
              trait_mean_by_snp = mean(trait)) %>%
    as.data.frame() %>%
    rowid_to_column("index") %>%
    column_to_rownames("snp")

# Calculate AIC differences: rowAIC - columnAIC
ml_AIC_diffs()
sem_AIC_diffs()

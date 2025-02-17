
## Functions to read in data

read_r2_matrix <- function(r2_path){
    r2_df = read.table(r2_path, sep="\t", )
    r2_matrix <- data.matrix(r2_df) 
    return(r2_matrix)
}

get_bim <- function(bim_path){
    bim = read.table(bim_path, sep = '\t')
    names(bim) = c('chr', 'SNP', 'cm', 'bp', 'A1', 'A2')
    return(bim)
    }

get_fam <- function(fam_path){
    fam = read.table(fam_path, sep = ' ')
    names(fam) = c('FID', 'IID', 'fatherID', 'motherID', 'sex', 'phenotype')
    return(fam)
    }

estimate_Ne <- function(mean_r2, S){
    adj1_r2 = mean_r2 * (S/(S-1))**2
    adj2_r2 = adj1_r2 - 0.0018 - 0.907/S - 4.44/(S**2)
    Ne_est = (0.308 + sqrt(.308**2 - 2.08*adj2_r2))/(2*adj2_r2)
    return(Ne_est)
}

get_Ne <- function(base_path){
    # load files
    r2_path  = paste(base_path, '.ld', sep ='')
    bim_path = paste(base_path, '.bim', sep ='')
    fam_path = paste(base_path, '.fam', sep ='')
    pop_mat = read_r2_matrix(r2_path)
    pop_bim = get_bim(bim_path)
    pop_fam = get_fam(fam_path)
    
    # get the sample size of the population from the .fam file
    S = nrow(pop_fam)
    
    # exclude loci on the same chromosome
    for (CH in 1:26){
        my_idx = which(pop_bim$chr==CH)
        pop_mat[my_idx, my_idx] <- NA
    }
    # get just the upper triangle of the square matrix
    r2_vals = pop_mat[upper.tri(x = pop_mat, diag = FALSE)]
    # remove NA values
    r2_vals = r2_vals[!is.na(r2_vals)]
    
    mean_r2 = mean(r2_vals)
    # Non -bias correcteted Ne estimate
    Ne_basic = 1.0/(3*mean_r2 - 3.0/S)
    
    # Bias corrected for low sample size (S<30)
    Ne_est = estimate_Ne(mean_r2=mean_r2, S=S)
    
    #print(c(Ne_est, Ne_basic))
    
    # return the bias-corrected estimate and the
    #  r2 matrix used in the calculation (with within-chromsome r2 values masked)
    return (list(Ne_est = Ne_est, r2_matrix = pop_mat))
}


# Options for the notebook 
options(jupyter.plot_mimetypes = "image/png") 
options(repr.plot.width = 6, repr.plot.height = 6)

read_pca_output <- function(path){
    pca_df = read.table(path)
    names(pca_df) = c('Population', 'Individual', 'PC1', 'PC2', 'PC3')
    # enforce the order of populations
    pca_df$Population <- factor(pca_df$Population , levels =c('Nome_ODD','Nome_EVEN', 'Koppen_ODD', 'Koppen_EVEN', 'Puget_ODD', 'Puget_EVEN'))
    pca_df = pca_df[order(pca_df$Population),]
    return(pca_df)
}

plot_pca_basic <- function(pca_df, title){
    plot(pca_df$PC1 , pca_df$PC2, col = pca_df$Population, pch = 16, main = title)
    legend(x="topleft", legend = levels(pca_df$Population), fill = palette()[1:6])
}

# unused - example how to make a scatteplot in ggplot2
plot_pca_ggplot <- function(pca_df){
    library(ggplot2)
    p = ggplot(data = pca_df, aes(x = PC1, y = PC2))
    p = p + geom_point(aes(colour = Population))
    p = p + theme_classic()
    return(p)
}

pca_initial = read_pca_output('./work/pink_data.initial.eigenvec')
head(pca_initial)

pca_clean = read_pca_output('./work/pink_salmon.clean.eigenvec')
head(pca_clean)

# Make some plots and save them to a file

png("./plots/PCA.pink_salmon.initial.png")
plot_pca_basic(pca_initial, title = 'pink_salmon.initial')
dev.off()

png("./plots/PCA.pink_salmon.clean.png")
plot_pca_basic(pca_clean, title = 'pink_salmon.clean')
dev.off()

print("Done plotting PCAs")



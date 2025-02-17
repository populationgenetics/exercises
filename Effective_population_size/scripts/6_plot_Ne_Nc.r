
# Load in the Ne and Nc estimates

Ne = read.table("./work/Ne_estimates.txt",sep = '\t', header = TRUE)
Nc = read.table("./data/Nc_estimates.txt",sep = '\t', header = TRUE)

Ne

Nc

## Use the merge command to join them

estimates = merge(Ne, Nc)
estimates

estimates$ratio = estimates$Ne_est / estimates$Nc_est
# reorder to match input order
estimates$Pop <- factor(estimates$Pop , levels =c('Nome_ODD','Nome_EVEN', 'Koppen_ODD', 'Koppen_EVEN', 'Puget_ODD', 'Puget_EVEN'))
estimates = estimates[order(estimates$Pop),]
estimates

for_barplot = data.matrix(t(estimates[,c('Ne_est', 'Nc_est')]))
colnames(for_barplot) = estimates$Pop

for_ratio_barplot = data.matrix(t(estimates[,'ratio']))
colnames(for_ratio_barplot) = estimates$Pop


png('./plots/Ne_estimates.png')
par(mar=c(10,4,4,2))
barplot(for_barplot['Ne_est',], col = "white", beside = TRUE, las=2, #axes = FALSE, 
       main = "Ne estimates for each population",ylab = 'Ne')
dev.off()

png('./plots/Ne_and_Nc_estimates.png')
par(mar=c(10,4,4,2))

barplot(for_barplot, col = c("white","black"), beside = TRUE, las=2, axes = FALSE, 
       main = "Ne and Nc estimates for each population",ylab = 'Size')
axis(side = 2, at = c(100, 10000, 500000, 1000000, 1500000))
legend("top",
  c("Ne_est","Nc_est"),
  fill = c("white","black")
)
dev.off()

# same plot with a log y axis
png('./plots/Ne_and_Nc_estimates_log-scaled.png')
par(mar=c(10,4,4,2))
barplot(for_barplot, col = c("white","black"), beside = TRUE, las=2, 
        log = 'y', axes = FALSE, ylim = c(100,1400000), 
       main = "Ne and Nc estimates for each population",ylab = 'Size (log scaled)')
axis(side = 2, at = c(100, 10000, 500000, 1000000, 1500000))
legend("top",
  c("Ne_est","Nc_est"),
  fill = c("white","black")
)
dev.off()

png('./plots/Ne-Nc_ratios.png')
par(mar=c(10,4,4,2))
barplot(for_ratio_barplot, col = "gray", beside = TRUE, las=2, #axes = FALSE, 
       main = "Ne/Nc ratios for each population",ylab = 'Ne/Nc ratio')
dev.off()

source('./scripts/R_functions.r')

Ne_Nome_ODD = get_Ne('./work/Nome_ODD')
Ne_Nome_EVEN = get_Ne('./work/Nome_EVEN')
Ne_Koppen_ODD = get_Ne('./work/Koppen_ODD')
Ne_Koppen_EVEN = get_Ne('./work/Koppen_EVEN')
Ne_Puget_ODD = get_Ne('./work/Puget_ODD')
Ne_Puget_EVEN = get_Ne('./work/Puget_EVEN')

png('./plots/LD_Nome_ODD.png', width=nrow(Ne_Nome_ODD$r2_matrix),height=nrow(Ne_Nome_ODD$r2_matrix))
image(Ne_Nome_ODD$r2_matrix, axes = FALSE, col = rev(heat.colors(256)))
dev.off()

png('./plots/LD_Nome_EVEN.png', width=nrow(Ne_Nome_EVEN$r2_matrix),height=nrow(Ne_Nome_EVEN$r2_matrix))
image(Ne_Nome_EVEN$r2_matrix, axes = FALSE, col = rev(heat.colors(256)))
dev.off()

png('./plots/LD_Koppen_ODD.png', width=nrow(Ne_Koppen_ODD$r2_matrix),height=nrow(Ne_Koppen_ODD$r2_matrix))
image(Ne_Koppen_ODD$r2_matrix, axes = FALSE, col = rev(heat.colors(256)))
dev.off()

png('./plots/LD_Koppen_EVEN.png', width=nrow(Ne_Koppen_EVEN$r2_matrix),height=nrow(Ne_Koppen_EVEN$r2_matrix))
image(Ne_Koppen_EVEN$r2_matrix, axes = FALSE, col = rev(heat.colors(256)))
dev.off()

png('./plots/LD_Puget_ODD.png', width=nrow(Ne_Puget_ODD$r2_matrix),height=nrow(Ne_Puget_ODD$r2_matrix))
image(Ne_Puget_ODD$r2_matrix, axes = FALSE, col = rev(heat.colors(256)))
dev.off()

png('./plots/LD_Puget_EVEN.png', width=nrow(Ne_Puget_EVEN$r2_matrix),height=nrow(Ne_Puget_EVEN$r2_matrix))
image(Ne_Puget_EVEN$r2_matrix, axes = FALSE, col = rev(heat.colors(256)))
dev.off()



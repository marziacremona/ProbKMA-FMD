# Set working directory to the folder of this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

require(fda)
source("../../probKMA-FMD_functions.r")

cer <- function(P,Q){
  # classification error rate (0=perfect agree 1=complete disagreement)
  if(length(P)!=length(Q)){stop('le due partizioni devono avere la stessa lunghezza')}
  cer.comp <- 0
  for(i in 1:(length(P)-1))
    for(j in (i+1):length(P))
      cer.comp <- cer.comp + abs((P[i]==P[j])-(Q[i]==Q[j]))
  cer.comp <- cer.comp/choose(length(P),2)
  return(cer.comp)
}


#############################
### BERKLEY GROWTH CURVES ###
#############################

load("growth_smoothed.RData")




###################
### prob k-mean ###
###################

###################################################################################
### No alignment, no elongation, L2 distance between entire curve derivatives (d1)
###################################################################################

c_min <- 101 # minimum motif length (entire curves)

# Try 10 different initializations, select the best clustering
J_d1 <- rep(NA, 10)
probKMA_results_d1_all <- vector('list', 10)
#for(i in 1:10){
#   probKMA_results_d1_all[[i]] <- probKMA(Y0 = Y0, Y1 = Y1, K = 2, c = c_min, c_max = c_min, diss = 'd1_L2', iter4clean=1000, worker_number = 1)
#   J_d1[i] <- probKMA_results_d1_all[[i]]$J_iter[probKMA_results_d1_all[[i]]$iter]
#}
#save(J_d1, probKMA_results_d1_all, file="./results/probKMA_d1_results.RData")

load("./results/probKMA_d1_results.RData")

J_d1
probKMA_results_d1 <- probKMA_results_d1_all[[which.min(J_d1)]]

probKMA_labels_d1 <- colSums(t(probKMA_results_d1$P > 0.5) * 1:2) # Dicotomize membership based on probability>0.5
probKMA_P_d1 <- round(probKMA_results_d1$P[, 2], 2)
col_P_d1 <- rgb(colorRamp(c("red", "blue"))( probKMA_P_d1 ), maxColorValue = 255)

# Comparison with sex classification
table(real_labels, probKMA_labels_d1)
cer(real_labels, probKMA_labels_d1)

# Look at probabilities of misclassified boys and girls
boys_mis <- which((real_labels == "boys") & (probKMA_labels_d1 == 2))
girls_mis <- which((real_labels == "girls") & (probKMA_labels_d1 == 1))

pdf('./results/probKMA_d1_misclassified.pdf', width = 7, height = 5.5)
boxplot(list(Boys = probKMA_P_d1[real_labels == "boys"], Girls = probKMA_P_d1[real_labels == "girls"]), col = c("cyan3", "magenta"), ylab = "Probabilistic membership",main = "probKMA d1", cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
points(rep(1, length(boys_mis)), probKMA_P_d1[boys_mis], pch = 16, col = col_P_d1[boys_mis])
points(rep(2, length(girls_mis)), probKMA_P_d1[girls_mis], pch = 16, col = col_P_d1[girls_mis])
dev.off()

# Plot
probKMA_align_d1 <- rep(NA, length(probKMA_labels_d1))
probKMA_align_d1[probKMA_labels_d1 == 1] <- probKMA_results_d1$S_clean[probKMA_labels_d1 == 1, 1]
probKMA_align_d1[probKMA_labels_d1 == 2] <- probKMA_results_d1$S_clean[probKMA_labels_d1 == 2, 2]
probKMA_lengths_d1 <- unlist(lapply(probKMA_results_d1$V0, nrow))
Y0_mat <- Reduce(cbind, Y0)
Y1_mat <- Reduce(cbind, Y1)

pdf('./results/probKMA_d1_results.pdf', width = 7, height = 5.5)
matplot(x, Y0_mat, type = 'l', col = 'darkgray', xlab = 'Age (years)', ylab = 'Height (cm)', main = 'probKMA d1', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
for(i in seq_along(probKMA_labels_d1)){
  index <- probKMA_align_d1[i] + seq_len(probKMA_lengths_d1[probKMA_labels_d1[i]]) - 1
  index <- index[(index > 0) & (index <= length(x))]
  lines(x[index], Y0[[i]][index, ], col = col_P_d1[i], lty = i, lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
}
legend('bottomright', legend = c("Cluster 1", "Cluster 2"), col = c("blue", "red"), lty = 1, lwd = 1.5, cex = 1.5, bty = "n")
matplot(x, Y1_mat, type = 'l', col = 'darkgray', xlab = 'Age (years)', ylab = 'Growth (cm/year)', main = 'probKMA d1', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
for(i in seq_along(probKMA_labels_d1)){
  index <- probKMA_align_d1[i] + seq_len(probKMA_lengths_d1[probKMA_labels_d1[i]]) - 1
  index <- index[(index > 0) & (index <= length(x))]
  lines(x[index], Y1[[i]][index, ], col = col_P_d1[i], lty = i, lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
}
legend('topright', legend = c("Cluster 1", "Cluster 2"), col = c("blue", "red"), lty = 1, lwd = 1.5, cex = 1.5, bty = "n")

# Plot of "uncertain" memberships vs cluster 1 vs cluster 2
i_cluster1 <- which((probKMA_P_d1 >= 0.6))
i_cluster2 <- which((probKMA_P_d1 <= 0.4))
i_uncertain <- which((probKMA_P_d1 > 0.4) & (probKMA_P_d1 < 0.6))

matplot(x, Y0_mat, type = 'l', col = 'darkgray', xlab = 'Age (years)', ylab = 'Height (cm)', main = 'probKMA d1 - Cluster 1', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
for(i in i_cluster1){
  index <- probKMA_align_d1[i] + seq_len(probKMA_lengths_d1[probKMA_labels_d1[i]]) - 1
  index <- index[(index > 0) & (index <= length(x))]
  lines(x[index], Y0[[i]][index, ], col = col_P_d1[i], lty = i, lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
}
matplot(x, Y0_mat, type = 'l', col = 'darkgray', xlab = 'Age (years)', ylab = 'Height (cm)', main = 'probKMA d1 - Cluster 2', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
for(i in i_cluster2){
  index <- probKMA_align_d1[i] + seq_len(probKMA_lengths_d1[probKMA_labels_d1[i]]) - 1
  index <- index[(index > 0) & (index <= length(x))]
  lines(x[index], Y0[[i]][index, ], col = col_P_d1[i], lty = i, lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
}
matplot(x, Y0_mat, type = 'l', col = 'darkgray', xlab = 'Age (years)', ylab = 'Height (cm)', main = 'probKMA d1 - Uncertain', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
for(i in i_uncertain){
  index <- probKMA_align_d1[i] + seq_len(probKMA_lengths_d1[probKMA_labels_d1[i]]) - 1
  index <- index[(index > 0) & (index <= length(x))]
  lines(x[index], Y0[[i]][index, ], col = col_P_d1[i], lty = i, lwd = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
}
matplot(x, Y1_mat, type = 'l', col = 'darkgray', xlab = 'Age (years)', ylab = 'Growth (cm/year)', main = 'probKMA d1 - Cluster 1', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
for(i in i_cluster1){
  index <- probKMA_align_d1[i] + seq_len(probKMA_lengths_d1[probKMA_labels_d1[i]]) - 1
  index <- index[(index > 0) & (index <= length(x))]
  lines(x[index], Y1[[i]][index, ], col = col_P_d1[i], lty = i, lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
}
matplot(x, Y1_mat, type = 'l', col = 'darkgray', xlab = 'Age (years)', ylab = 'Growth (cm/year)', main = 'probKMA d1 - Cluster 2', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
for(i in i_cluster2){
  index <- probKMA_align_d1[i] + seq_len(probKMA_lengths_d1[probKMA_labels_d1[i]]) - 1
  index <- index[(index > 0) & (index <= length(x))]
  lines(x[index], Y1[[i]][index, ], col = col_P_d1[i], lty = i, lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
}
matplot(x, Y1_mat, type = 'l', col = 'darkgray', xlab = 'Age (years)', ylab = 'Growth (cm/year)', main = 'probKMA d1 - Uncertain', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
for(i in i_uncertain){
  index <- probKMA_align_d1[i] + seq_len(probKMA_lengths_d1[probKMA_labels_d1[i]]) - 1
  index <- index[(index > 0) & (index <= length(x))]
  lines(x[index], Y1[[i]][index, ], col = col_P_d1[i], lty = i, lwd = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
}
dev.off()
#####



########################################
### prob k-mean with local alignment ###
########################################

########################################################################
### Shift alignment, L2 distance between derivative curves (d1), c=51
########################################################################

c_min <- 51 # Allow a shift of around 8.5 years

# Try 10 different initializations, select the best clustering
J_d1_shift_c51 <- rep(NA, 10)
probKMA_results_d1_shift_c51_all <- vector('list', 10)
#for(i in 1:10){
#  probKMA_results_d1_shift_c51_all[[i]] <- probKMA(Y0 = Y0, Y1 = Y1, K = 2, c = c_min, c_max = c_min, diss = 'd1_L2', iter4clean=1000, worker_number = 1)
#  J_d1_shift_c51[i] <- probKMA_results_d1_shift_c51_all[[i]]$J_iter[probKMA_results_d1_shift_c51_all[[i]]$iter]
#}
#save(J_d1_shift_c51, probKMA_results_d1_shift_c51_all, file="./results/probKMA_d1_shift_c51_results.RData")

load("./results/probKMA_d1_shift_c51_results.RData")

J_d1_shift_c51

####### Use probKMA dichotomization of membership based on distance <= median(distances)
probKMA_results_d1_shift_c51 <- probKMA_results_d1_shift_c51_all[[which.min(J_d1_shift_c51)]]
probKMA_labels_d1_shift_c51 <- colSums(t(probKMA_results_d1_shift_c51$P_clean) * 1:2) 

pdf('./results/probKMA_d1_shift_c51_distances.pdf', width = 7, height = 5.5)
hist(probKMA_results_d1_shift_c51$D, breaks = 30, xlab = "Distances portions to cluster centers", main = "probKMA d1 c=8.5 years", cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
abline(v = median(probKMA_results_d1_shift_c51$D), col = "red", lty = 2, lwd = 2)
dev.off()

pdf('./results/probKMA_d1_shift_c51_silhouette.pdf', width = 7, height = 8.5)
par(cex.main = 2, cex.axis = 1.5, cex.lab = 1.5)
probKMA_silhouette(probKMA_results_d1_shift_c51)
dev.off()
probKMA_P_d1_shift_c51 <- round(probKMA_results_d1_shift_c51$P[, 1], 2)
col_P_d1_shift_c51 <- rgb(colorRamp(c("red", "blue"))( probKMA_P_d1_shift_c51 ), maxColorValue = 255)

# Comparison with sex classification
table(real_labels, probKMA_labels_d1_shift_c51)

# Plot
probKMA_align_d1_shift_c51 <- probKMA_results_d1_shift_c51$S_clean
probKMA_lengths_d1_shift_c51 <- unlist(lapply(probKMA_results_d1_shift_c51$V0, nrow))
Y0_mat <- Reduce(cbind, Y0)
Y1_mat <- Reduce(cbind, Y1)

pdf('./results/probKMA_d1_shift_c51_results.pdf', width = 7, height = 5.5)
matplot(x, Y0_mat, type = 'l', col = 'darkgray', xlab = 'Age (years)', ylab = 'Height (cm)', main = 'probKMA d1 c=8.5 years', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
for(i in seq_along(probKMA_labels_d1_shift_c51)){
  if(probKMA_labels_d1_shift_c51[i] %in% c(1, 2)){
    index <- probKMA_align_d1_shift_c51[i, probKMA_labels_d1_shift_c51[i]] + 
      seq_len(probKMA_lengths_d1_shift_c51[probKMA_labels_d1_shift_c51[i]]) - 1
    index <- index[(index > 0) & (index <= length(x))]
    lines(x[index], Y0[[i]][index, ], col = col_P_d1_shift_c51[i], lty = i, lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
  }else if(probKMA_labels_d1_shift_c51[i] == 3){
    index <- probKMA_align_d1_shift_c51[i, 1] + seq_len(probKMA_lengths_d1_shift_c51[1]) - 1
    index <- index[(index > 0) & (index <= length(x))]
    lines(x[index], Y0[[i]][index, ], col = col_P_d1_shift_c51[i], lty = i, lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
    index <- probKMA_align_d1_shift_c51[i, 2] + seq_len(probKMA_lengths_d1_shift_c51[2]) - 1
    index <- index[(index > 0) & (index <= length(x))]
    lines(x[index], Y0[[i]][index, ], col = col_P_d1_shift_c51[i], lty = i, lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
  }
}
legend('bottomright', legend = c("Cluster 1", "Cluster 2"), col = c("blue", "red"), lty = 1, lwd = 1.5, cex = 1.5, bty = "n")
matplot(x, Y1_mat, type = 'l', col = 'darkgray', xlab = 'Age (years)', ylab = 'Growth (cm/year)', main = 'probKMA d1 c=8.5 years', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
for(i in seq_along(probKMA_labels_d1_shift_c51)){
  if(probKMA_labels_d1_shift_c51[i] %in% c(1, 2)){
    index <- probKMA_align_d1_shift_c51[i, probKMA_labels_d1_shift_c51[i]] + 
      seq_len(probKMA_lengths_d1_shift_c51[probKMA_labels_d1_shift_c51[i]]) - 1
    index <- index[(index > 0) & (index <= length(x))]
    lines(x[index], Y1[[i]][index, ], col = col_P_d1_shift_c51[i], lty = i, lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
  }else if(probKMA_labels_d1_shift_c51[i] == 3){
    index <- probKMA_align_d1_shift_c51[i, 1] + seq_len(probKMA_lengths_d1_shift_c51[1]) - 1
    index <- index[(index > 0) & (index <= length(x))]
    lines(x[index], Y1[[i]][index, ], col = col_P_d1_shift_c51[i], lty = i, lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
    index <- probKMA_align_d1_shift_c51[i, 2] + seq_len(probKMA_lengths_d1_shift_c51[2]) - 1
    index <- index[(index > 0) & (index <= length(x))]
    lines(x[index], Y1[[i]][index, ], col = col_P_d1_shift_c51[i], lty = i, lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
  }
}
legend('topright', legend = c("Cluster 1", "Cluster 2"), col = c("blue", "red"), lty = 1, lwd = 1.5, cex = 1.5, bty = "n")

# Plot of cluster 1 vs cluster 2 aligned
i_cluster1 <- which(probKMA_labels_d1_shift_c51 %in% c(1, 3))
i_cluster2 <- which(probKMA_labels_d1_shift_c51 %in% c(2, 3))

Y0_cluster1 <- c()
for(i in i_cluster1)
  Y0_cluster1 <- cbind(Y0_cluster1, Y0_mat[ , i][probKMA_align_d1_shift_c51[i, 1] + seq_len(probKMA_lengths_d1_shift_c51[1]) - 1])
Y0_cluster1_mean <- rowMeans(Y0_cluster1, na.rm = TRUE)
Y0_cluster2 <- c()
for(i in i_cluster2)
  Y0_cluster2 <- cbind(Y0_cluster2, Y0_mat[ , i][probKMA_align_d1_shift_c51[i, 2] + seq_len(probKMA_lengths_d1_shift_c51[2]) - 1])
Y0_cluster2_mean <- rowMeans(Y0_cluster2, na.rm = TRUE)
matplot(x[median(probKMA_align_d1_shift_c51[i_cluster1, 1]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[1])], 
        Y0_cluster1, type = 'l', col = "blue", xlim = range(x), ylim = range(cbind(Y0_cluster1, Y0_cluster2)),
        xlab = 'Age (years)', ylab = 'Height (cm)', main = 'probKMA d1 c=8.5 years - Cluster 1', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
lines(x[median(probKMA_align_d1_shift_c51[i_cluster1, 1]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[1])], 
      Y0_cluster1_mean, col = "black", lwd = 3)
matplot(x[median(probKMA_align_d1_shift_c51[i_cluster2, 2]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[2])], 
        Y0_cluster2, type = 'l', col = "red", xlim = range(x), ylim = range(cbind(Y0_cluster1, Y0_cluster2)),
        xlab = 'Age (years)', ylab = 'Height (cm)', main = 'probKMA d1 c=8.5 years - Cluster 2', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
lines(x[median(probKMA_align_d1_shift_c51[i_cluster2, 2]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[2])], 
      Y0_cluster2_mean, col = "black", lwd = 3)
# together
matplot(x[median(probKMA_align_d1_shift_c51[i_cluster1, 1]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[1])], 
        Y0_cluster1, type = 'l', col = "blue", xlim = range(x), ylim = range(cbind(Y0_cluster1, Y0_cluster2)),
        xlab = 'Age (years)', ylab = 'Height (cm)', main = 'probKMA d1 c=8.5 years', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
lines(x[median(probKMA_align_d1_shift_c51[i_cluster1, 1]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[1])], 
      Y0_cluster1_mean, col = "black", lwd = 3)
matplot(x[median(probKMA_align_d1_shift_c51[i_cluster2, 2]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[2])], 
        Y0_cluster2, type = 'l', col = "red", xlim = range(x), ylim = range(cbind(Y0_cluster1, Y0_cluster2)),
        xlab = 'Age (years)', ylab = 'Height (cm)', main = 'probKMA d1 c=8.5 years', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, add = TRUE)
lines(x[median(probKMA_align_d1_shift_c51[i_cluster2, 2]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[2])], 
      Y0_cluster2_mean, col = "black", lwd = 3)
legend('bottomright', legend = c("Cluster 1", "Cluster 2"), col = c("blue", "red"), lty = 1, lwd = 2, cex = 1.5, bty = "n")

Y1_cluster1 <- c()
for(i in i_cluster1)
  Y1_cluster1 <- cbind(Y1_cluster1, Y1_mat[ , i][probKMA_align_d1_shift_c51[i, 1] + seq_len(probKMA_lengths_d1_shift_c51[1]) - 1])
Y1_cluster1_mean <- rowMeans(Y1_cluster1, na.rm = TRUE)
Y1_cluster2 <- c()
for(i in i_cluster2)
  Y1_cluster2 <- cbind(Y1_cluster2, Y1_mat[ , i][probKMA_align_d1_shift_c51[i, 2] + seq_len(probKMA_lengths_d1_shift_c51[2]) - 1])
Y1_cluster2_mean <- rowMeans(Y1_cluster2, na.rm = TRUE)
matplot(x[median(probKMA_align_d1_shift_c51[i_cluster1, 1]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[1])], 
        Y1_cluster1, type = 'l', col = "blue", xlim = range(x), ylim = range(cbind(Y1_cluster1, Y1_cluster2)),
        xlab = 'Age (years)', ylab = 'Growth (cm/year)', main = 'probKMA d1 c=8.5 years - Cluster 1', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
lines(x[median(probKMA_align_d1_shift_c51[i_cluster1, 1]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[1])], 
      Y1_cluster1_mean, col = "black", lwd = 3)
matplot(x[median(probKMA_align_d1_shift_c51[i_cluster2, 2]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[2])], 
        Y1_cluster2, type = 'l', col = "red", xlim = range(x), ylim = range(cbind(Y1_cluster1, Y1_cluster2)), 
        xlab = 'Age (years)', ylab = 'Growth (cm/year)', main = 'probKMA d1 c=8.5 years - Cluster 2', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
lines(x[median(probKMA_align_d1_shift_c51[i_cluster2, 2]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[2])], 
      Y1_cluster2_mean, col = "black", lwd = 3)
#together
matplot(x[median(probKMA_align_d1_shift_c51[i_cluster1, 1]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[1])], 
        Y1_cluster1, type = 'l', col = "blue", xlim = range(x), ylim = range(cbind(Y1_cluster1, Y1_cluster2)),
        xlab = 'Age (years)', ylab = 'Growth (cm/year)', main = 'probKMA d1 c=8.5 years', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
lines(x[median(probKMA_align_d1_shift_c51[i_cluster1, 1]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[1])], 
      Y1_cluster1_mean, col = "black", lwd = 3)
matplot(x[median(probKMA_align_d1_shift_c51[i_cluster2, 2]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[2])], 
        Y1_cluster2, type = 'l', col = "red", xlim = range(x), ylim = range(cbind(Y1_cluster1, Y1_cluster2)), 
        xlab = 'Age (years)', ylab = 'Growth (cm/year)', main = 'probKMA d1 c=8.5 years', lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, add = TRUE)
lines(x[median(probKMA_align_d1_shift_c51[i_cluster2, 2]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[2])], 
      Y1_cluster2_mean, col = "black", lwd = 3)
legend('topright', legend = c("Cluster 1", "Cluster 2"), col = c("blue", "red"), lty = 1, lwd = 2, cex = 1.5, bty = "n")

matplot(cbind(x[median(probKMA_align_d1_shift_c51[i_cluster1, 1]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[1])], 
              x[median(probKMA_align_d1_shift_c51[i_cluster2, 2]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[2])]), 
        cbind(Y0_cluster1_mean, Y0_cluster2_mean), col = c("blue", "red"), type='l', lwd = 2, xlim = range(x), ylim = range(cbind(Y0_cluster1, Y0_cluster2)),
        xlab = 'Age (years)', ylab = 'Height (cm)', main = 'probKMA d1 c=8.5 years - Averages', cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
legend('bottomright', legend = c("Cluster 1", "Cluster 2"), col = c("blue", "red"), lty = 1, lwd = 2, cex = 1.5, bty = "n")
matplot(cbind(x[median(probKMA_align_d1_shift_c51[i_cluster1, 1]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[1])], 
              x[median(probKMA_align_d1_shift_c51[i_cluster2, 2]) - 1 + seq_len(probKMA_lengths_d1_shift_c51[2])]), 
        cbind(Y1_cluster1_mean, Y1_cluster2_mean), col = c("blue", "red"), type='l', lwd = 2, xlim = range(x), ylim = range(cbind(Y1_cluster1, Y1_cluster2)),
        xlab = 'Age (years)', ylab = 'Growth (cm/year)', main = 'probKMA d1 c=8.5 years - Averages', cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
legend('topright', legend = c("Cluster 1", "Cluster 2"), col = c("blue", "red"), lty = 1, lwd = 2, cex = 1.5, bty = "n")
dev.off()
######


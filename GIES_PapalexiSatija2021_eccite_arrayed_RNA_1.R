library("pcalg")
## Load predefined data
working_residuals <- read.csv("D:/SummerIntern/Datasets/PapalexiSatija2021_eccite_arrayed_RNA_working_residuals.csv", header = TRUE)
residuals_matrix <- as.matrix(sapply(working_residuals, as.numeric))
residuals_matrix <- residuals_matrix[,-1]
set.seed(56)
p <- ncol(residuals_matrix)
n <- nrow(residuals_matrix)
Perturbation_Targets <- read.csv("D:/SummerIntern/Datasets/PapalexiSatija2021_eccite_arrayed_RNA_perturbation_information.csv", header = TRUE)
gGtrue <- randomDAG(p, prob = 0.3)
targets <- list(integer(0),1)
target_genes <- c("IRF1")
target.index <- c(rep(1,n))
for (i in 1:n){
  for (j in 1: length(target_genes)){
    if (Perturbation_Targets[i,2] == target_genes[j]){
      target.index[i] = j+1
    }
  }
}

scPerturb_data <- list(x = residuals_matrix, 
                       targets = targets, 
                       target.index = target.index, 
                       g = gGtrue)

## Define the score (BIC)
score <- new("GaussL0penIntScore", scPerturb_data$x, scPerturb_data$targets, scPerturb_data$target.index)
## Estimate the essential graph
gies.fit <- gies(score)
## Plot the estimated essential graph and the true DAG
plot(gies.fit$essgraph, main = "Estimated CPDAG for leukaemia genes")

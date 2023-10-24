library("pcalg")
library(igraph)
library(readxl)
library("caret")
library(mltools)
library(SID)
start_time = Sys.time()
## Load predefined data
working_residuals <- read.csv("D:/SummerIntern/Datasets/MAGIC_PapalexiSatija2021_eccite_RNA_working_residuals2_exp1.csv", header = TRUE)
residuals_matrix <- as.matrix(sapply(working_residuals, as.numeric))
residuals_matrix <- residuals_matrix[,-1]
set.seed(5)
p <- ncol(residuals_matrix)
n <- nrow(residuals_matrix)
Perturbation_Targets <- read.csv("D:/SummerIntern/Datasets/PapalexiSatija2021_eccite_RNA_perturbation_information_exp1.csv", header = TRUE)
gGtrue <- randomDAG(p, prob = 0.3)
targets <- list(integer(0),7,8,9,10)
target_genes <- c('STAT1', 'MYC', 'STAT3', 'JAK2')
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
est_igraph <- as(gies.fit$essgraph, "graphNEL")
est_ig <- igraph.from.graphNEL(est_igraph)
plot(est_ig,, main = "Est-pathway 1")
sparse_matrix <- as_adjacency_matrix(est_ig)
est_adj <- as.matrix(sparse_matrix)
end_time = Sys.time()
running_time = end_time - start_time
print(running_time)
adj_mat <- read_excel("D:/SummerIntern/Datasets/pathway information-JAK STAT.xlsx")
adj_mat <- as.matrix(sapply(adj_mat, as.numeric))
adj_mat <- adj_mat[,-1]
rownames(adj_mat) <- colnames(adj_mat)
graph_true <- graph_from_adjacency_matrix(adj_mat)
plot(graph_true, main = "True-pathway 1")
index <- as.vector(colnames(adj_mat))
x1 <- est_adj[index,index]
est_adj <- x1
sid <- structIntervDist(est_adj, adj_mat, output = FALSE, spars = FALSE)
shd <- hammingDist(est_adj, adj_mat)
mcc <- mcc(as.vector(est_adj),as.vector(adj_mat))
result <- confusionMatrix(factor(est_adj), factor(adj_mat), positive = "1")
TN <- result$table[1,1]
FN <- result$table[1,2]
FP <- result$table[2,1]
TP <- result$table[2,2]
FDR <- FP/(TP+FP)
TPR <- TP/(TP+FN)
print(result)
cat("FDR =",FDR)
cat("TPR = ",TPR)
print(mcc)
# 
kf=(commandArgs(TRUE))
print(kf)

library(R.matlab)
library(igraph)

### update according to file  location ###
##########################################

#filenameR <- dir("../../Results/","*.RData")

# the subset for which we have the stoichiometry
#filename_DONE <- dir("../../Results/stoichiometric_coupling/","*.mat")

#filename_DONE <- gsub(".mat",".RData",filename_DONE)
#filename_DONE <- gsub("stoichiometric_coupling","concordant",filename_DONE)

#filenameR = intersect(filenameR,filename_DONE)

###########################################

#print(length(filenameR))
#print(filenameR[k])

###results for the Excel tables -- geneated species by species (given by the species order in the filename vector)

load(file = paste("../../Results/", kf, sep = ""))
OverviewTable = rep(0, 22)

Y <- as.matrix(s$Results.balanced[,,1]$MODEL.r[[1]][[1]][,,1]$Y)
met_lkm <- rep(0, length(mets))
for (i in 1:nrow(Y))
{
	p <- which(Y[i,] > 0)
	if (length(which(p %in% final_lres[[1]])) == length(p)) # if all complex with metabolite i are in one kinetic module
	{
		met_lkm[i] <- 1
	}
}
length(which(met_lkm == 1))

pEC <- s$Results.balanced[,,1]$MODEL.r[[1]][[1]][,,1]$rxnECNumber
ECnum <- rep("", length(pEC))
for (k in 1:length(pEC))
{
	p <- unlist(pEC[[k]])
	if (length(p) > 0)
	{
		pp <- unlist(strsplit(p, ""))
		ECnum[k] <- paste(pp[1], pp[2], pp[3], sep = "")
	}
}

diff_one_met <- c()
diff_two_met <- c()
for (j in 1:length(final_lres)) # for kinetic modules
{
if (length(final_lres[[j]]) > 1) # if there is more than one complex assigned
{
Y_sub <- Y[,final_lres[[j]]]
for (k in 1:(ncol(Y_sub) - 1)) # for two complexes i and j  in the kinetic module
{
	for (m in (k+1):ncol(Y_sub))
	{
		pp <- Y_sub[,k] - Y_sub[,m] 
		q <- which(pp != 0) # if they differ in exactly one metabolite
		if (length(q) == 1)
		{                                   # metabolite they differ, complex i, complex j, kinetic module
			diff_one_met <- rbind(diff_one_met, c(q[1], final_lres[[j]][k], final_lres[[j]][m], j))
		}
		qq1 <- length(which(pp != 0)) == 2
		qq2 <- length(setdiff(which(Y_sub[,k] != 0), which(pp != 0))) == 1
		qq3 <- length(setdiff(which(Y_sub[,m] != 0), which(pp != 0))) == 1
		if (qq1 && qq2 && qq3)
		{
			diff_two_met <- rbind(diff_two_met, c(q[1], q[2], final_lres[[j]][k], final_lres[[j]][m], j))
		}
	}
}
}
}

metNames <- unlist(s$Results.balanced[,,1]$MODEL.r[[1]][[1]][,,1]$metNames)

# #stoichiometric coupling
print('stoichiometric coupling')
q <- readMat(paste("../../Results/stoichiometric_coupling/", gsub("RData","mat",gsub("concordant","stoichiometric_coupling",kf)), sep = ""))
# #stoichiometric coupling of reactions
A_mat <- abs(q$coupling.matrix)
A_mat[which(A_mat == 2)] <- 0
A_mat[which(A_mat == 3)] <- 0
A_mat[which(A_mat == 4)] <- 0
gsc <- graph.adjacency(A_mat, mode = "undirected", diag = F)
pg <- components(gsc)
length(pg$csize) # number of components
max(pg$csize) # maximum number of complexes in a component
mean(pg$csize)
sd(pg$csize)
length(which(pg$csize == 1)) # number of clusters with a single complex
# #stoichiometric coupling of complexes (project reactions to complexes)
A_complex <- matrix(0, vcount(g), vcount(g))
diff_two_met_full_coup <- c()
for (i in 1:(nrow(A_mat)-1))
{
	for (j in (i+1):nrow(A_mat))
	{
		if (A_mat[i,j] == 1)
		{
			comp_i <- which(adjm[,i] == -1)
			comp_j <- which(adjm[,j] == -1)
			A_complex[comp_i, comp_j] <- A_complex[comp_j, comp_i] <- 1
		}
	}
}
g_scomplex <- graph.adjacency(A_complex, mode = "undirected", diag = F)
pgs <- components(g_scomplex)
OverviewTable[9] = length(pgs$csize)
OverviewTable[11] = max(pgs$csize)
OverviewTable[12] = mean(pgs$csize)
OverviewTable[13] = sd(pgs$csize)
OverviewTable[10] = length(which(pgs$csize == 1))

# 
print('singles / doubles from full coupling')
diff_one_met_full_coup <- c()
diff_two_met_full_coup <- c()
for (j in 1:length(pgs$no)) # for coupling modules
{
  if (length(pgs$csize[j]) > 1) # if there is more than one complex assigned
  {
    Y_sub <- Y[,which(pgs$membership==j)]
    for (k in 1:(ncol(Y_sub) - 1)) # for two complexes i and j  in the kinetic module
    {
      for (m in (k+1):ncol(Y_sub))
      {
        pp <- Y_sub[,k] - Y_sub[,m] 
        q <- which(pp != 0) # if they differ in exactly one metabolite
        if (length(q) == 1)
        {                                   # metabolite they differ, complex i, complex j, kinetic module
          diff_one_met_full_coup <- rbind(diff_one_met_full_coup, c(q[1], which(pgs$membership==j)[k], which(pgs$membership==j)[m], j))
        }
        qq1 <- length(which(pp != 0)) == 2
        qq2 <- length(setdiff(which(Y_sub[,k] != 0), which(pp != 0))) == 1
        qq3 <- length(setdiff(which(Y_sub[,m] != 0), which(pp != 0))) == 1
        if (qq1 && qq2 && qq3)
        {
          diff_two_met_full_coup <- rbind(diff_two_met_full_coup, c(q[1], q[2], which(pgs$membership==j)[k], which(pgs$membership==j)[m], j))
        }
      }
    }
  }
}


#robustness metab / pairs
print('get metab/ pairs ...')
OverviewTable[14] = length(unique(diff_one_met[,1]))
if (length(unique(diff_two_met[,1:2])) > 0)
{

if (is.null(nrow(unique(diff_two_met[,1:2])))) {
OverviewTable[15] = 1
} else {
OverviewTable[15] = nrow(unique(diff_two_met[,1:2]))
}

#clusters of coupled metabolites
Amet <- matrix(0, nrow(Y), nrow(Y))
for (i in 1:nrow(diff_two_met))
{
	Amet[diff_two_met[i, 1], diff_two_met[i, 2]] <- Amet[diff_two_met[i, 2], diff_two_met[i, 1]] <- 1
}
gmet <- graph.adjacency(Amet, mode = "undirected", diag = F)
pgm <- components(gmet)
OverviewTable[16] = length(pgm$csize)
OverviewTable[18] = max(pgm$csize)
OverviewTable[19] = mean(pgm$csize)
OverviewTable[20] = sd(pgm$csize)
OverviewTable[17] = length(which(pgm$csize == 1))
}
q1 <- as.matrix(metNames[unique(diff_one_met[,1])])
qq <- unique(diff_two_met[,1:2])
if (OverviewTable[15] == 1){
q2 = metNames[qq]
} else {
q2 <- cbind(metNames[qq[,1]], metNames[qq[,2]])
}
write.table(q1, file = paste("../../Results/MetSingle/MetSingle_", gsub(".RData","",kf), ".csv", sep = ""))
write.table(q2, file = paste("../../Results/MetDouble/MetDouble_", gsub(".RData","",kf), ".csv", sep = ""))


#kinetic coupling
print('get kinetic coupling ...')
OverviewTable[6] = max(unlist(lapply(final_lres, length)))
OverviewTable[7] = mean(unlist(lapply(final_lres, length)))
OverviewTable[8] = sd(unlist(lapply(final_lres, length)))
OverviewTable[5] = single_all
OverviewTable[4] = length(final_lres)

#basic model properties
OverviewTable[1] = vcount(g)
OverviewTable[2] = ecount(g)
OverviewTable[3] = dim(Y)[1]

# #how many reactions are with substrate in giant module
 gm <- which.max(unlist(lapply(final_lres, length)))
 final_lres[[gm]] = final_lres[[gm]][which(is.na(final_lres[[gm]]) == FALSE)]
 A <- s$Results.balanced[,,1]$MODEL.r[[1]][[1]][,,1]$A
 reactions_in_giant <- c()
 if (length(final_lres[[gm]]) > 0){
 for (i in 1:length(final_lres[[gm]]))
 {
 	reactions_in_giant <- c(reactions_in_giant, which(A[final_lres[[gm]][i],] == -1))
 }
 OverviewTable[23] = length(reactions_in_giant)
}

write.table(OverviewTable, file = paste("../../Results/Overview/Overview_", gsub(".RData","",kf), ".csv", sep = ""))


# #how many reactions are with substrate in giant module
gm <- which.max(unlist(lapply(final_lres, length)))
final_lres[[gm]] = final_lres[[gm]][which(is.na(final_lres[[gm]]) == FALSE)]
A <- s$Results.balanced[,,1]$MODEL.r[[1]][[1]][,,1]$A
reactions_in_giant <- c()
for (i in 1:length(final_lres[[gm]]))
{
	reactions_in_giant <- c(reactions_in_giant, which(A[final_lres[[gm]][i],] == -1))
}
length(reactions_in_giant)
write.table(reactions_in_giant, file = paste("../../Results/Reactions_Giant/Reactions_Giant_", gsub(".RData","",kf), ".csv", sep = ""))
 
print('Done')
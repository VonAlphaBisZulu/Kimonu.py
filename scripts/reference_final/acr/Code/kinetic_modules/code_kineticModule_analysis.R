##program to calculate kinetic modules
k=(commandArgs(TRUE))
#k=as.numeric(k)
print(k)
print('------------')

library(R.matlab)
library(igraph)
setwd('../../Results/concordant/')

#filename <- dir(".","*.mat")

print(k)
s <- readMat(k)
adjm <- s$Results.balanced[,,1]$MODEL.r[[1]][[1]][,,1]$A
g <- graph.empty(n = nrow(adjm), directed = TRUE)
for (j in 1:ncol(adjm))
{
	g <- add_edges(g, c(which(adjm[,j] == -1), which(adjm[,j] == 1)))
}
g <- delete_edges(g, which(s$Elementary.Kinetic.Param == 0))
clust <- components(g, "strong")

terminal_c <- rep(1, clust$no) #1 means terminal
for (i in 1:clust$no){
	vc <- which(clust$membership == i)
	sn <- unlist(ego(g, 1, vc, mode = "out"))
	sdsn <- setdiff(sn, vc)
	if (length(sdsn) > 0)
		terminal_c[i] <- 0 		
}

terminal <- rep(1, length(clust$membership))
for (i in 1:length(terminal))
{
	terminal[i] <- terminal_c[clust$membership[i]]
}

##autonomous
autonom <- function(x, g){ ##x contains the union of concordance module and the balanced module
	#phase I
	f <- TRUE
	while (f)
	{
		res <- rep(1, length(x))
		for (i in 1:length(x))
		{
			s <- unlist(ego(g, 1, x[i], mode = "in"))
			check <- setdiff(s, x)
			if (length(check) > 0)
				res[i] <- 0 	
		}
		x <- x[which(res == 1)]
		if ((length(which(res == 0)) == 0) || (length(x) == 0))
			f <- FALSE
	}
	#phase II
	if (length(x) > 0)
	{
		dg <- degree(g, mode = "out")
		res <- rep(1, length(x))
		for (i in 1:length(x))
		{
			if (dg[x[i]] == 0)
				res[i] <- 0
		}
		x <- x[which(res == 1)]			
	}
	#phase III
	x1 <- x
	if (length(x) > 0)
		f <- TRUE
	while (f)
	{
		res <- rep(1, length(x))
		for (i in 1:length(x))
		{
			s <- unlist(ego(g, 1, x[i], mode = "out"))
			check <- setdiff(s, x)
			if (length(check) > 0)
				res[i] <- 0 	
		}
		x <- x[which(res == 1)]
		if ((length(which(res == 0)) == 0) || (length(x) == 0))
			f <- FALSE
	}
	phaseIII <- setdiff(x1, x)  
	#phase IV
	phaseIV <- c()
	if (length(x) > 0)
	{
		g_x <- induced_subgraph(g, x, impl = "copy_and_delete")
		clust <- components(g, mode = "strong")
	
		terminal_c <- rep(1, clust$no) #1 means terminal
		for (i in 1:clust$no){
			vc <- which(clust$membership == i)
			sn <- unlist(ego(g, 1, vc, mode = "out"))
			sdsn <- setdiff(sn, vc)
			if (length(sdsn) > 0)
				terminal_c[i] <- 0 		
		}

		terminal <- rep(1, length(clust$membership))
		for (i in 1:length(terminal))
		{
			terminal[i] <- terminal_c[clust$membership[i]]
		}
		phaseIV <- x[terminal[x] == 0] ##which complexes are non-terminal
	}
	return(union(phaseIII, phaseIV))
}

mets <- unlist(s$Results.balanced[,,1]$MODEL.r[[1]][[1]][,,1]$mets)
complexes <- unlist(s$Results.balanced[,,1]$MODEL.r[[1]][[1]][,,1]$complexes)

lres <- list()
# the balanced class is the one holding B; it need not be the largest class
maxj <- 1
if (!is.null(s$B))
{
	for (i in 1:length(s$class.with.balanced))
	{
		if (setequal(s$class.with.balanced[[i]][[1]][,1], as.vector(s$B)))
			maxj <- i
	}
}

for (i in 1:length(s$class.with.balanced))#i = maxj contains the balanced
{
	if (i != maxj)
	{
		possible <- union(s$class.with.balanced[[i]][[1]][,1], s$class.with.balanced[[maxj]][[1]][,1])
		p <- autonom(possible, g)
		#pnt <- intersect(p, which(terminal == 0))
		if (length(p) > 0)
			lres[[length(lres)+1]] <- p #pnt
	}
}

clustw <- components(g, "weak")
res_blc <- rep(0, max(clustw$membership)) #0 means that the linkage class is not composed of only balanced complexes
for (i in 1:max(clustw$membership))
{
	p <- which(clustw$membership == i)
	if (length(setdiff(p, s$class.with.balanced[[maxj]][[1]][,1])) == 0)
	{
		res_blc[i] <- 1
	}
}
blc_one <- which(res_blc == 1)
for (i in 1:length(blc_one))
{
	lres[[length(lres)+1]] <- which(clustw$membership == blc_one[i])
}

mdiff <- matrix(0, length(lres), length(lres))
for (i in 1:(length(lres) - 1))
{
	for (j in (i+1):length(lres))
	{
		mdiff[i, j] <- mdiff[j, i] <- length(intersect(lres[[i]], lres[[j]]))
	}
}

mg <- graph_from_adjacency_matrix(mdiff, mode = "undirected")
clustmg <- components(mg)

final_lres <- list()
complexes_single <- rep(1, length(complexes))
for (i in 1:clustmg$no)
{
	p <- which(clustmg$membership == i)
	final_lres[[length(final_lres) + 1]] <- lres[[p[1]]]
	complexes_single[lres[[p[1]]]] <- 0
	if (length(p) > 1)
	{
		for (j in 2:length(p))
		{
			final_lres[[length(final_lres)]] <- union(final_lres[[length(final_lres)]], lres[[p[j]]])
			complexes_single[lres[[p[j]]]] <- 0
		}
	}		
}
p <- which(complexes_single == 1)
for (i in 1:length(p))
{
	final_lres[[length(final_lres)+1]] <- p[i]	
}

##interface reactions
eg <- get.edges(g, es = E(g))
c12 <- c()
for (j in 1:length(final_lres))
{
	c1 <- which(eg[,1] %in% final_lres[[j]])
	c2 <- which(eg[,2] %in% final_lres[[j]])		
	c12 <- union(c12, intersect(c1, c2))
}


##look at results
#largest kinetic module
max_coupled_substrate <- max(unlist(lapply(final_lres, length)))
#substrate complexes
substrate_complexes <- length(complexes) - length(which(degree(g, mode = "out") == 0))
single_all <- length(which(unlist(lapply(final_lres, length)) == 1))
total_reactions <- ecount(g)
total_interface <- ecount(g) - length(c12)
deg <- degree(g, mode = "out")	
 
save(list = ls(), file = paste("../", gsub(".mat",".RData",k), sep = ""))
print('Done')

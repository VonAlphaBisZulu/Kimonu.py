library(R.matlab); library(igraph)
# the R code of the reference: a clone made by scripts/reference/fetch_upstream.sh (UPSTREAM_ALGORITHM). The test patches
# the original line, so run it on the original code: bash scripts/reference/fetch_upstream.sh d729161 <dir>, then
# UPSTREAM_ALGORITHM=<dir>; the fork's 8a6ed17 has the fix already (commit 28350a9).
U <- Sys.getenv("UPSTREAM_ALGORITHM", "scripts/reference/Upstream_Algorithm")
L <- readLines(file.path(U, "Code/kinetic_modules/code_kineticModule_analysis.R"))
a <- grep("^autonom <- function", L); b <- grep("^\treturn\\(union\\(phaseIII, phaseIV\\)\\)", L)
src <- L[a:(b+1)]
eval(parse(text = paste(src, collapse="\n")))                       # verbatim -> autonom
fixed <- sub("phaseIV <- x\\[terminal == 0\\]", "phaseIV <- x[terminal[x] == 0]", src)
fixed <- sub("^autonom <- function", "autonom_fixed <- function", fixed)
stopifnot(any(grepl("terminal\\[x\\]", fixed)))
eval(parse(text = paste(fixed, collapse="\n")))                     # one line changed

s <- readMat(commandArgs(TRUE)[1])
adjm <- s$Results.balanced[,,1]$MODEL.r[[1]][[1]][,,1]$A
g <- graph.empty(n = nrow(adjm), directed = TRUE)
for (j in 1:ncol(adjm)) g <- add_edges(g, c(which(adjm[,j]==-1), which(adjm[,j]==1)))
cwb <- s$class.with.balanced
sz <- sapply(cwb, function(x) length(x[[1]][,1])); maxj <- which.max(sz)
bal <- cwb[[maxj]][[1]][,1]

run <- function(fn, tag) {
  lres <- list(); nas <- 0
  for (i in seq_along(cwb)) { if (i==maxj) next
    p <- fn(union(cwb[[i]][[1]][,1], bal), g)
    nas <- nas + sum(is.na(p)); p <- p[!is.na(p)]
    if (length(p)>0) lres[[length(lres)+1]] <- p }
  n <- length(lres); l <- sapply(lres, length); uni <- length(unique(unlist(lres)))
  md <- matrix(0,n,n)
  for (i in 1:(n-1)) for (j in (i+1):n) md[i,j] <- md[j,i] <- length(intersect(lres[[i]],lres[[j]]))
  cm <- components(graph_from_adjacency_matrix(md, mode="undirected"))
  fl <- sapply(1:cm$no, function(i) length(unique(unlist(lres[which(cm$membership==i)]))))
  cat(sprintf("%-22s NAs %5d | sets %3d total %5d union %4d overlap %4d | blocks %4d largest %4d\n",
      tag, nas, n, sum(l), uni, sum(l)-uni, cm$no, max(fl)))
}
run(autonom,       "as published")
run(autonom_fixed, "terminal[x] == 0")
cat("Kimonu                 NAs     - | sets 257 total   594 union  594 overlap    0 | blocks  250 largest   19\n")

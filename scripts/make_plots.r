# Load the required packages
library(data.table)
library(qqman)

# Source function for qqplot (in case if the package GenABEL for qqplot not working, https://github.com/cran/GenABEL/blob/master/R/estlambda.R)
source("source_dir/estlambda.R")

# Function for qqplot
qqplotpval <- function(P, filename=NULL)
{
  l <- estlambda(P, method="median")
  nom <- paste("lambda = ", round(l$estimate, 3), sep="")
  if(!is.null(filename))
  {
    png(filename)
  }
  estlambda(P, method="median", plot=TRUE, main=nom)
  if(!is.null(filename))
  {
    dev.off()
  }
}

# Read file, the result of GWAS run after step 2
a <- fread("out.mlma")

# qqplot
qqplotpval(a$p, "out_qqplot.png")

# manhattan plot
png("out_manhattan.png",units = "in", width=12, height=5, res = 72)
manhattan(a, bp='bp', chr='Chr', p='p')
dev.off()

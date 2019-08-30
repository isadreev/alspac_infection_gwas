---
  title: Generating CMV data
author: Ildar Sadreev
date: 29/08/2019
---
  
  
  
  
  ## Objective
  
  - Read in the data available for CMV
- Maximise sample size by getting the value for each individual recorded at the most recent time point
- Create binary variable based on OD cutoff of 1
- Rank transform the continuous variable
- Obtain relevant covariates
- Write files for GWAS analysis



## Read in the data


```{r}
# open the phenotypic file
b <- get(load('File_pheno.Rdata'))

# Remove the withdrawn concent (NAs).
bb <- b[complete.cases(b), ]

# read the file with Gender to create covariates file
a <- read.table('Gender.txt', header = TRUE, sep = "\t")

# Merge the two files for Gender and phenotype based on ID
b1 <- merge(bb, a, by.bb= "aln", by.a = "aln")
```


How does OD compare against std_OD?
  
  
  ```{r}
plot(cmv_od_F7 ~ cmv_std_F7, data=b1)
```

How much similarity is there across age groups?
  
  
  ```{r}
plot(cmv_od_F7 ~ cmv_od_F11, data=b1)
plot(cmv_od_F7 ~ cmv_std_TF3, data=b1)
plot(cmv_od_F11 ~ cmv_std_TF3, data=b1)

```


## Selecting the optical density for the latest available age and assigning this age to covariates

There are three time points for children with CMV measures... How many from each time point?
  
  
  ```{r}
t=c()
n=c()
r=c()
q=c()

for (i in 1:nrow(b1)) {
  # od
  res1 = -1
  # age
  res2 = 0
  # od std
  res3 = -1
  # frequency of visits
  res4 = 0
  for (j in c("cmv_od_F7","cmv_od_F11","cmv_od_TF3")) {
    if (b1[i,j]>0) {
      res1 = b1[i,j]
      res4=res4+1
      if (j=="cmv_od_F7") {
        res2=7
        res3=b1[i,"cmv_std_F7"]
      }
      if (j=="cmv_od_F11") {
        res2=11
        res3=b1[i,"cmv_std_F11"]
      }
      if (j=="cmv_od_TF3") {
        res2=15
        res3=b1[i,"cmv_std_TF3"]
      }
    }
  }
  t[i]=res1
  n[i]=res2
  r[i]=res3
  q[i]=res4
}
c = cbind(b1,"od_max"=t)
c = cbind(c,"od_std_max"=r)
# Add threshold to obtain binary data for od std (od std >= 1 -> 1 else od std < 1 -> 0)
c = cbind(c,"od_std_max_bin"=r)
c[c[,"od_std_max_bin"]<1,"od_std_max_bin"] <- 0
c[c[,"od_std_max_bin"]>=1,"od_std_max_bin"] <- 1
# Add covariates: age max, age^2, age*gender (where 0 - female, 1 - male), age^2*gender), frequency of visits
c = cbind(c,"age_max"=n)
c = cbind(c,"age_max2"=n*n)
c = cbind(c,"age_max_gender"=n*(as.numeric(c$gender)-1))
c = cbind(c,"age_max2_gender"=n*n*(as.numeric(c$gender)-1))
c = cbind(c,"freq"=q)

# Remove non-complete cases (-1)
# This is the file with all the data
c1 <- c[c$od_max > 0, ]

```

How many records for each time point?
  
  ```{r}
nrow(c[c$cmv_od_F7>-1,])
nrow(c[c$cmv_od_F11>-1,])
nrow(c[c$cmv_od_TF3>-1,])

```


Calculate variance

```{r}
var(c1$od_std_max)
```

How many cases vs controls when data is binary?
  
  
  ```{r}
table(c1$od_std_max_bin)

```




Rank transform the data



```{r}
rank_transform <- function(x, s=NULL, m=NULL)
{
  if((is.null(s) & !is.null(m)) | (is.null(m) & !is.null(s)))
  {
    stop("s and m must both be null or not null")
  }
  out <- rank(x) - 0.5
  out[is.na(x)] <- NA
  mP <- 0.5/max(out, na.rm = T)
  out <- out/(max(out, na.rm = T) + 0.5)
  out <- qnorm(out)
  if(!is.null(s) & !is.null(m))
  {
    out <- out * s + m
  }
  out
}

c1$od_max_rt <- rank_transform(c1$od_max)
```

Plot histograms for std and rank transformed

```{r}

png("Pheno_std_rt.png") 
par(mfrow=c(2,1))
hist(c1$od_std_max, breaks=100, col="lightblue")
abline(v = 1, col="red", lwd=3, lty=2)
hist(c1$od_max_rt, breaks=100, col="lightblue")
dev.off()

```


## Write out the data

Use the format for GCTA: http://cnsgenomics.com/software/gcta/#GREMLanalysis
  
  We duplicate ID for each sample

Prepare the phenotype data for continious od


```{r}
c2 = c1[,c("alnqlet","alnqlet","od_max_rt")]

write.table(c2, file = "Pheno_cont.phen", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

# Prepare the phenotype data for binary od
c3 = c1[,c("alnqlet","alnqlet","od_std_max_bin")]

write.table(c3, file = "Pheno_bin.phen", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

```


Prepare covariates file with Gender (discrete)


```{r}
cv_1 = c1[,c("alnqlet","alnqlet","gender")]

write.table(cv_1, file = "Covariates_discrete.covar", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")



```

Prepare covariates file with age (quantitative)


```{r}
cv_2 = c1[,c("alnqlet","alnqlet","age_max","age_max2","age_max_gender","age_max2_gender","freq")]

write.table(cv_2, file = "Covariates_quant.qcovar", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
```


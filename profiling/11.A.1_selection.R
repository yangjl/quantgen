
parents <- rnorm(1000)
hist(parents, breaks=100)

### select top 5%
i =0.05
p <- head(sort(parents, decreasing=TRUE), round(length(parents)*i,0))
k <- (rep(p, each=length(p)) + rep(p, times=length(p)))/2

## R: response to selection, it is the difference of mean phenotypic value between the offspring of the selection
## parents and the whole of the parental generation before selection
R <- mean(k) - mean(parents)

## S: selection differential, it is the mean phenotypic values of the individuals selected as parents expressed
## as a deviation from the population mean
S <- mean(p) - mean(parents)

h2 <- R/S

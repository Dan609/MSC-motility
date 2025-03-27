library(readxl)
wj1 <- read_excel("~/Downloads/wj1.xlsx")
wj1 <- read_excel("~/Downloads/wj1_all.xlsx")
df2 <- read_excel("~/Downloads/df2.xlsx")
fetmsc <- read_excel("~/Downloads/FetMSC.xlsx")
xpa <- read_excel("~/Downloads/XPA.xlsx")
#View(wj1)
#str(wj1)


# Firstly, we conduct Hartigans’ dip test for unimodality or multimodality with dip.test() function available in diptest package (Maechler, 2021)
# Maechler, M. (2021). diptest: Hartigan’s Dip Test Statistic for Unimodality – Corrected. R package version 0.76-0.

data <- wj1[,1] # choose your column
data <- df2[,1]
data <- fetmsc[,4]
data <- xpa[,3]

data <- dplyr::pull(data)
data <- na.omit(data)


hist(data, # histogram
     breaks = 15,
     col="lightblue", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "Mean speed, microns per hour",
     main = "DF-2 passages  8") # MSCWJ-1 passages  6, 15, 26, 31
                              # DF-2 passages  8, 14, 21, 28
                               # FetMSC passages  12, 16, 20, 26
                                # XPA passages  3, 13, 19

lines(density(data), # density plot
      lwd = 2, # thickness of line
      col = "red")




library(diptest)
dip.test(data)

# Secondly, we use is.modal() function available in LaplacesDemon package (Statisticat, 2021).
# Statisticat, LLC. (2021). LaplacesDemon: Complete Environment for Bayesian Inference. Bayesian-Inference.com. R package version 16.1.6.

library(LaplacesDemon)
is.unimodal(data)
is.multimodal(data)
is.bimodal(data)
is.trimodal(data)
Modes(data)$modes















library(mousetrap)
bimodality_coefficient(data)

# The bimodality coefficient is obtained as 0.6014164. If this value is larger than 0.555 (stated in Pfister et al., 2013), it indicates bimodality of data.
# Pfister, R., Schwarz, K.A., Janczyk, M., Dale, R., Freeman, J.B. (2013). Good things peak in pairs: A note on the bimodality coefficient. Frontiers in Psychology, 4, 700.


#  we obtain the parameters of distributions for binomial data in R. We include cutoff package (Choisy, 2015) available in github. We need to install this package from github with install_github() function available in devtools package (Wickham, 2021).
# Wickham, H., Hester, J., Chang, W. (2021). devtools: Tools to Make Developing R Packages Easier. R package version 2.4.2.

# library(devtools)
# install_github("choisy/cutoff")
library(cutoff)

out <- em(data,"normal","normal")
confint(out)

# We use cutoff() function available in cutoff package (Choisy, 2015) to obtain bimodal data in R. It also returns confidence interval for the estimation.
cutoff(out)

# You could fit a Gaussian finite mixture model. Note that this makes the very strong assumption that your data are drawn from one or more true normals. As both @whuber and @NickCox point out in the comments, without a substantive interpretation of these data—supported by well-established theory—to support this assumption, this strategy should be considered exploratory as well.
# Now lets fit a Gaussian finite mixture model. In R, you can use the Mclust package to do this:
library(mclust)
x.gmm <-  Mclust(data)
summary(x.gmm)


x.gmm.1 = Mclust(data, G=1)
logLik(x.gmm.1)

logLik(x.gmm)-logLik(x.gmm.1)

1-pchisq(0, df=3)  # [1] 1.294187e-05
 
x.gmm$parameters
x.gmm.1$parameters

# https://stats.stackexchange.com/questions/138223/how-to-test-if-my-distribution-is-multimodal
# Some people don't feel comfortable using a parametric test here (although if the assumptions hold, I don't know of any problem). One very broadly applicable technique is to use the Parametric Bootstrap Cross-fitting Method (I describe the algorithm here). We can try applying it to these data:
# https://stats.stackexchange.com/questions/78973/zero-inflate-models-vs-generalized-mixture-model/78977#78977
set.seed(7809)
B = 10000;    x2.d = vector(length=B);    x1.d = vector(length=B)
for(i in 1:B){
  x2      = c(rnorm(68, mean=12346.98, sd=sqrt( 4514863)), 
              rnorm(52, mean=23322.06, sd=sqrt(24582180)) )
  x1      = rnorm( 120, mean=17520.91, sd=sqrt(43989870))
  x2.d[i] = Mclust(x2, G=2)$loglik - Mclust(x2, G=1)$loglik
  x1.d[i] = Mclust(x1, G=2)$loglik - Mclust(x1, G=1)$loglik
}
x2.d = sort(x2.d);  x1.d = sort(x1.d)
summary(x1.d)

# Formalizing this approach leads to the test given in Silverman (1981), "Using kernel density estimates to investigate modality", JRSS B, 43, 1. Schwaiger & Holzmann's silvermantest package implements this test, and also the calibration procedure described by Hall & York (2001), "On the calibration of Silverman's test for multimodality", Statistica Sinica, 11, p 515, which adjusts for asymptotic conservatism. Performing the test on your data with a null hypothesis of unimodality results in p-values of 0.08 without calibration and 0.02 with calibration. I'm not familiar enough with the dip test to guess at why it might differ.

# kernel density estimate for x using Sheather-Jones 
# method to estimate b/w:
x <- data
density(x, kernel="gaussian", bw="SJ") -> dens.SJ
# tweak b/w until mode just disappears:
density(x, kernel="gaussian", bw=3160) -> prox.null
# fill matrix with simulated samples from the proximal 
# null:
x.sim <- matrix(NA, nrow=length(x), ncol=10)
for (i in 1:10){
  x.sim[ ,i] <- rnorm(length(x), sample(x, size=length(x), 
                                        replace=TRUE), prox.null$bw)
}
# perform Silverman test without Hall-York calibration:
require(silvermantest)
silverman.test(x, k=1, M=10000, adjust=F)
# perform Silverman test with Hall-York calibration:
silverman.test(x, k=1, M=10000, adjust=T)

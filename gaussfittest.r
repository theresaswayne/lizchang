# gaussfittest.r
# test/demonstration of gaussian fitting using nls

require(tidyverse)
require(stringr) # for string harvesting

# using random data generated from a normal (Gaussian) density distribution
# x=seq(-4,4,0.01) 
# set.seed(123) # ensure identical results --- change seed to generate different numbers
# # normal distrib with mean of 0.4, sd of 2
# #y=2*dnorm(x, 0.4, 1.5) + runif( length(x) , min = -0.01, max = 0.01)
# y <- dnorm(x, mean=0.4 , sd=2)
# df=data.frame(x,y)
# #nls(y ~ k*dnorm(x, mu,sigma), data = df, start = list(k=2,mu=0.4,sigma=2))
# #df_fit <- nls(y ~ k*dnorm(x, mu,sigma), data = df, start = list(k=1,mu=0,sigma=1))
# df_fit <- nls(y ~ a + (b-a)*exp(-((x-c)^2)/(2*(d^2))), data = df, start = list(a=0, b=0.1, c=0.4, d=2))
# df_mean <- coef(df_fit)[3]
# p <- coef(df_fit)
# plot(df$x, df$y)
# curve(p["a"] + (p["b"]-p["a"])*exp(-((x-p["c"])^2)/(2*(p["d"]^2))), from = -4, to = 4, lwd=2, col="Red", add=TRUE)

# using real data

# specific file
#gp <- read_csv("Hapi/test output/otsu-2025926(17h06)/Masked (by Sum of ordered + disordered) GP images/Histograms/J774-CL90min-20uM-2chs.nd2GP Histogram(masked by Sum of ordered + disordered).csv")

# user chooses file each time
gpPath <- file.choose()
gp <- read_csv(gpPath)

# select the relevant columns
# select raw, normalized, smoothed, or unsmoothed counts
gpdf <- data.frame(x = gp$`GP values`, y = gp$`Counts (Smoothed Normalized)`)

# guess at the mean and other params 
# guesses for a, b, d based on fitting GP data in Fiji
# a, b will be proportional to the counts
mu_guess <- gpdf$x[which.max(gpdf$y)]
#a_guess <- 13000
#b_guess <- 300000
#d_guess <- 0.07
a_guess <- 0
b_guess <- 0.02
d_guess <- 0.2

#gp_fit <- nls(y ~ k*dnorm(x, mu,sigma), data = gpdf, start = list(k=1,mu=0,sigma=0.5))
# d is sqrt of variance, I think
gp_fit <- nls(y ~ a + (b-a)*exp(-((x-c)^2)/(2*(d^2))), data = gpdf, start = list(a=a_guess, b=b_guess, c=mu_guess, d=d_guess))
#gp_coeffs <- coef(gp_fit)
#gp_mean <- gp_coeffs[4]

# plot the fitted function in red, guessed function in blue, and the raw data
q <- coef(gp_fit)
plot(gpdf$x, gpdf$y)
curve(q["a"] + (q["b"]-q["a"])*exp(-((x-q["c"])^2)/(2*(q["d"]^2))), from = -1, to = 1, lwd=2, col="Red", add=TRUE)
curve(a_guess + (b_guess-a_guess)*exp(-((x-mu_guess)^2)/(2*(d_guess^2))), from = -1, to = 1, lwd=2, col="Blue",add=TRUE)
# retrieve the estimated mean (c) and its standard error
gp_mean <- summary(gp_fit)$coefficients[3,1]
gp_se <- summary(gp_fit)$coefficients[3,2]

# create a data row for output
# dummy for ROI name
# pattern J774-CL90min-5uM-2chs_Cyto_1_GP Histogram (masked by Sum of ordered + disordered)
gpFile <- basename(gpPath)
nameParsed <- str_split_1(gpFile, "_")
cellID <- nameParsed[3]
region <- nameParsed[2]

resultHeaders <- c("Filename", "Cell_ID", "Region", "GP_peak", "GP_StdErr")
resultValues <- c(gpFile, cellID, region, gp_mean, gp_se)
resultTable <- data.frame(rbind(resultHeaders, resultValues))
names(resultTable) <- resultTable[1,]
resultTable <- resultTable[-1,]

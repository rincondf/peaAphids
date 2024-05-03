Aphid_2022FCl <- read.csv("./Data/year2022.csv")
Aphid_2020FCl <- read.csv("./Data/year2020.csv")
Aphid_2019FCl <- read.csv("./Data/year2019.csv")
Aphid_2020FCl2 <- read.csv("./Data/year2020Val.csv") # Validation


require(MASS)
require(bbmle)

# Modeling only degree-days

plot(Aphid_2022FCl$AphidDDS, (Aphid_2022FCl$CumAphid), xlim = c(100, 1500), xlab = "Degree-days", ylab = "Cumulative proportion")
points(Aphid_2020FCl$Aphids_DD, (Aphid_2020FCl$CumAphid), col  = "blue")
points(Aphid_2019FCl$Aphids_DD, (Aphid_2019FCl$CumAphid), col  = "red")

degd <- c(Aphid_2022FCl$AphidDDS, Aphid_2020FCl$Aphids_DD, Aphid_2019FCl$Aphids_DD) # vector with degree-days for all years
PrAphid <- c((Aphid_2022FCl$PrAphid), (Aphid_2020FCl$PrAphid), (Aphid_2019FCl$PrAphid)) # vector with proportion of aphids (out of the total per site) per site per sampling time for all years
PrCAphid <- c((Aphid_2022FCl$CumAphid), (Aphid_2020FCl$CumAphid), (Aphid_2019FCl$CumAphid)) # vector with cumulative proportion of aphids (out of the total per site) per site per sampling time for all years


aphDDs <- rep(degd, round(PrAphid * 1000)) # variable repeating degree-days weighted by proportions of aphids. This is required for parameter estimation of a gamma pdf

modDDs <- fitdistr(aphDDs, "gamma")

MLL_fun <- function(sh, rt) {
  -sum(dgamma(dt, shape = sh, rate = rt, log = TRUE))
}

modelDDs <- mle2(MLL_fun, start = list(sh = coef(modDDs)[1], rt = coef(modDDs)[2]), data = list(dt = aphDDs))
summary(modelDDs)

lines(seq(200, 1400), pgamma(seq(200, 1400), shape = coef(modelDDs)[1], rate = coef(modelDDs)[2])) # adds the line with the gamma model
points(Aphid_2020FCl2$Aphids_DD, (Aphid_2020FCl2$CumAphid), col = "brown", lwd = 2) # adds the "validation" points (from pea crops, not included in parameter estimation)


# Modeling only day length

# Important: I had to tweak the daylength variable. The original one (dl) is circular. I'm using the cumulative difference between lengths in hours (dlD), so it does not go back after reaching the max daylength

plot(Aphid_2022FCl$dlD, (Aphid_2022FCl$CumAphid), xlim = c(14, 17), xlab = "Cumulative change in day length", ylab = "Cumulative proportion")
points(Aphid_2020FCl$dlD, (Aphid_2020FCl$CumAphid), col  = "blue")
points(Aphid_2019FCl$dlD, (Aphid_2019FCl$CumAphid), col  = "red")

light <- c(Aphid_2022FCl$dlD, Aphid_2020FCl$dlD, Aphid_2019FCl$dlD) # vector with degree-days for all years
aphlight <- rep(light, round(PrAphid * 1000)) # variable repeating daylengths weighted by proportions of aphids. This is required for parameter estimation of a gamma pdf

modlgt <- fitdistr(aphlight, "gamma")

modelLGT <- mle2(MLL_fun, start = list(sh = coef(modlgt)[1], rt = coef(modlgt)[2]), data = list(dt = aphlight))
summary(modelLGT)

lines(seq(14.5, 17, 0.1), pgamma(seq(14.5, 17, 0.1), shape = coef(modelLGT)[1], rate = coef(modelLGT)[2])) # adds the line with the gamma model
points(Aphid_2020FCl2$dlD, (Aphid_2020FCl2$CumAphid), col = "brown", lwd = 2) # adds the "validation" points (from pea crops, not included in parameter estimation)

require(piecewiseSEM)

data_complete <- data.frame(AphidC = PrAphid, cums = PrCAphid, DDs = degd, light  = light)

sem1 <- psem(glm(AphidC ~ DDs + light, data = data_complete, family = "binomial"))
summary(sem1)

sem2 <- psem(gam(cums ~ DDs + light, data = data_complete, family = "binomial"))
summary(sem2)


require(randomForest)

mod_RF <- randomForest(cums ~ DDs + light, data = data_complete, na.action = na.omit)
mod_RF
varImpPlot(mod_RF)

mod_RF1 <- randomForest(cums ~ DDs, data = data_complete, na.action = na.omit)
mod_RF1

mod_RF2 <- randomForest(cums ~ light, data = data_complete, na.action = na.omit)
mod_RF2


PC_mod <- prcomp(as.matrix(data_complete[, -c(1, 2)]), scale. = TRUE, retx = TRUE)
PC_mod
fviz_eig(PC_mod)
get_eigenvalue(PC_mod)


plot(PC_mod$x[, 1]+3, PrCAphid)

PC1_comb <- rep(PC_mod$x[, 1]+3, round(PrAphid * 1000))

modPC1 <- fitdistr(PC1_comb, "gamma")

modelComb <- mle2(MLL_fun, start = list(sh = coef(modPC1)[1], rt = coef(modPC1)[2]), data = list(dt = PC1_comb))
summary(modelComb)
lines(seq(0, 6, 0.1), pgamma(seq(0, 6, 0.1), shape = coef(modelComb)[1], rate = coef(modelComb)[2]))

anova(modelDDs, modelLGT, modelComb)
AICtab(modelDDs, modelLGT, modelComb)

Aphid_2022FCl <- read.csv("./Data/year2022.csv")
Aphid_2020FCl <- read.csv("./Data/year2020.csv")
Aphid_2019FCl <- read.csv("./Data/year2019.csv")
Aphid_2020FCl2 <- read.csv("./Data/year2020Val.csv")


require(MASS)

# Modeling only degree-days

plot(Aphid_2022FCl$AphidDDS, (Aphid_2022FCl$CumAphid), xlim = c(100, 1500), xlab = "Degree-days", ylab = "Cumulative proportion")
points(Aphid_2020FCl$Aphids_DD, (Aphid_2020FCl$CumAphid), col  = "blue")
points(Aphid_2019FCl$Aphids_DD, (Aphid_2019FCl$CumAphid), col  = "red")

degd <- c(Aphid_2022FCl$AphidDDS, Aphid_2020FCl$Aphids_DD, Aphid_2019FCl$Aphids_DD) # vector with degree-days for all years
PrAphid <- c((Aphid_2022FCl$PrAphid), (Aphid_2020FCl$PrAphid), (Aphid_2019FCl$PrAphid)) # vector with proportion of aphids (out of the total per site) per site per sampling time for all years

aphDDs <- rep(degd, round(PrAphid * 1000)) # variable repeating degree-days weighted by proportions of aphids. This is required for parameter estimation of a gamma pdf

modDDs <- fitdistr(aphDDs, "gamma")

lines(seq(200, 1400), pgamma(seq(200, 1400), shape = coef(modDDs)[1], rate = coef(modDDs)[2])) # adds the line with the gamma model
points(Aphid_2020FCl2$Aphids_DD, (Aphid_2020FCl2$CumAphid), col = "brown", lwd = 2) # adds the "validation" points (from pea crops, not included in parameter estimation)


# Modeling only day length

plot(Aphid_2022FCl$dlD, (Aphid_2022FCl$CumAphid), xlim = c(14, 17), xlab = "Cumulative change in day length", ylab = "Cumulative proportion")
points(Aphid_2020FCl$dlD, (Aphid_2020FCl$CumAphid), col  = "blue")
points(Aphid_2019FCl$dlD, (Aphid_2019FCl$CumAphid), col  = "red")

light <- c(Aphid_2022FCl$dlD, Aphid_2020FCl$dlD, Aphid_2019FCl$dlD) # vector with degree-days for all years
aphlight <- rep(light, round(PrAphid * 1000)) # variable repeating daylengths weighted by proportions of aphids. This is required for parameter estimation of a gamma pdf

modlgt <- fitdistr(aphlight, "gamma")

lines(seq(14.5, 17, 0.1), pgamma(seq(14.5, 17, 0.1), shape = coef(modlgt)[1], rate = coef(modlgt)[2])) # adds the line with the gamma model
points(Aphid_2020FCl2$dlD, (Aphid_2020FCl2$CumAphid), col = "brown", lwd = 2) # adds the "validation" points (from pea crops, not included in parameter estimation)


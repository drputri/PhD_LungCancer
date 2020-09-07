setwd("/mnt/exports/shared/home/dputri/PhD_LungCancer/")

load(file = "Data/ALLDATA.RData")

clin <- c("Group", "Age", "Packyears","BMI", 
          "COPD", "Cardiac", "Smoking", "Coagulation")
meta <- names(ALL.train)[grepl("var", names(ALL.train), ignore.case = TRUE)]
var <- c(clin, meta)

ALL.train[, "Smoking"] <- factor(ALL.train[, "Smoking"],
                                 levels = c("Never", "Active", "Stopped"))
table(ALL.train[, "Smoking"])
ALL.test[, "Smoking"] <- factor(ALL.test[, "Smoking"],
                                levels = c("Never", "Active", "Stopped"))
table(ALL.test[, "Smoking"])


# Create the matrix for lasso input
y <- ALL.train[, "Group"]
X <- ALL.train[var]
X <- subset(X, select = -c(Group))
X <- as.matrix.data.frame(X)

penalty.fctr <- rep(0, ncol(X))
names(penalty.fctr) <- colnames(X)
penalty.fctr[meta] <- 1


# Create response for test data
X.test <- ALL.test[var]
X.test <- subset(X.test, select = -c(Group))
X.test <- as.matrix.data.frame(X.test)

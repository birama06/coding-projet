#Project in Statistical computing
#Project Title: Development of Nonparametric Kruskal-Wallis Test in R
#Write a  R script to perform the nonparametric Kruskal Wallis test. The output variable is the probability value of the test.

kruskalwallis.test <- function(x, g) {
  
  if (length(x) != length(g)) {                                              # To check the length of x and g 
    stop("Length of data vector and grouping variable must be the same.")
  }
  
  # clean data to handle the NA
  complete_cases <- complete.cases(x, g) # This function is used to remove the rows which contain the NA values
  x <- x[complete_cases]
  g <- g[complete_cases]
  
  # Check if x is numeric and g is a factor
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (!is.factor(g)) { # if g is not a factor it will convert g as a factor
    g <- as.factor(g)
  }
  
  
  if (length(unique(g)) < 3) {   # In Kruskal-Wallis test is done if we have more than two samples
    stop("There must be at least three distinct groups in the grouping variable.")
  }
  
  # Calculation of ranks
  ranks <- rank(x)   # The rank of x is calculated using rank function 
  
  # Reorganization of ranks by groups
  ranks <- split(ranks, g)     # I allocated the rank of each value to it corresponding group 
  
  # Total number of observations
  N <- length(x)
  
  # Number of levels in the factor
  p <- length(unique(g))    # Total number of sample or group
  
  # Degrees of freedom
  df <- p - 1
  
  # Calculation of ties correction factor z
  tied <- table(x)                          # check for tied values 
  z <-1- (sum(tied^3 - tied) / (N^3 - N))    # the formula of correction factor z
  
  # Sum of ranks for each group
  sum_ranks <- sapply(ranks, sum)
  
  # Calculation of the Hobs statistic
  Hobs <- (12 / (N * (N + 1))) * sum((sum_ranks^2) / sapply(ranks, length))
  Hobs <- Hobs - 3 * (N + 1)
  
  # Adjust Hobs for tied values 
  if (any(tied > 1)) {
    Hobs <- Hobs / z
  }
  
  # Calculation of p-value
  p_value <- 1 - pchisq(Hobs, df = df)
  if (p_value>=0.05){
    print("H0 is accepted at 95 percent")
  }else{
    print("H0 is rejected at 0.05 level")
  }
  
  # Return the results
  cat("p_value =", p_value, "\n", "df =", df, "\n", "H =", Hobs,"\n")
}


#Example_1 of using my function

d <- read.csv("Class.csv", header = T)
d
attach(d)
names(d)
kruskalwallis.test(GR, CZ)

#Example_2
dat1  <- read.csv("Clas.csv", header = T)
dat1
attach(dat1)
names(dat1)
kruskalwallis.test(M, Z)


# Development of Nonparametric Kruskal-Wallis Test in R and MATLAB

## Project Overview
This project implements the Kruskal-Wallis H-test in R and MATLAB to perform non-parametric statistical testing. The function calculates the H statistic and p-value for a given dataset, allowing users to test the null hypothesis that samples originate from the same distribution.

The output variable is the p-value of the test, along with degrees of freedom and the H statistic.

---

## Features
- Handles missing data (NA values) automatically
- Accepts numeric vectors as data and categorical factors as grouping variables
- Checks for at least three groups, as required by Kruskal-Wallis
- Applies tie correction when necessary
- Provides clear output indicating whether the null hypothesis is accepted or rejected at the 5% significance level

---

## Function Usage (R)

```R
kruskalwallis.test <- function(x, g) {
  # Check length of x and g
  if (length(x) != length(g)) stop("Length of data vector and grouping variable must be the same.")
  
  # Remove missing values
  complete_cases <- complete.cases(x, g)
  x <- x[complete_cases]
  g <- g[complete_cases]
  
  # Ensure correct data types
  if (!is.numeric(x)) stop("x must be numeric.")
  if (!is.factor(g)) g <- as.factor(g)
  
  if (length(unique(g)) < 3) stop("There must be at least three distinct groups.")
  
  # Calculate ranks
  ranks <- rank(x)
  ranks <- split(ranks, g)
  
  N <- length(x)
  p <- length(unique(g))
  df <- p - 1
  
  # Ties correction
  tied <- table(x)
  z <- 1 - (sum(tied^3 - tied) / (N^3 - N))
  
  sum_ranks <- sapply(ranks, sum)
  Hobs <- (12 / (N * (N + 1))) * sum((sum_ranks^2) / sapply(ranks, length)) - 3 * (N + 1)
  
  if (any(tied > 1)) Hobs <- Hobs / z
  
  p_value <- 1 - pchisq(Hobs, df = df)
  
  if (p_value >= 0.05) {
    print("H0 is accepted at 95% confidence level")
  } else {
    print("H0 is rejected at 0.05 level")
  }
  
  cat("p_value =", p_value, "\n", "df =", df, "\n", "H =", Hobs, "\n")
}

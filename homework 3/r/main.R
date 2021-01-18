# Dejong
# Schwefel
# Rastrigin
# Michalewicz

# HC-FI
# HC-BI
# SA

outputFile = "output_ga"

n <- length(scan(file = outputFile, nlines = 1L, what = character()))
y <- scan(file = outputFile)

getResult <- function() {
  functionTypes = c( "br17", "ftv33", "ftv47", "ft53", "ftv70", "kro124", "ftv170", "rbg323", "rbg358", "rbg443")
  functionexpected = c(39, 1286, 1776, 6905, 1950, 36230, 2755, 1326, 1163, 2720)
  
  for (funct in 1:10) {
    
    print(functionTypes[funct])
    cat("Expected: ", functionexpected[funct], "\n")
    end = (funct - 1) * 30 + 30
    start = end - 9
  
    populationInfo(y[start:end])
    cat('\n')
    
  }
  
  cat('\n')
  
}

populationInfo <- function(c) {
  s = mean(c)
  d = sd(c)
  m = min(c)
  cat("Minimum: ", m, "\n")
  cat("Mean: ", s,"\n");
  cat("Standard Deviation: ", d, "\n")
  
}

getResult()





outputFileOne = "output_ga_one"

n <- length(scan(file = outputFileOne, nlines = 1L, what = character()))
y <- scan(file = outputFileOne)


plot(1:n, y, main="Cost Evolution(GA)", xlab="Generation", ylab="Cost", type = "l")

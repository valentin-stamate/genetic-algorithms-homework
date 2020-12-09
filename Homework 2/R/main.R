# Dejong
# Schwefel
# Rastrigin
# Michalewicz

# HC-FI
# HC-BI
# SA

outputFile = "output_aux.out"

n <- length(scan(file = outputFile, nlines = 1L, what = character()))
y <- scan(file = outputFile)

getResult <- function() {
  functionTypes = c("Dejong", "Schwefel", "Rastrigin", "Michalewicz")
  
  for (funct in 1:4) {
    
    print(functionTypes[funct])
    
    end = (funct - 1) * 10 + 10
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


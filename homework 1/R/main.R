# Dejong
# Schwefel
# Rastrigin
# Michalewicz

# HC-FI
# HC-BI
# SA

outputFile = "output_30.out"

n <- length(scan(file = outputFile, nlines = 1L, what = character()))
y <- scan(file = outputFile)

getResult <- function() {
  functionTypes = c("Dejong", "Schwefel", "Rastrigin", "Michalewicz")
  functionMethods = c("HillClimbing - Frist Improvement", "HillClimbing - Best Improvement", "Simulated Annealing")
  
  for (funct in 1:4) {
    
    print(functionTypes[funct])
    
    for (method in 1:3) {
      end = ( (funct - 1) * 3) * 30 + method * 30
      start = end - 29
      
      print(functionMethods[method])
      populationInfo(y[start:end])
      cat('\n')
    }
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

# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(ivcrtest)

test_check("ivcrtest")
# test_check_compatibility.R

# Load necessary packages
library(testthat)
library(wooldridge)
mydata <- wooldridge::card
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

d<-nrow(mydata)
attach(mydata)
x <- mydata$educ
y <- mydata$lwage
fitx=0
fity=0
fit=0

fity<-lm(y~exper +expersq  +black +south +smsa +smsa66+reg661+ reg662 +reg663+ reg664+ reg665+
           reg666 +reg667 +reg668, data=mydata)
y<-resid(fity)

fitx <-lm(x~exper +expersq  +black +south +smsa +smsa66+reg661+ reg662 +reg663+ reg664+ reg665+
            reg666 +reg667 +reg668, data=mydata)
x<-resid(fitx)

# ######################################
#  #########  prepare  IVs
z0=0
z<-fatheduc
fitz<-lm(z~exper +expersq  +black +south +smsa +smsa66+reg661+ reg662 +reg663+ reg664+ reg665+
           reg666 +reg667 +reg668, data=mydata)
z0<-resid(fitz)

k1<-min(length(z0),length(x),length(y))

z0<-predict(lm(z0[1:k1]~x[1:k1]+y[1:k1],data=mydata))

s <- min(length(y),length(x), length(z0) )
z <- z0[1:s]

df <- data.frame(x[1:s],y[1:s],z[1:s])
names(df) <- c("x","y","z")

# Run the function
test_result <- check_compatibility(
  df= df,
  n = n,
  k = -1,
  alpha = 0.05,
  rxu_range = c(0, 0.8) # smaller B for faster testing
)


# Print the result
print(test_result)

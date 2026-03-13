source("http://thegrantlab.org/misc/cdc.R")
View(cdc)
tail(cdc$weight, 20)
plot(cdc$height, cdc$weight,
     xlab = "Height (inches)",
     ylab = "Weight (pounds)")

cor(cdc$height, cdc$weight)

weight_kg <- cdc$weight * 0.454 
height_m <- cdc$height * 0.0254
bmi <- weight_kg / (height_m^2)

plot(cdc$height, bmi,
     xlab = "Height (inches",
     ylab = "BMI")

cor(cdc$height, bmi)

sum(bmi >= 30)

plot(cdc[1:100, "height"], cdc[1:100, "weight"],
     xlab="Height (inches)", ylab="Weight (pounds)")

obese <- bmi >= 30
table(cdc$gender[obese])


# Running plot tests:

print(getReqMassNaturalGas(100, -50, 969))

pList <- c()
m0List <- c()


for (p in seq(from = 0, to = 1200, by = 10)) {
  pList <- append(pList, p)
  m0List <- append(m0List, getReqMassNaturalGas(p, -50.0, 969.41))

}

plot(x = pList, y = m0List, main = "Medium Temp, Variable Pressure", xlab = "Pressure (psig)", ylab = "Projected Mass", type = "l")


# Reset x & y lists.

pList <- c()
mList <- c()

for (p in seq(from = 100, to = 1500, by = 10)) {
  pList <- append(pList, p)
  mList <- append(mList, getReqMassNaturalGas(p, 212.0, 969.41))
}

plot(x = pList, y = mList, main = "High Temp, Variable Pressure", xlab = "Pressure (psig)", ylab = "Projected Mass (kg)")

# Reset x & y lists.

tList <- c()
mList <- c()

for (t in seq(from = 10, to = 400, by = 10)) {
  tList <- append(tList, t)
  mList <- append(mList, getReqMassNaturalGas(100, t, 969.41))
}

plot(x = tList, y = mList, main = "Low Pressure, Variable Temp", xlab = "Temperature (F)", ylab = "Projected Mass (kg)")

# Reset x & y lists.

tList <- c()
mList <- c()

for (t in seq(from = 10, to = 400, by = 10)) {
  tList <- append(tList, t)
  mList <- append(mList, getReqMassNaturalGas(1200, t, 969.41))
}

plot(x = tList, y = mList, main = "High Pressure, Variable Temp", xlab = "Temperature (F)", ylab = "Projected Mass (kg)")




getY <- function(a, b, c, d, x) {

  return (a * (x ^ 3) + b * (x ^ 2) + c * x + d)

}

nextApprox <- function(a, b, c, d, interval) {
  x <- interval[1]
  minX <- interval[1]
  minY <- getY(a, b, c, d, minX)
  yPrev <- minY

  stepSize <- (interval[2] - interval[1]) / 10.0

  while (x < interval[2]) {
    y <- getY(a, b, c, d, x)
    if (sign(y) != sign(yPrev)) {
      minX <- x - stepSize
      minY <- yPrev
    }
    else if (abs(y) <= abs(minY)) {
      minX <- x
      minY <- y
    }
    yPrev <- getY(a, b, c, d, x)
    x <- x + stepSize
  }

  return (minX)

}

getCubicRootSimplified <- function(a, b, c, d) {

  interval <- 10
  minX <- 0

  while (interval >= 0.00000000000001) {
    minX <- nextApprox(a, b, c, d, c(minX, minX + interval))
    interval <- interval / 10
  }

#comment to compare newton's method effect

  return (minX)

}

getCubicRoot <- function(a, b, c, d) {

  interval <- 10
  minX <- 0

  while (interval >= 0.0000000000001) {
    minX <- nextApprox(a, b, c, d, c(minX, minX + interval))
    interval <- interval / 10
  }

#comment to compare newton's method effect
  tempMinX <- minX
  posMinX <- minX
  for (i in seq(from = 1, to = 10)) {
    x <- tempMinX
    y <- getY(a, b, c, d, x)
    dY <- getY(0, 3 * a, 2 * b, c, x)
    tempMinX <- x - (y / dY)
    if (tempMinX > minX - 0.001 & tempMinX < minX + 0.001) {
      posMinX <- tempMinX
    }
  }

  minX <- posMinX

  return (minX)

}

getMolarVolume <- function(p, t, pC, tC, w) {

  R <- 0.08206
  omegaA <- (9 * ((2 ^ (1 / 3)) - 1)) ^ (-1)
  omegaB <- (2 ^ (1 / 3) - 1) / 3
  tR <- t / tC
  a <- omegaA * (R ^ 2) * (tC ^ 2) / pC
  b <- omegaB * R * tC / pC
  alpha <- (1 + (0.48508 + 1.55171 * w - 0.15613 * (w ^ 2)) * (1 - (tR ^ 0.5))) ^ 2

  getP <- function(mV) {
    return ((R * t / (mV - b)) - (a * alpha / (mV * (mV + b))))
  }

  nextApprox <- function(interval) {
    x <- interval[1]
    closestX <- interval[1]
    closestY <- getP(closestX)
    yPrev <- closestY

    stepSize <- (interval[2] - interval[1]) / 10.0

    while (x < interval[2]) {
      y <- getP(x)
      if (sign(y - p) != sign(yPrev - p)) {
        closestX <- x - stepSize
        closestY <- yPrev
      }
      yPrev <- getP(closestX)
      x <- x + stepSize
    }
    print(getP(interval[1]))
    print(getP(interval[2]))

    return (closestX)
  }

  # Root finder
  iLen <- 10
  currentX <- b + 0.00000001

  while (iLen >= 0.000000001) {
    currentX <- nextApprox(c(currentX, currentX + iLen))
    iLen <- iLen / 10.0
  }

  mV <- currentX

  return (mV)
}

getReqMassNitrogen <- function(p, t, volume) {

  p <- p / 14.6959488
  t <- 5 * (t - 32) / 9
  t <- t + 273.15
  volume <- volume * 3.785

  pC <- 33.5
  tC <- 126
  molarMass <- 14.00674
  w <- 0.04

  molarVolume <- getMolarVolume(p, t, pC, tC, w)
  mass <- molarMass * volume / molarVolume
  return (mass / 1000)
}

getReqMassNaturalGas <- function(p, t, volume) {

  p <- p / 14.6959488
  t <- 5 * (t - 32) / 9
  t <- t + 273.15
  volume <- volume * 3.785

  pC <- 45.79
  tC <- 190.55
  molarMass <- 16.0426
  w <- 0.011

  molarVolume <- getMolarVolume(p, t, pC, tC, w)
  mass <- molarMass * volume / molarVolume
  return (mass / 1000)
}




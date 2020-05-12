# Root for Cech complex, general
CechGeneral1 = function(x) {
  l <- 1/2*(r-1+sqrt((1-x)^2-r^2)-x)
  rhol <- 1-sqrt((1-x)^2-r^2)+l
  rdeltacsq <- min(2*x, 1/2*(r^2-l^2+x*(2-x)+rhol*(2-rhol)))
  return(
    x + sqrt(
      r^2 - l^2 + x*(2-x) - ((1-x)^2 - r^2 + l^2 + (1-rhol)^2) *
      (1/sqrt(1-rdeltacsq) - 1)
    )-r
  )
}
CechGeneral2 = function(x) {
  rmin <- sqrt(
    1 - x*(2-x) - (2-r^2-x*(2-x))^2/4
  )
  rbsq <- 2*(rmin^2+x*(2-x)) / (1+sqrt(1-(rmin^2+x*(2-x))))
  rdeltabsq <- min(2*x, rbsq/2)
  return(
    (sqrt(
      rbsq - (2 - rbsq) * (1/sqrt(1-rdeltabsq) - 1)
    ) + 2*x)/sqrt(2) - rmin
  )
}
for (r in seq(from = 0.510, to = 0.514, by = 0.001)) {
  print(c(uniroot(f = CechGeneral1, interval = c(0, 0.1))[["root"]],
        uniroot(f = CechGeneral2, interval = c(0, 0.1))[["root"]]))
}


# Root for Cech complex, noiseless
CechNoiseless1 = function(x) {
  l <- 1/2*(r-1+sqrt(1-r^2)-x)
  rhol <- 1-sqrt(1-r^2)+l
  rdeltacsq <- min(x^2, 1/2*(r^2-l^2+rhol*(2-rhol)))
  return(
    x + sqrt(
      r^2 - l^2 - (1 - r^2 + l^2 + (1-rhol)^2) *
        (1/sqrt(1-rdeltacsq) - 1)
    )-r
  )
}
CechNoiseless2 = function(x) {
  rmin <- sqrt(
    1 - (2-r^2)^2/4
  )
  rbsq <- 2*(rmin^2) / (1+sqrt(1-rmin^2))
  rdeltabsq <- min(x^2, rbsq/2)
  return(
    (sqrt(
      rbsq - (2 - rbsq) * (1/sqrt(1-rdeltabsq) - 1)
    ) + 2*x)/sqrt(2) - rmin
  )
}
for (r in seq(from = 0.546, to = 0.55, by = 0.001)) {
  print(c(uniroot(f = CechNoiseless1, interval = c(0, 0.1))[["root"]],
        uniroot(f = CechNoiseless2, interval = c(0, 0.2))[["root"]]))
}


# Root for Cech complex, asymptotic
CechAsymptotic1 = function(x) {
  l <- 1/2*(r-1+sqrt((1-x)^2-r^2))
  rhol <- 1-sqrt((1-x)^2-r^2)+l
  rdeltacsq <- min(x*(2-x), 1/2*(r^2-l^2+x*(2-x)+rhol*(2-rhol)))
  return(
    sqrt(
      r^2 - l^2 + x*(2-x) - ((1-x)^2 - r^2 + l^2 + (1-rhol)^2) *
        (1/sqrt(1-rdeltacsq) - 1)
    )-r
  )
}
CechAsymptotic2 = function(x) {
  rmin <- sqrt(
    1 - x*(2-x) - (2-r^2-x*(2-x))^2/4
  )
  rbsq <- 2*(rmin^2+x*(2-x)) / (1+sqrt(1-(rmin^2+x*(2-x))))
  rdeltabsq <- min(x*(2-x), rbsq/2)
  return(
    (sqrt(
      rbsq - (2 - rbsq) * (1/sqrt(1-rdeltabsq) - 1)
    ) + 2*x)/sqrt(2) - rmin
  )
}
for (r in seq(from = 0.469, to = 0.473, by = 0.001)) {
  print(c(uniroot(f = CechAsymptotic1, interval = c(0, 0.1))[["root"]],
        uniroot(f = CechAsymptotic2, interval = c(0, 0.2))[["root"]]))
}


# Root for Vietoris-Rips complex, general
RipsGeneral1 = function(x) {
  r <- sqrt(1/2*(1-x)^2)
  rmin <- sqrt(
    1 - x*(2-x) - (2-r^2-x*(2-x))^2/4
  )
  rbsq <- 2*(r^2+x*(2-x)) / (1+sqrt(1-(r^2+x*(2-x))))
  rdeltabsq <- min(2*x, rbsq/2)
  return(
    (sqrt(
      rbsq - (2 - rbsq) * (1/sqrt(1-rdeltabsq) - 1)
    ) + 2*x) - 2*r
  )
}
RipsGeneral2 = function(x) {
  r <- sqrt(1/2*(1-x)^2)
  rmin <- sqrt(
    1 - x*(2-x) - (2-r^2-x*(2-x))^2/4
  )
  rbsq <- 2*(rmin^2+x*(2-x)) / (1+sqrt(1-(rmin^2+x*(2-x))))
  rdeltabsq <- min(2*x, rbsq/2)
  return(
    (sqrt(
      rbsq - (2 - rbsq) * (1/sqrt(1-rdeltabsq) - 1)
      ) + 2*x)/sqrt(2) - rmin
  )
}
uniroot(f = RipsGeneral1, interval = c(0, 0.2))
uniroot(f = RipsGeneral2, interval = c(0, 0.1))


# Root for Vietoris-Rips complex, noiseless
RipsNoiseless1 = function(x) {
  sqrt(1/2)-1/2*sqrt(
    1 / (2+sqrt(3))
  )-x
}
RipsNoiseless2 = function(x) {
  r <- sqrt(1/2*(1-x)^2)
  rmin <- sqrt(
    1 - (2-r^2)^2/4
  )
  rbsq <- 2*(rmin^2) / (1+sqrt(1-(rmin^2)))
  rdeltabsq <- min(x^2, rbsq/2)
  return(
    (sqrt(
      rbsq - (2 - rbsq) * (1/sqrt(1-rdeltabsq) - 1)
    ) + 2*x)/sqrt(2) - rmin
  )
}
uniroot(f = RipsNoiseless1, interval = c(0, 1))
uniroot(f = RipsNoiseless2, interval = c(0, 0.2))



# Root for Vietoris-Rips complex, asymptotic
RipsAsmptotic1 = function(x) {
  sqrt(1/2)*(1-x)-1/2*sqrt(
    2*(1/4*(1-x)^2+x*(2-x) / (1+sqrt(1-(1/4*(1-x)^2+x*(2-x)))) )
  )
}
RipsAsmptotic2 = function(x) {
  r <- sqrt(1/2*(1-x)^2)
  rmin <- sqrt(
    1 - x*(2-x) - (2-r^2-x*(2-x))^2/4
  )
  rbsq <- 2*(rmin^2+x*(2-x)) / (1+sqrt(1-(rmin^2+x*(2-x))))
  rdeltabsq <- min(x*(2-x), rbsq/2)
  return(
    (sqrt(
      rbsq - (2 - rbsq) * (1/sqrt(1-rdeltabsq) - 1)
    ) + 2*x)/sqrt(2) - rmin
  )
}
uniroot(f = RipsAsmptotic1, interval = c(0, 1))
uniroot(f = RipsAsmptotic2, interval = c(0, 0.1))



## Comparison

# Cech complex, general, from Niyogi et al. (2008)
3-sqrt(8)

# Cech complex, general, from Attali et al. (2013)
(-3+sqrt(22))/13

# Vietoris-Rips complex, general, from Attali et al. (2013)
(2*sqrt(2-sqrt(2))-sqrt(2))/(2+sqrt(2))
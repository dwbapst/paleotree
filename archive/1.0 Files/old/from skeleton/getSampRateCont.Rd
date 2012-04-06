\name{getSampRateCont}
\alias{getSampRateCont}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSampRateCont(timeData, n_tbins = 1, grp1 = NA, grp2 = NA, threshold = 0.1, est_only = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timeData}{
%%     ~~Describe \code{timeData} here~~
}
  \item{n_tbins}{
%%     ~~Describe \code{n_tbins} here~~
}
  \item{grp1}{
%%     ~~Describe \code{grp1} here~~
}
  \item{grp2}{
%%     ~~Describe \code{grp2} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
}
  \item{est_only}{
%%     ~~Describe \code{est_only} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (timeData, n_tbins = 1, grp1 = NA, grp2 = NA, threshold = 0.1, 
    est_only = F) 
{
    make_ft <- function(dur, FO, LO, n_tbins, grp1, grp2, g1s, 
        g2s) {
        qr_predict <- qr_predict_multpar(FO, LO, n_tbins, grp1, 
            grp2, g1s, g2s)
        function(par) {
            qr <- qr_predict(par)
            q <- qr[, 1]
            r <- qr[, 2]
            ft <- ifelse(dur == 0, q/(r + q), q * r * exp(-q * 
                dur)/(r + q))
            -sum(log(ft))
        }
    }
    qr_predict_multpar <- function(FO, LO, n_tbins, grp1, grp2, 
        g1s, g2s) {
        if (n_tbins > 1) {
            if (g1s > 0) {
                if (g2s > 0) {
                  function(par) {
                    n_tb <- n_tbins - 1
                    t_raw <- c(0.5, par[1:n_tb])
                    t_prop <- t_raw/sum(t_raw)
                    t_bl <- t_prop * (max(FO) - min(LO))
                    tbin <- c(max(FO), max(FO) - cumsum(t_bl))[1:n_tbins]
                    mqr <- matrix(par[(n_tb + 1):(length(par))], 
                      , 2, byrow = T)
                    mqrt <- matrix(mqr[1:n_tbins, ], , 2, byrow = T)
                    mqrg1 <- mqr[(n_tbins + 1):(n_tbins + g1s), 
                      ]
                    mqrg2 <- mqr[(n_tbins + g1s + 1):(n_tbins + 
                      g1s + g2s), ]
                    tcat <- sapply(FO, function(x) sum(tbin >= 
                      x))
                    qrt <- sapply(tcat, function(x) mqrt[x, ])
                    qrg1 <- sapply(grp1, function(x) mqrg1[x, 
                      ])
                    qrg2 <- sapply(grp2, function(x) mqrg2[x, 
                      ])
                    t((qrt + qrg1 + qrg2)/3)
                  }
                }
                else {
                  function(par) {
                    n_tb <- n_tbins - 1
                    t_raw <- c(0.5, par[1:n_tb])
                    t_prop <- t_raw/sum(t_raw)
                    t_bl <- t_prop * (max(FO) - min(LO))
                    tbin <- c(max(FO), max(FO) - cumsum(t_bl))[1:n_tbins]
                    mqr <- matrix(par[(n_tb + 1):(length(par))], 
                      , 2, byrow = T)
                    mqrt <- matrix(mqr[1:n_tbins, ], , 2, byrow = T)
                    mqrg1 <- mqr[(n_tbins + 1):(n_tbins + g1s), 
                      ]
                    tcat <- sapply(FO, function(x) sum(tbin >= 
                      x))
                    qrt <- sapply(tcat, function(x) mqrt[x, ])
                    qrg1 <- sapply(grp1, function(x) mqrg1[x, 
                      ])
                    t((qrt + qrg1)/2)
                  }
                }
            }
            else {
                function(par) {
                  n_tb <- n_tbins - 1
                  t_raw <- c(0.5, par[1:n_tb])
                  t_prop <- t_raw/sum(t_raw)
                  t_bl <- t_prop * (max(FO) - min(LO))
                  tbin <- c(max(FO), max(FO) - cumsum(t_bl))[1:n_tbins]
                  mqrt <- matrix(par[(n_tb + 1):(length(par))], 
                    , 2, byrow = T)
                  tcat <- sapply(FO, function(x) sum(tbin >= 
                    x))
                  tcount <- sapply(sort(unique(tcat)), function(x) sum(tcat == 
                    x))
                  t(sapply(tcat, function(x) mqrt[x, ]))
                }
            }
        }
        else {
            if (g1s > 0) {
                if (g2s > 0) {
                  function(par) {
                    mqr <- matrix(par, , 2, byrow = T)
                    mqrg1 <- mqr[1:g1s, ]
                    mqrg2 <- mqr[(g1s + 1):(g1s + g2s), ]
                    qrg1 <- sapply(grp1, function(x) mqrg1[x, 
                      ])
                    qrg2 <- sapply(grp2, function(x) mqrg2[x, 
                      ])
                    t((qrg1 + qrg2)/2)
                  }
                }
                else {
                  function(par) {
                    mqr <- matrix(par, , 2, byrow = T)
                    t(sapply(grp1, function(x) mqr[x, ]))
                  }
                }
            }
            else {
                function(par) {
                  matrix(par, length(FO), 2, byrow = T)
                }
            }
        }
    }
    FO <- timeData[, 1]
    LO <- timeData[, 2]
    dur <- (-LO) - (-FO)
    dur[dur < threshold] <- 0
    if (length(grp1) > 1) {
        g1s <- length(unique(grp1))
    }
    else {
        g1s <- 0
    }
    if (length(grp2) > 1) {
        g2s <- length(unique(grp2))
    }
    else {
        g2s <- 0
    }
    support_ft <- make_ft(dur, FO, LO, n_tbins, grp1, grp2, g1s, 
        g2s)
    npar <- (n_tbins - 1) + (ifelse(n_tbins > 1, 2 * n_tbins, 
        0) + ifelse(g1s > 0, g1s * 2, 0) + ifelse(g2s > 0, g2s * 
        2, 0))
    if (npar == 0) {
        npar <- 2
    }
    par_lim <- c(1e-04, 10)
    par_init <- rep(0.5, npar)
    par_min <- rep(par_lim[1], npar)
    par_max <- rep(par_lim[2], npar)
    answer <- optim(par_init, support_ft, method = "L-BFGS-B", 
        lower = par_min, upper = par_max, control = list(maxit = 1e+06))
    par <- answer$par
    if (est_only) {
        qr_predict <- qr_predict_multpar(FO, LO, n_tbins, grp1, 
            grp2, g1s, g2s)
        qr <- qr_predict(par)
        colnames(qr) <- c("qMax", "rMax")
        res <- qr
    }
    else {
        mes <- answer$message
        SMax <- (-answer$value)
        n <- length(dur)
        aicc <- (2 * npar) - (2 * SMax) + ((2 * npar * ((2 * 
            npar) + 1))/(n - (2 * npar) - 1))
        title <- paste("Analysis with", n_tbins, "time bins and", 
            sum(c(g1s > 0, g2s > 0)), "groupings (", g1s, "and", 
            g2s, "States),with", npar, "parameters")
        title <- c(title, "Note that par output is rate components, NOT avg rate")
        if (n_tbins > 1) {
            n_tb <- n_tbins - 1
            t_raw <- c(0.5, par[1:n_tb])
            t_prop <- t_raw/sum(t_raw)
            t_bl <- t_prop * (max(FO) - min(LO))
            t_ends <- c(max(FO), max(FO) - cumsum(t_bl))
            mqr <- matrix(par[(n_tb + 1):(length(par))], , 2, 
                byrow = T)
            Pp <- mqr[, 2]/(mqr[, 1] + mqr[, 2])
            mqr <- cbind(mqr, Pp)
            colnames(mqr) <- c("qMax", "rMax", "Comp")
            mqrt <- mqr[1:n_tbins, ]
            if (g1s > 0) {
                mqrg1 <- mqr[(n_tbins + 1):(n_tbins + g1s), ]
                if (g2s > 0) {
                  mqrg2 <- mqr[(n_tbins + 1 + g1s):(n_tbins + 
                    g1s + g2s), ]
                  res <- list(Title = title, SMax = SMax, AICc = aicc, 
                    mqrt = mqrt, mqrg1 = mqrg1, mqrg2 = mqrg2, 
                    t_ends = t_ends, binlen = t_bl, message = mes)
                }
                else {
                  res <- list(Title = title, SMax = SMax, AICc = aicc, 
                    mqrt = mqrt, mqrg1 = mqrg1, t_ends = t_ends, 
                    binlen = t_bl, message = mes)
                }
            }
            else {
                res <- list(Title = title, SMax = SMax, AICc = aicc, 
                  mqrt = mqrt, t_ends = t_ends, binlen = t_bl, 
                  message = mes)
            }
        }
        else {
            mqr <- matrix(par, , 2, byrow = T)
            Pp <- mqr[, 2]/(mqr[, 1] + mqr[, 2])
            mqr <- cbind(mqr, Pp)
            colnames(mqr) <- c("qMax", "rMax", "Comp")
            if (g1s > 0) {
                mqrg1 <- mqr[1:g1s, ]
                if (g2s > 0) {
                  mqrg2 <- mqr[(g1s + 1):(g1s + g2s), ]
                  res <- list(Title = title, SMax = SMax, AICc = aicc, 
                    mqrg1 = mqrg1, mqrg2 = mqrg2, message = mes)
                }
                else {
                  res <- list(Title = title, SMax = SMax, AICc = aicc, 
                    mqrg1 = mqrg1, message = mes)
                }
            }
            else {
                rMax <- mqr[, 2]
                qMax <- mqr[, 1]
                Pp <- rMax/(rMax + qMax)
                res <- list(Title = title, c(qMax = qMax, rMax = rMax, 
                  Comp = Pp), SMax = -SMax, AICc = aicc, message = mes)
            }
        }
    }
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

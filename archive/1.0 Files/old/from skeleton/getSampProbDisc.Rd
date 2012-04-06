\name{getSampProbDisc}
\alias{getSampProbDisc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSampProbDisc(timeData, n_tbins = 1, grp1 = NA, grp2 = NA, est_only = F, throttle = 1)
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
  \item{est_only}{
%%     ~~Describe \code{est_only} here~~
}
  \item{throttle}{
%%     ~~Describe \code{throttle} here~~
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
function (timeData, n_tbins = 1, grp1 = NA, grp2 = NA, est_only = F, 
    throttle = 1) 
{
    getPp <- function(R, q) {
        res <- numeric()
        for (T in 1:10000) {
            res[T] <- (1 - (1 - R)^T) * (exp(-q * (T - 1)) - 
                exp(-q * T))
        }
        sum(res)
    }
    ft_T <- function(R, q, t) {
        Pp <- sum(sapply(1:(15 * 10^throttle), function(T) (1 - 
            (1 - R)^T) * (exp(-q * (T - 1)) - exp(-q * T))))
        sum(sapply(t:(120 * throttle), function(T) (exp(-q * 
            (T - 1)) - exp(-q * T)) * (R^(ifelse(t > 2, 2, 1)) * 
            (T - t + 1) * ((1 - R)^(T - t))/Pp)))
    }
    make_ft <- function(dur, FO, LO, n_tbins, grp1, grp2, g1s, 
        g2s) {
        qR_predict <- qR_predict_multpar(FO, LO, n_tbins, grp1, 
            grp2, g1s, g2s)
        function(par) {
            qR <- qR_predict(par)
            q <- qR[, 1]
            R <- qR[, 2]/10
            ft <- sapply(1:length(dur), function(x) ft_T(R[x], 
                q[x], dur[x]))
            -sum(log(ft))
        }
    }
    qR_predict_multpar <- function(FO, LO, n_tbins, grp1, grp2, 
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
    dur <- apply(timeData, 1, diff) + 1
    timeData1 <- max(timeData) - timeData + 1
    FO <- timeData1[, 1]
    LO <- timeData1[, 2]
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
        lower = par_min, upper = par_max, control = list(maxit = 10000, 
            trace = 1))
    par <- answer$par
    if (est_only) {
        qR_predict <- qR_predict_multpar(FO, LO, n_tbins, grp1, 
            grp2, g1s, g2s)
        qR <- qR_predict(par)
        colnames(qR) <- c("qMax", "RMax")
        res <- qR
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
        title <- c(title, "Note that par output is rate/prob components, NOT avg rate/prob")
        if (n_tbins > 1) {
            n_tb <- n_tbins - 1
            t_raw <- c(0.5, par[1:n_tb])
            t_prop <- t_raw/sum(t_raw)
            t_bl <- t_prop * (max(FO) - min(LO))
            t_ends <- c(max(FO), max(FO) - cumsum(t_bl))
            mqR <- matrix(par[(n_tb + 1):(length(par))], , 2, 
                byrow = T)
            Pp <- sapply(1:nrow(mqR), function(x) getPp(mqR[x, 
                2], mqR[x, 1]))
            mqR <- cbind(mqR, Pp)
            colnames(mqR) <- c("qMax", "RMax", "Comp")
            mqRt <- mqR[1:n_tbins, ]
            if (g1s > 0) {
                mqRg1 <- mqR[(n_tbins + 1):(n_tbins + g1s), ]
                if (g2s > 0) {
                  mqRg2 <- mqR[(n_tbins + 1 + g1s):(n_tbins + 
                    g1s + g2s), ]
                  res <- list(Title = title, SMax = SMax, AICc = aicc, 
                    mqRt = mqRt, mqRg1 = mqRg1, mqRg2 = mqRg2, 
                    t_ends = t_ends, binlen = t_bl, message = mes)
                }
                else {
                  res <- list(Title = title, SMax = SMax, AICc = aicc, 
                    mqRt = mqRt, mqRg1 = mqRg1, t_ends = t_ends, 
                    binlen = t_bl, message = mes)
                }
            }
            else {
                res <- list(Title = title, SMax = SMax, AICc = aicc, 
                  mqRt = mqRt, t_ends = t_ends, binlen = t_bl, 
                  message = mes)
            }
        }
        else {
            mqR <- matrix(par, , 2, byrow = T)
            Pp <- sapply(1:nrow(mqR), function(x) getPp(mqR[x, 
                2], mqR[x, 1]))
            mqR <- cbind(mqR, Pp)
            colnames(mqR) <- c("qMax", "RMax", "Comp")
            if (g1s > 0) {
                mqRg1 <- mqR[1:g1s, ]
                if (g2s > 0) {
                  mqRg2 <- mqR[(g1s + 1):(g1s + g2s), ]
                  res <- list(Title = title, SMax = SMax, AICc = aicc, 
                    mqRg1 = mqRg1, mqRg2 = mqRg2, message = mes)
                }
                else {
                  res <- list(Title = title, SMax = SMax, AICc = aicc, 
                    mqRg1 = mqRg1, message = mes)
                }
            }
            else {
                RMax <- mqR[, 2]
                qMax <- mqR[, 1]
                Pp <- sapply(1:nrow(mqR), function(x) getPp(mqR[x, 
                  2], mqR[x, 1]))
                res <- list(Title = title, c(qMax = qMax, RMax = RMax, 
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

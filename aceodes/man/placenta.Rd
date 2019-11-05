\name{Placenta}
\alias{PlacentaStart}
\alias{PlacentaMod01}
\alias{PlacentaNael}
\alias{PlacentaNsel}
\alias{PlacentaEst01}
\alias{PlacentaLimits}
\alias{PlacentaFinal}
\alias{PlacentaTerminal}
\alias{PlacentaEvals}
\alias{PlacentaDefaultEvals}
\alias{PlacentaSolver}

\title{
Functions for the placenta example
}
\encoding{UTF-8}
\description{
Functions for reproducing the results from the placenta example in Section 5 of Overstall et al (2018). These functions implement approximations to four different utility functions; generate starting, terminal and final designs, provide grids of points satisfying 
constraints on the design to be used in conjunction with the \code{\link[acebayes]{ace}} function in the \code{\link[acebayes]{acebayes}} package, and generate values from the probabilistic solution.
}
\usage{
PlacentaStart(M)

PlacentaMod01(d, B)

PlacentaNael(d, B)

PlacentaNsel(d, B)

PlacentaEst01(d, B)

PlacentaLimits(i, j, d)

PlacentaTerminal(M, rep, utility = "NSEL", scaled = TRUE)

PlacentaEvals(M, rep = NULL, utility = "NSEL")

PlacentaDefaultEvals(M, utility = "NSEL")

PlacentaFinal(M, utility = "NSEL", scaled = TRUE)

PlacentaSolver(grid, x2u2, B, model.uncertainty = FALSE)

}

\arguments{
  \item{M}{
Number of placentas. In the example, in Overstall et al (2018), \code{M = 2,...,7}.
}
	\item{d}{
A \code{K} by \code{2} matrix giving the design. Each element is scaled to [-1,1] and \code{K = M + N/2} with \code{N = 8} being the number of time points. The first \code{M} rows control the concentrations of non-radioactive serine outside (first column) and inside (second column) the placenta. The remaining \code{N/2} rows give the \code{N} time points, where the first column gives the first \code{N/2} and the second column the last \code{N/2}. 
}
	\item{B}{
Monte Carlo sample size.
}
 \item{i,j}{
Arguments specifying the row (\code{i}) and column (\code{j}) of the design in current iteration of the ACE algorithm. These arguments are for the 
function \code{PlacentaLimits}. Typically, there is no need to call this function directly to reproduce the results in Overstall et al (2018).
}
\item{utility}{
Argument specifying the utility. It can be one of \code{"NSEL"} (the default), \code{"EST01"}, \code{"NAEL"} or \code{"MOD01"}.
}
\item{rep}{
Repetition number. ACE is run twenty times from twenty different starting designs. This argument specifies the repetition number and should be an integer from 1 to 20.
}
\item{x2u2}{
Vector length 2 giving the initial amounts of non-radioactive serine outside and inside the vesicle, at which to evaluate the probabilistic solver.
}
\item{grid}{
A grid on the interval [0,600]s at which to evaluate the probabilistic solver.
}
\item{model.uncertainty}{
Boolean argument specifying whether to account for model uncertainty in the \code{PlacentaSolver} function. If \code{TRUE}, then the 1st and 2nd reaction rates are the same. If \code{FALSE} (the default), then the reaction rates can be different.
}
\item{scaled}{
Boolean argument controlling whether the design returned by \code{placentafinal} is scaled (\code{TRUE}; the default) or not (\code{FALSE}). Designs are found by ACE having scaled the elements to the interval [-1,1]. By choosing \code{scaled = FALSE} the returned design will be on the original scale of [0,60]mins for time points and [0,1000] for concentrations of serine.
}
}
\details{
Functions for reproducing the results from the placenta example in Section 5 of Overstall et al (2018). These functions implement approximations to four different utility functions; generate starting, terminal and final designs, provide grids of points satisfying 
constraints on the design to be used in conjunction with the \code{\link[acebayes]{ace}} function in the \code{\link[acebayes]{acebayes}} package, and generate values from the probabilistic solution.

For more details on this example and the underlying methodology, see Overstall et al (2018). For more details on the ACE algorithm and its implementation in \code{acebayes}, see Overstall & Woods (2017) and
Overstall et al (2017), respectively.

}

\value{
The \code{PlacentaStart} function will return a \code{K} by \code{2} matrix giving the initial design where each element is scaled to [-1,1] and \code{K = M + N/2}. The first \code{M} rows control the concentrations of non-radioactive serine outside (first column) and inside (second column) the placenta. The remaining \code{N/2} rows give the \code{N} time points, where the first column gives the first \code{N/2} and the second column the last  \code{N/2}. This design will satisfy constraints on the design.

If the argument \code{scaled} is \code{TRUE}, then the \code{PlacentaTerminal} and \code{PlacentaFinal} functions will return \code{K} by \code{2} matrices giving the terminal and final designs, respectively, found by \code{\link[acebayes]{ace}} in Overstall et al (2018), where each element is scaled to [-1,1]. If the argument \code{scaled} is \code{FALSE}, the a list is returned with two elements: 1) an \code{M} by \code{2} matrix giving the concentrations of non-radioactive serine outside (first column) and inside (second column) the placenta; and 2) a vector of \code{N} elements giving the time points.

The functions \code{PlacentaMod01}, \code{PlacentaNael}, \code{PlacentaNsel} and \code{PlacentaEst01} return \code{B} Monte Carlo approximations to the Model 0-1, NAEL, NSEL and Estimation 0-1 utility functions, respectively.

If \code{rep} is specified then \code{PlacentaEvals} will return a vector of length \code{20} being the twenty double loop Monte Carlo approximations to the expected utility for the terminal designs found from one run of the ACE algorithm, as found by \code{\link[acebayes]{ace}} in Overstall et al (2018).  If \code{rep} is unspecified then \code{PlacentaEvals} will return a \code{20} by \code{20} matrix where the ith row is the same as a call to \code{PlacentaEvals} with rep = i. \bold{Note that \code{rep} should be an integer from 1 to 20, otherwise an error will be returned.}

The function \code{PlacentaDefaultEvals} will return a vector of length \code{20} being the twenty double loop Monte Carlo approximations to the expected utility for the original design for \code{M = 7} placentas.

The function \code{Placentalimits} will return a candidate grid of points over which the Gaussian process emulator to the expected utility will be maximised over the \code{i,j}th coordinate of the design.
}

\references{
Overstall, A.M. & Woods, D.C. (2017). Bayesian design of experiments using
approximate coordinate exchange. \emph{Technometrics}, \bold{59}, 458-470.

Overstall, A.M., Woods, D.C. & Adamou, M. (2017). acebayes: An R Package for Bayesian Optimal 
Design of Experiments via Approximate Coordinate Exchange. \emph{arXiv:1705.08096}.

Overstall, A.M., Woods, D.C. & Parker, B.M. (2018). Bayesian optimal design for ordinary differential 
equation models with application in biological science. \emph{arXiv:1509.04099}.
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods & Ben M. Parker 
}

\examples{
# Set M=7 placentas
M<-7

# Set a seed for reproducibility
set.seed(1)
# Create a starting design.
start.d<-PlacentaStart(M)

# Calculate a nested Monte Carlo approximation to the expected NSEL utility for 
# B = 200.
mean(PlacentaNsel(d=start.d,B=200))
}



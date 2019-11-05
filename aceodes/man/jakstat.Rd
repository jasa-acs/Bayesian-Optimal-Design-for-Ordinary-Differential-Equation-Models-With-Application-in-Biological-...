\name{Jakstat}
\alias{JakstatStart}
\alias{JakstatNael}
\alias{JakstatNsel}
\alias{JakstatSig}
\alias{JakstatLimits}
\alias{JakstatFinal}
\alias{JakstatTerminal}
\alias{JakstatEvals}
\alias{JakstatDefaultEvals}
\alias{JakstatSolver}

\title{
Functions for the JAK-STAT example
}
\encoding{UTF-8}
\description{
Functions for reproducing the results from the JAKSTAT example in Section 4.4 of Overstall et al (2018). These functions implement approximations to three different utility functions; generate starting, terminal and final designs, provide grids of points satisfying 
constraints on the design to be used as the \code{limits} argument of the \code{\link[acebayes]{ace}} function in the \code{\link[acebayes]{acebayes}} package, reproduce approximations to the expected utility for all terminal designs and original design under each utility, and generate values from the probabilistic solution.
}
\usage{
JakstatStart()

JakstatNael(d, B)

JakstatNsel(d, B)

JakstatSig(d, B)

JakstatLimits(i, j, d)

JakstatTerminal(rep, utility = "SIG", scaled = TRUE)

JakstatEvals(rep = NULL, utility = "SIG")

JakstatDefaultEvals(utility = "SIG")

JakstatFinal(utility = "SIG", scaled = TRUE)

JakstatSolver(grid, B)

}

\arguments{
	\item{d}{
An \code{17} by \code{1} matrix giving the design. Each element is scaled to [-1,1]. The first row corresponds to the time at which \eqn{G_4} is measured and the remaining \code{16} rows correspond to the times at which \eqn{G_1} and \eqn{G_2} are measured.
}
	\item{B}{
Monte Carlo sample size.
}
 \item{i,j}{
Arguments specifying the row (\code{i}) and column (\code{j}) of the design in current iteration of the ACE algorithm. These arguments are for the 
function \code{JakstatLimits}. Typically, there is no need to call this function directly to reproduce the results in Overstall et al (2018).
}
\item{rep}{
Repetition number. ACE is run twenty times from twenty different starting designs. This argument specifies the repetition number and should be an integer from 1 to 20.
}
\item{utility}{
Argument specifying the utility. It can be one of \code{"SIG"} (the default), \code{"NSEL"} or \code{"NAEL"}.
}
\item{scaled}{
Boolean argument controlling whether the design returned by \code{jakstatfinal} is scaled (\code{TRUE}; the default) or not (\code{FALSE}). Designs are found by ACE having scaled the elements to the interval [-1,1]. By choosing \code{scaled = FALSE} the returned design will be on the original scale of [0,60]s.
}
\item{grid}{
A grid on the interval [0,60]s at which to evaluate the probabilistic solver.
}
}
\details{
Functions for reproducing the results from the JAKSTAT example in Section 4.4 of Overstall et al (2018). These functions implement approximations to three different utility functions; generate starting, terminal and final designs, provide grids of points satisfying 
constraints on the design to be used as the \code{limits} argument of the \code{\link[acebayes]{ace}} function in the \code{\link[acebayes]{acebayes}} package, reproduce approximations to the expected utility for all terminal designs and original design under each utility, and generate values from the probabilistic solution.

For more details on this example and the underlying methodology, see Overstall et al (2018). For more details on the ACE algorithm and its implementation in \code{acebayes}, see Overstall & Woods (2017) and
Overstall et al (2017), respectively.

}

\value{
The \code{JakstatStart} function will return a \code{17} by \code{1} matrix giving an initial design where each element is scaled to [-1,1]. The first row corresponds to the time at which \eqn{G_4} is measured and the remaining \code{16} rows correspond to the times at which \eqn{G_1} and \eqn{G_2} are measured. This design will satisfy constraints on the design.

The \code{JakstatTerminal} and \code{JakstatFinal} functions will return \code{17} by \code{1} matrices giving the terminal and final designs, respectively, found by \code{\link[acebayes]{ace}} in Overstall et al (2018), where each element is scaled to [-1,1] (if \code{scaled = TRUE}).

If \code{rep} is specified then \code{JakstatEvals} will return a vector of length \code{20} being the twenty double loop Monte Carlo approximations to the expected utility for the terminal designs found from one run of the ACE algorithm, as found by \code{\link[acebayes]{ace}} in Overstall et al (2018).  If \code{rep} is unspecified then \code{JakstatEvals} will return a \code{20} by \code{20} matrix where the ith row is the same as a call to \code{JakstatEvals} with rep = i. \bold{Note that \code{rep} should be an integer from 1 to 20, otherwise an error will be returned.}

The function \code{JakstatDefaultEvals} will return a vector of length \code{20} being the twenty double loop Monte Carlo approximations to the expected utility for the original design.

The functions \code{JakstatNael}, \code{JakstatNsel} and \code{JakstatSig} return \code{B} Monte Carlo approximations to the NAEL, NSEL and SIG utility functions, respectively.

The function \code{JakstatLimits} will return a candidate grid of points over which the Gaussian process emulator to the expected utility will be maximised over the \code{i,j}th coordinate of the design.
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
# Set a seed for reproducibility
set.seed(1)
# Create a starting design.
start.d<-JakstatStart()

# Calculate a nested Monte Carlo approximation to the expected SIG utility for 
# B = 200.
mean(JakstatSig(d=start.d,B=200))
}



\name{Compartmental}
\alias{CompartmentalNaelProb}
\alias{CompartmentalNselProb}
\alias{CompartmentalSigProb}
\alias{CompartmentalNaelExact}
\alias{CompartmentalNselExact}
\alias{CompartmentalSigExact}
\alias{CompartmentalLimits}
\alias{CompartmentalTerminal}
\alias{CompartmentalFinal}
\alias{CompartmentalEvals}

\title{
Functions for the compartmental model example
}
\encoding{UTF-8}
\description{
Functions for reproducing the results from the compartmental model example in Section 4.2 of Overstall et al (2018). These functions implement approximations to three different utility functions; generate final and terminal designs, provide grids of points satisfying 
constraints on the design to be used as the \code{limits} argument of the \code{\link[acebayes]{ace}} function in the \code{\link[acebayes]{acebayes}} package, and reproduce approximations to the expected utility for all terminal designs under each utility and solver.
}
\usage{
CompartmentalNaelProb(d, B)

CompartmentalNaelExact(d, B)

CompartmentalNselProb(d, B)

CompartmentalNselExact(d, B)

CompartmentalSigProb(d, B)

CompartmentalSigExact(d, B)

CompartmentalLimits(i, j, d)

CompartmentalTerminal(rep, utility = "SIG", scaled = TRUE, solver = "probabilistic")

CompartmentalEvals(rep, utility = "SIG", solver = "probabilistic")

CompartmentalFinal(utility = "SIG", scaled = TRUE, solver = "probabilistic")

}

\arguments{
	\item{d}{
A \code{15} by \code{1} matrix giving the design. Each element is scaled to [-1,1] and corresponds to a time point. 
}
	\item{B}{
Monte Carlo sample size.
}
 \item{i,j}{
Arguments specifying the row (\code{i}) and column (\code{j}) of the design in current iteration of the ACE algorithm. These arguments are for the 
function \code{CompartmentalLimits}. Typically, there is no need to call this function directly to reproduce the results in Overstall et al (2018).
}
\item{utility}{
Argument specifying the utility. It can be one of \code{"SIG"} (the default), \code{"NSEL"} or \code{"NAEL"}.
}
\item{rep}{
Repetition number. ACE is run twenty times from twenty different starting designs. This argument specifies the repetition number and should be an integer from 1 to 20.
}
\item{scaled}{
Boolean argument controlling whether the designs returned by \code{CompartmentalTerminal} and \code{CompartmentalFinal} are scaled (\code{TRUE}; the default) or not (\code{FALSE}). Designs are found by ACE having scaled the elements to the interval [-1,1]. By choosing \code{scaled = FALSE} the returned designs will be on the original scale of [0,24] hours.
}
\item{solver}{
Argument controlling whether the designs returned by \code{CompartmentalTerminal} and \code{CompartmentalFinal} or evaluations of the expected utility were under the probabilistic solution to the ODEs (\code{solver = "probabilistic"}) or exact (\code{solver = "exact"}).
}
}
\details{
Functions for reproducing the results from the compartmental model example in Section 4.2 of Overstall et al (2018). These functions implement approximations to three different utility functions; generate final and terminal designs, provide grids of points satisfying 
constraints on the design to be used as the \code{limits} argument of the \code{\link[acebayes]{ace}} function in the \code{\link[acebayes]{acebayes}} package, and reproduce approximations to the expected utility for all terminal designs under each utility and solver.

For more details on this example and the underlying methodology, see Overstall et al (2018). For more details on the ACE algorithm and its implementation in \code{\link[acebayes]{acebayes}}, see Overstall & Woods (2017) and
Overstall et al (2017), respectively.
}

\value{
The functions \code{CompartmentalNaelExact}, \code{CompartmentalNselExact} and \code{CompartmentalSigExact} return \code{B} Monte Carlo approximations to the NAEL, NSEL and SIG utility functions, respectively, under the exact solution.

The functions \code{CompartmentalNaelProb}, \code{CompartmentalNselProb} and \code{CompartmentalSigProb} return \code{B} Monte Carlo approximations to the NAEL, NSEL and SIG utility functions, respectively, under the probabilistic solution.

The function \code{CompartmentalLimits} will return a candidate grid of points over which the Gaussian process emulator to the expected utility will be maximised over the \code{i,j}th coordinate of the design.

If \code{rep} is specified then \code{CompartmentalEvals} will return a vector of length \code{20} being the twenty double loop Monte Carlo approximations to the expected utility for the terminal designs found from one run of the ACE algorithm, as found by \code{\link[acebayes]{ace}} in Overstall et al (2018).  If \code{rep} is unspecified then \code{CompartmentalEvals} will return a \code{20} by \code{20} matrix where the ith row is the same as a call to \code{CompartmentalEvals} with rep = i. \bold{Note that \code{rep} should be an integer from 1 to 20, otherwise an error will be returned.}

The \code{CompartmentalTerminal} and \code{CompartmentalFinal} functions will return \code{15} by \code{1} matrices giving the terminal and final design, respectively, found by \code{\link[acebayes]{ace}} in Overstall et al (2018), where each element is scaled to [-1,1] (if \code{scaled = TRUE}).
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
start.d<-matrix(runif(15)*2-1,ncol=1)

# Calculate a nested Monte Carlo approximation to the expected SIG utility for 
# B = 200.
mean(CompartmentalSigProb(d=start.d,B=200))
}



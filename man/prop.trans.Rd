% Generated by roxygen2 (4.0.2): do not edit by hand
\name{prop.trans}
\alias{prop.trans}
\title{Transformation of Proportion Data (prop.trans)}
\usage{
prop.trans(your.object)
}
\arguments{
\item{Numeric}{column position x in the dataframe}
}
\value{
Returns transformed column vector
}
\description{
Freeman Tukey Arcsin Transformation For Proportions
}
\details{
This technique is widely used to Stabilize Variance during model building.  The proportion should be arrived from a COUNT data only.  It should not be applied under scenarios like where the Proportion is arrived from Continuous Data.  Should not be applied for Profit Margins.
}


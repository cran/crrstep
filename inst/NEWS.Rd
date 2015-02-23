\name{NEWS}
\alias{NEWS}

\title{New for Package \pkg{crrstep}}

\section{Changes in version 2014-07.16}{
  \itemize{
  \item Corrected a bug in backward selection procedure when one variable name is   contained in another variable name (e.g. "Mitosis" and "AtyMitosis").
  \item Correceted a bug in forward selection procedure when a as.factor() variable is to be added to the model.   
  }
}

\section{Changes in version 2015-2.1}{
  \itemize{
  \item Corrected a bug in backward selection procedure in how scope.min was checked to terminate the selection procedure.
  \item Added checks to make sure that `scope.min' formula input is correct
  \item Simplified how std error was computed
  \item Fixed a bug related to outputting results in forward selection
  }
}


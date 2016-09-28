\name{Harmannews}
\title{Harman News}
\encoding{UTF-8}

\section{Version 1.2.0}{\itemize{

\item Dynamically resizes legends depending upon the number of batches in
pcaPlot()
\item A custom prcomp function to get the appropriate scores. The standard R
`prcomp` function did not work in instances where the number of samples
(matrix columns) was greater than the number of assays (matrix rows), so a
special case was created for less assays than samples.
In this case we need to use u' instead of v from the SVD.
\item prcompPlot now has an argument for scaling, which defaults to FALSE.
Previously, scaling was always TRUE. This new default makes prcompPlots plots
agree with the 'original' plots of plotting harmanresults objects.
\item Extensive updating of the vignette with a new comparison to ComBat.
}}


\section{Version 1.0.2}{\itemize{

\item
First public version on Bioconductor.

}}

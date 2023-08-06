# ImputeBench

Author: Robin Richter

Institution: Deutsches Zentrum für Neurodegenerative Erkrankungen (DZNE) e.V.

Licence: MIT Licence

Note: Please feel free to ask questions or report bugs.

The ImputeBench package includes two functions to benchmark imputation methods either on provided, or, on simulated data. The benchmarking, data simulation
and missingness pattern simulation protocols have been described in the accompanying paper 

R. Richter, J.F. Tavares, A. Miloschewski, M.M.B. Breteler and S. Mukherjee (2022) ImputeBench - Benchmarking Single Imputation Methods, submitted.

Moreover, ImputeBench includes a number of plot functions to visualize observed, missing and imputed data distributions using ggplot2 or plotly. The Vignette.Rmd file can be knit into the vignette of ImputeBench (to adjust the time of knitting change 
64  sz = 500
65  sz.2 = 500
66  reps = 10
in the .Rmd-file.).

The functionality of ImputeBench is explained in detail in the accompanying vignette that is available in this repository. The compiled vignette is available under https://figshare.com/articles/online_resource/ImputeBench_Vignette_html/23896677.

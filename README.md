semper
======

### What is this?

### Install

From CRAN (though, paper is not yet on CRAN)

    install.packages("semper")

Devel version

    install.packages("devtools")
    devtools::install_github("brandmaier/semper")

    library("semper")
    #> Loading required package: parallel
    #> Loading required package: MASS
    #> Loading required package: OpenMx
    #> OpenMx is not compiled to take advantage of computers with multiple cores.
    #> 
    #> Attaching package: 'semper'
    #> The following object is masked from 'package:stats':
    #> 
    #>     power

### Examples

    require("semper")
    model <- lgcm(timepoints = 0:4,
     intercept.variance = 10,               
     slope.variance = 1,residual.variance = 50)

    cat("ECR: ", ecr(model),"\n") 
    #> ECR:  0.2857143
    cat("Effective error: ", 
    effectiveError(model),"\n") 
    #> Effective error:  2.5
    cat("GRR: ", grr(model))
    #> GRR:  0.1666667

![](./inst/flag_yellow_low.jpg) | ABC

This project has received funding from the European Union’s Horizon 2020
research and innovation programme under grant agreement No [number].

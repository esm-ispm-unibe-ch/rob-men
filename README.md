ROB-MEN
=======
ROB-MEN is a web application for evaluating the risk of bias due to missing evidence in network meta-analysis built with R Shiny and Bootstrap.

Currently, the app is available [online](https://cinema.ispm.unibe.ch/rob-men/), thanks to funding from the Swiss National Science Foundation. It's also possible to run ROB-MEN locally in the following way.

Requirements
------------
* [R 4.2.3](https://www.r-project.org/) or above
* [BUGSnet](https://bugsnetsoftware.github.io/)
* [JAGS](https://sourceforge.net/projects/mcmc-jags/)
* R packages: `BUGSnet`, `bslib`, `data.table`, `dplyr`, `DT`, `meta@5.2-0`, `netmeta@2.1-0`, `R2jags`, `shiny`, `tidyverse`
* modern web browser

Usage
-----
0. Install all required software (specifically, respecting the recommended versions of `meta` and `netmeta`)
1. Run R from the root directory of the project
2. Execute the following:
    ```
    load(shiny)
    runApp()
    ```
3. The app will open directly in your browser (or in the built-in browser for RStudio users)

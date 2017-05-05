-   [Intalling rROMADash](#intalling-rromadash)
-   [Using rROMADash](#using-rromadash)
-   [Using the docker container](#using-the-docker-container)

This package provides a shiny interface that can be used to complement
the [rRoma R package](https://github.com/Albluca/rRoma).

Intalling rROMADash
===================

The rROMADash package relies on the `rRoma` package, which needs to be
installed via `devtools`

    if(!requireNamespace("devtools")){
      install.packages("devtools")
    }
    devtools::install_github("Albluca/rRoma")

After installing `rRoma`, it is possible to install `rROMADash` by
typing

    if(!requireNamespace("devtools")){
      install.packages("devtools")
    }
    devtools::install_github("Albluca/rRomaDash")

The other packages required by the interface will be installed
automatically by R if not already available on the system.

Using rROMADash
===============

To launch the interface it is sufficient to type to following command:

    library(rRomaDash)
    rRomaDash()

It is also possible to enable interactivity (via the `plotly` package),
by typing

    library(rRomaDash)
    rRomaDash(Interactive = TRUE)

Note that the interactivity is still experimental and may not work
properly on certain systems.

Using the docker container
==========================

A docker image containing the RStusio web server and all the necessary
packages needed to run rRoma and rRomaDash is also available. After
installing [docker](http://www.docker.com), the image can be obtained
using

    docker pull albluca/rroma

After installation, the image can be started with the command

    docker run -p 8787:8787 albluca/rroma

At this point, the RStudio web interface will be availabe. In Unix-like
systems, the interface will be available at the address
<http://localhost:8787>. It is then possible to access the interface by
using "rstudio" as the username and password.

On windows, it may be necessary to replace localhost with the address of
the machine. If any problem is encountered see the instruction
[here](https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image).

After starting the interface, it is possible to write r code as in the
desktop version of RStudio.

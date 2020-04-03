# snphub

SnpHub is a Shiny-based server framework for retrieving, analyzing and visualizing the large genomic variations data in a lab.

[Our homepage](http://guoweilong.github.io/SnpHub/)

For more **details**, check [our tutorial](https://esctrionsit.github.io/snphub_tutorial/)

For **Docker-encapsulated version**, check [here](https://github.com/esctrionsit/snphub4docker)

To get a **quick deploy**, see [here](https://esctrionsit.github.io/snphub_tutorial/content/Setup/quick_deploy.html)

[Here](http://wheat.cau.edu.cn/Wheat_SnpHub_Portal/) are our **live demos**.

![](SnpHub.jpg)

## Environment request

To run the SnpHub, make sure the following software programs are **already** installed:
- samtools (≥ **v1.4**)
- bcftools (≥ **v1.8**)
- seqkit
- tabix (≥ **v1.6**)

Also, the following R packages are also prerequisites:
- ggplot2
- crayon
- ggmap
- dplyr
- rjson
- shiny
- pegas
- maps
- vcfR
- ape
- DT

## Demo

Two steps are needed to run on demo data set:

- Using `git` to clone SnpHub to local.
```sh
git clone https://github.com/esctrionsit/snphub
```

- Using `R` to run on demo data set. *(Or copy the SnpHub into your `shiny-server app folder`)*
```sh
R -e "shiny::runApp('./snphub', port=5000, host='0.0.0.0')"
```

*Make sure that `samtools`, `bcftools`, `seqkit` and `tabix` are added to system PATH. Otherwise, `tool application paths` part in `setup.conf` is needed to change to fit.*

## Setup on own data set

### Fulfill config file

Basic config file is named `setup.conf`, fulfill it and then run command:

``` sh
Rscript ./setup.R
```

Check [configuration](https://esctrionsit.github.io/snphub_tutorial/content/Setup/configuration.html) for more details about `setup.conf`.

Check [file format](https://esctrionsit.github.io/snphub_tutorial/content/Setup/file-formats.html) for more details about required file formats.

If it is your **first time** to set up `SnpHub`, you need to **delete** `advanced_config.R`, and then **rename** `advanced_config_O.R` to `advanced_config.R`.

For more **details**, check [our tutorial](https://esctrionsit.github.io/snphub_tutorial/)
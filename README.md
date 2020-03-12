# snphub

SnpHub is a Shiny-based server framework for retrieving, analyzing and visualizing the large genomic variations data in a lab.

For more **details**, check [our tutorial](https://esctrionsit.github.io/snphub_tutorial/)

For **Docker-encapsulated version**, check [here](https://github.com/esctrionsit/snphub4docker)

To get a **quick start**, see [here](https://esctrionsit.github.io/snphub_tutorial/content/Setup/overview.html)

[Here](http://wheat.cau.edu.cn/Wheat_SnpHub_Portal/) are our **live demos**.

![](SnpHub.jpg)

## Environment request

To run the SnpHub, make sure these softwares are **already** installed:
- samtools
- bcftools
- seqkit
- tabix

Alse, these R packages are also **needed**:
- ggplot2
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

Only two steps are needed to run on demo data set:
- Using command `git clone https://github.com/esctrionsit/snphub` to clone SnpHub to local.
- Using command like `R -e "shiny::runApp('./snphub', port=5000, host='0.0.0.0')"` to run on demo data set. *(Or copy the SnpHub into your `shiny-server app folder`)*

*Make sure that `samtools`, `bcftools`, `seqkit` and `tabix` are added to system PATH. Otherwise, `tool application paths` part in `setup.conf` is needed to change to fit.*

## Setup on own data set

There are two config files, which are `setup.conf` and `advanced_config.R`.

To setup on your own data set, you would need to:
- Edit the `setup.conf` file, make sure all the paths are correct.
- Delete the `advanced_config.R`, and rename the `advanced_config_O.R` as `advanced_config.R`. *(Otherwise SnpHub will still run on demo data set)*
- Use shell command `Rscript setup.R` to setup on your own data.
- Wait for SnpHub to finish.

For more **details**, check [our tutorial](https://esctrionsit.github.io/snphub_tutorial/)
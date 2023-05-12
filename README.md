# iprscanR 🔎
(pronounced _ipr-scanner_) <br>
An R package for accessing InterProScan5 API | for protein scanning & characterization

## Installation
`devtools::install_github("jravilab/iprscanr")` or <br>
`remotes::install_github("jravilab/iprscanr")`

## Usage
```
library(iprscanr)
submit_ipr(path2seq = here("inst/extdata/ex-in-CAA75348.1.faa"),
           outfolder = here("inst/extdata"),
           email = "jravilab.msu@gmail.com")
```

## Documentation
Learn more about the in-dev package [here](http://jravilab.github.io/iprscanr/vignettes/index.html).

How to a build a package? Material from [GLBIO 2023](//github.com/2023-glbio)

## License
[MIT License](https://github.com/JRaviLab/rinterpro/blob/main/LICENSE.md)

## Authors
Jacob Krol, Faisal Alquaddoomi, Janani Ravi

*General correspondence should be addressed to JR at janani.ravi@cuanschutz.edu.

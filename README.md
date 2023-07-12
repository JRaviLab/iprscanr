# iprscanR 🔎
(pronounced _ipr-scanner_) <br>
An R package for accessing InterProScan5 API | for protein scanning & characterization

## Installation
`devtools::install_github("jravilab/iprscanr")` or <br>
`remotes::install_github("jravilab/iprscanr")`

## Usage
```
library(iprscanr)
# use example data included with package
submit_ipr(path2seq = system.file(file.path("extdata", "ex-in-CAA75348.1.faa"), package = "iprscanr"),
           outfolder = system.file("extdata", package = "iprscanr"),
           email = "jravilab.msu@gmail.com",
           .applications = c("PfamA", "Phobius")
)

# print a character vector of all InterProScan applications
print(APPLICATIONS)
```

## Documentation
Learn more about the in-dev package [here](http://jravilab.github.io/iprscanr/vignettes/index.html).

How to a build a package? Material from [GLBIO 2023](//github.com/2023-glbio)

## License
[MIT License](https://github.com/JRaviLab/rinterpro/blob/main/LICENSE.md)

## Authors
Jacob Krol, Faisal Alquaddoomi, Janani Ravi

*General correspondence should be addressed to JR at janani.ravi@cuanschutz.edu.


![](https://github.com/JamesColl-NOAA/ePydRo/blob/dev/man/figures/3bcf273d-2154-42a0-8c19-6d6954d3bca2.png)

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://choosealicense.com/licenses/mit/)
<!-- badges: end -->

# Executive Summary

The [eHydro database](https://www.arcgis.com/apps/dashboards/4b8f2ba307684cf597617bf1b6d2f85d) represents a gold standard in river bathymetry measurements, composed of a living record of acoustic soundings and depth surveys across navigable waterways at regular intervals and approved by the US army Corp of Engineers.  However. utilizing that data still requires a great deal of manual manipulation in order to meet the range of applications this data can augment.  The ePydRo (eHydro with both Python and R) software presented here solves much of this by normalizng the ever-growing database of eHydro measurements into a form amenable for power users with fine grain access controls.

# Publicly facing catalog

See the transformed and staged versions of these data [here](https://spatial.dev.water.noaa.gov/data/surface/nws-ehydro/)

# Installation

As the badges above indicate, this package is Active and experimental.  Breaking changes will result in new data structures

``` r
# install.packages("devtools")                  # You only need this if this is your very first time opening RStudio
remotes::install_github("NOAA-OWP/ePydRo")    # If it asks for package updates: press 1
library(data.table)
ePydRo::marco()                           # A "hello world" test
```

> Note: ePydRo is an opinionated processing and staging utility for the already quite expansive [USACE Hydrographic Survey (eHydro) survey database](https://www.arcgis.com/apps/dashboards/4b8f2ba307684cf597617bf1b6d2f85d).  Users who are interested in directly working with that data may be better off starting from that authoritative source.  This repository stages that data in accessible and analysis ready forms for consumption across the hydrofabric and adjacent efforts.

## Tutorials

## Explainers

### What's ePydRo?
Short answer?:
This is a suite of tools and scripts written in both Python and R that are used to manipulate the eHydro database into a form that can be used for intimidate integration into the NOAA master modeling surface.

Slightly longer version:
Measuring the landscape is a complex science with several critical methodological inflection points that bake assumptions into the data that is read off the survey pole, recorded in the computer, pushed through the internal processing pipeline, and delivered to a client for hosting to the end users.  At that point those end users (us) come along and discover that the data is not in a form immediately integrable into our workflow.  For some applications, addressing this can simply be a means of harmonizing the naming schema between the two databases but in the case of terrain and elevations, the process of "harmonizing that naming schema" involves complex datum transformations and re-projections that place that measured point in the correct spot for the intended space.

This application transforms the eHydro data from their recorded datum to a normalized universal datum and reshapes the data structure to a more efficient form amenable to cloud access patterns and partial reads as an analysis ready value used to plug into the living master modeling surface constructed by the hydrofabric team.

### How does ePydRo data flow?

![_9743a3fb7f8bf395e9f8ac2282cc5490.png](https://github.com/JamesColl-NOAA/ePydRo/blob/dev/man/figures/_9743a3fb7f8bf395e9f8ac2282cc5490.png)

## HowTo's
The general workflow in cookbook form:
```{r}
ehydro_root <- "~/data/ePydRo"
ePydRo::collect_ehydro(ehydro_root)
                                                                                 
# Surveys around puerto-rico-virgin-islands
epydro_tracker <- sf::st_read(file.path(ehydro_root,"epydro_SurveyJob.fgb",fsep = .Platform$file.sep))
pr_surveys <- epydro_tracker[AOI::geocode("Puerto Rico" ,bbox = TRUE) |> sf::st_buffer(dist = 7000),]
                                                                                 
ePydRo::scrape_and_transform(ehydro_root,row.names(pr_surveys)[3])
ePydRo::to_tiff_and_stage(ehydro_root,row.names(pr_surveys)[3],"~/data/ePydRo/puerto-rico-virgin-islands")
```

## Credits and references

Credit to the packages used in the development, testing, and deployment of this package not exclusive of the ones found in the [DESCRIPTION](https://github.com/NOAA-OWP/ePydRo/blob/dev/DESCRIPTION) file.  We are incredibly grateful to the [USACE eHydro](https://www.sam.usace.army.mil/Missions/Spatial-Data-Branch/) team for the phenomenal job they do orchestrating the nations hydrographic survey collection and the open dissemination pattern they use.

### For questions

[Jim Coll](james.coll@noaa.gov) (FIM Developer), [Mike Johnson](mike.johnson@noaa.gov) (Geospatial Science and Technology lead), [Fernando Salas](fernando.salas@noaa.gov) (Director, OWP Geospatial Intelligence Division)

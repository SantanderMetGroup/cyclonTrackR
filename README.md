# What is cyclonTrackR?

cyclonTrackR is an R package build on [climate4R bundle](http://www.meteo.unican.es/climate4r) for the identification and traking of cyclones using different criteria based on two synoptic variables, the sea level pressure (slp) and the vorticiy of the lower troposphere (850 hPa).

This package contains two main functions, [getCyclonCenters.R](http://github.com/SantanderMetGroup/cyclonTrackR/blob/master/getCyclonCenters.R) and [getCyclonTrack.R](http://github.com/SantanderMetGroup/cyclonTrackR/blob/master/getCyclonTrack.R), the first to identify the cyclone centers and the second to define the track using the cyclone centers obtained.

The recommended installation procedure is to use the `install_github` command from the devtools R package (see the installation info in the wiki):

```r
devtools::install_github("cyclonTrackR")
```
---
Reference and further information: 

Iparragirre I. (2018) Climate change projections of explosive cyclogenesis events frequency in the Iberian Peninsula.  **Master's Thesis: Physics, Instrumentation and the Environment - Universidad de Cantabria***, http://meteo.unican.es/en/theses/ 

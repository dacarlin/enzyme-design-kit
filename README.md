# Enzyme design tool shed 

A collection of useful tools for students and practioners of enzyme design. A live version of this app may be fount at http://bagel.genomecenter.ucdavis.edu. And yes, we should probably change that URL. 

## Version History

- 0.1, original implementation of a Michaelis-Menten "fitter"
- 0.2, rewrite in Flask and deploy to a Genome Center VM with help from MCL and RF
- 0.3, implemented a fitter for thermal stability studies 

## Milestones

- 1.0, generalized UI for multiple tools and launch at new URL 
- 1.1, add oligo design (oligos for Kunkels)
- 1.2, add computational deep mutational scanning 
- 1.3, add comparative modeling 

## Available tools 

As of version 1.0, these are the available tools. 

### Determine <em>k</em><sub>cat</sub> and K<sub>M</sub> from kinetic data at various substrate concentrations   

This app allows you to fit you experimental data to the Michaelis-Menten equation to obtain the <em>k</em><sub>cat</sub> and K<sub>M</sub>. 

### Determine protein melting temperature from heat challenge data 

This app allows you to determine 

## How to add a new app (Contributing) 

Glad you asked! The simple answer is: create a new Python file, write a Flask view function, and import it into "application.py". 


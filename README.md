# Enzyme Design Kit: tools for enzyme design 

A collection of Python tools for enzyme design. Enzyme design approaches use a wide variety of tools, including alignment software from bioinformatics, 3D visualization software like PyMOL and UCSF Chimera, as well as a wide variety of synthetic DNA tools. This project is a collection of useful tools for enzyme design.

*Project Structure.* Briefly, the project has three parts: the documentation, the Python implementations of enzyme design and synthetic metagenomic algorithms, and the user interfaces (implemented as web applications using Flask. 

## Contributing 

1. If you spot a bug, you can file an issue (see "Issues" tab above)
1. Contributions are welcome! Within the Siegel group, ask Alex to become a contributor after you create a GitHub account. Otherwise, please fork this repo and submit a pull request, or email questions! 

## Available tools 

### Data analysis tools for determining Michalis-Menten constants from experimental data 

This is a web app for performing nonlinear regression and producing statistical plots. It is used for fitting experimental data to the Michaelis-Menten equation to characterize the functional parameters of enzyme mutants. 

### Data analysis tools for determining protein melting temperature from experimental data 

Similar to above, used for fitting thermal stability data to the logistic equation. 

### Mutagenic oligo design for Kunkel and nicking mutagenesis 

Use a coding nucleotide sequence to generate oligos for Kunkel mutagenesis and nicking mutagenesis 


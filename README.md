# Bonds

![License][license]

This repository contains the [LaTeX][latex] sources for my manuscript intitled "Efficient evaluation of interatomic distances in large atomic scale models". 
It is prepared for publication in the [Living Journal of Computational Molecular Science][livecom]. 

It also contains the latest PDF version: [Efficient evaluation of interatomic distances in large atomic scale models][bonds]

## Build instructions

To build the PDF file of the manuscript:

  1. Clone this repository

  ```
  git clone https://github.com/Slookeur/Bonds
  ```

  2. Change directory to `Bonds`

  ```
  cd Bonds
  ```

  3. Use the associated script `viewt`

  ```
  ./viewt bonds-live
  ```

> [!IMPORTANT]
> It is required to have working TeX and LaTeX distributions to build the mansucript, either:
>  - [TeX Live][texlive] distribution for GNU/Linux
>  - [MacTeX][mactex] distribution for MacOS
>  - [MiKTeX][miktex] distribution for Windows
> It also requires all needed dependencies, in particular the [`livecoms`][livecoms_class] document class.


> [!WARNING]
> Issues in version 1.0:
> - The ̀\affil[1]{` instruction in `bonds-lievcom.tex` line **136** is always producing an error, and it is required to comment this line to build the PDF.
> - On page 6 of 25, now way to print all figures in the same column and then continue the discussion, for some reason(s) the third figure ("Corner of the pixel grid") is always inserted on second column.

This manuscript was written by [Dr. Sébastien Le Roux][slr], research engineer for the [CNRS][cnrs]

<p align="center">
  <a href="https://www.cnrs.fr/"><img width="100" src="https://www.cnrs.fr/themes/custom/cnrs/logo.svg" alt="CNRS logo" align="center"></a>
</p>

[Dr. Sébastien Le Roux][slr] works at the Institut de Physique et Chimie des Matériaux de Strasbourg [IPCMS][ipcms]

<p align="center">
  <a href="https://www.ipcms.fr/"><img width="100" src="https://www.ipcms.fr/uploads/2020/09/cropped-dessin_logo_IPCMS_couleur_vectoriel_r%C3%A9%C3%A9quilibr%C3%A9-2.png" alt="IPCMS logo" align="center"></a>
</p>

[license]:https://img.shields.io/badge/License-CC_BY_NC_4.0-blue
[livecom]:https://livecomsjournal.org/index.php/livecoms/index
[LaTeX]:https://www.latex-project.org/
[TeX]:https://www.tug.org/
[texlive]:http://www.tug.org/texlive/
[mactex]:http://www.tug.org/mactex/
[miktex]:http://miktex.org/
[bonds]:bonds-livecoms.pdf
[slr]:https://www.ipcms.fr/sebastien-le-roux/
[cnrs]:https://www.cnrs.fr/
[ipcms]:https://www.ipcms.fr/
[github]:https://github.com/
[livecoms_class]:https://github.com/livecomsjournal/article_templates

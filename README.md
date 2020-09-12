<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/MHH-RCUG/scrnaseq_app">
    <img src="www/MHH.png" alt="Logo">
  </a>

  <h2 align="center">scrnaseq_app</h2>

  <p align="center">
    Shiny app for visualisation of scRNASeq data. In the UI you can select genes and the app will generate several plots (Feature, Ridge, Violin and Dot-Plot).
    This project is based around <a href="https://github.com/ktrns/scrnaseq">scrnaseq</a>. This workflow and was developed by <a href="https://github.com/ktrns">Katrin Sameith</a> and <a href="https://github.com/andpet0101">Andreas Petzold</a> at the <a href="https://genomecenter.tu-dresden.de/about-us">Dresden-concept Genome Center (Dresden, Germany)</a>. Through collaboration with the <a href="https://www.mhh.de/genomics">Research Core Unit Genomics (Hannover, Germany)</a> the workflow has grown substantially and has been standardised.
    <br />
    <a href="https://github.com/MHH-RCUG/scrnaseq_app"><strong>Explore the docs »</strong></a>
    <br />
    <br />
     :bug:
    <a href="https://github.com/MHH-RCUG/scrnaseq_app/issues">Report Bug</a>
    ·  :question:
    <a href="https://github.com/MHH-RCUG/scrnaseq_app/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
  * [Built With](#built-with)
* [Usage](#usage)
  * [Upload .rds file](#upload-rds-file)
  * [Select genes](#select-genes)
    * [Select genes through UI](select-genes-through-ui)
    * [Select genes through an excel file with header](select-genes-through-an-excel-file-with-header)
    * [Select genes through an excel file without header](select-genes-through-an-excel-file-without-header)
* [Getting Started](#getting-started)
* [Contributing](#contributing)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)



<!-- ABOUT THE PROJECT -->
## About The Project

![screenshot](https://github.com/MHH-RCUG/scrnaseq_app/blob/master/www/screenshot.PNG)
The app reads a .rds file which contains the seurat-object. Therefore every data is available and the user can generate plots for every gene they select through the select list or an excel file (.xlsx). If the user has the selection of genes that are interesting, it is possible to download those plots in .pdf and .png format.

### Built With

* [R](https://www.r-project.org/)
* [Shiny from RStudio](https://shiny.rstudio.com/)
* [Seurat](https://satijalab.org/seurat/)



<!-- USAGE EXAMPLES -->
## Usage

To get the app running you need to follow these simple steps.

### Upload .rds file
First you need to upload a .rds file containing your Seurat object. If you don't have this file, you can use this [example](https://owncloud.gwdg.de/index.php/s/rRawkhIOVe1T5qi).

### Select genes

#### Select genes through UI
Once the .rds file is uploaded and processed, you can click on the field below `Genes` and select your genes by clicking or searching. To delete genes you have to click on the gene and hit the backspace key.

#### Select genes through an excel file with header
You can select genes by uploading an .xlsx file ([example](https://owncloud.gwdg.de/index.php/s/ZwY0iVPji6uBVKO)).  
The spreadsheet has to look something like this:

| Ensembl_IDs | Genes       |
|-----------|-------------|
| First ID  | First Gene  |
| Second ID | Second Gene |

The app will only read the first column where it expects Ensembl IDs, whereas the column `Genes` is optional.

#### Select genes through an excel file without header
If your .xlsx file contains no header please uncheck the checkbox. Otherwise the first row/ EnsemblID will not show up in the selection.

### Aspect Ratio
In order to change the aspect ratio of the plots you can change the pixels of the X and Y-Axis. By clicking on `Default settings for axes` the default settings will be restored: `X = 1024px`and `Y = 576px`

### Download
Once the plots are generated you can download them. The download is an archive which contains a **PDF and PNG** version of the selected plots.
If you want to rename your archive before downloading, you can do so in the given field.



<!-- GETTING STARTED -->
## Getting Started

### Installation

#### Singularity
If you are in the network of the MHH, you can use this [link](http://172.24.148.210:3838/)

#### Local
This is an example of how you may give instructions on setting up your project locally. To get a local copy up and running follow these simple example steps.

1. Head to [Releases](https://github.com/MHH-RCUG/scrnaseq_app/releases) and download the latest release
2. Unzip the project und run scrnaseq_app.R with [RStudio](https://rstudio.com/)
3. Install needed libraries


<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request



<!-- LICENSE 
## License

Distributed under the MIT License. See `LICENSE` for more information.
-->


<!-- CONTACT -->
## Contact

Project Link: [https://github.com/MHH-RCUG/scrnaseq_app](https://github.com/MHH-RCUG/scrnaseq_app)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

* [scrnaseq](https://github.com/ktrns/scrnaseq)
* [README-Template](https://github.com/othneildrew/Best-README-Template)  


[![Contributors](https://img.shields.io/github/contributors/MHH-RCUG/scrnaseq_app)](https://github.com/MHH-RCUG/scrnaseq_app/graphs/contributors)
[![GitHub issues](https://img.shields.io/github/issues/MHH-RCUG/scrnaseq_app)](https://github.com/MHH-RCUG/scrnaseq_app/issues)
[![GitHub forks](https://img.shields.io/github/forks/MHH-RCUG/scrnaseq_app)](https://github.com/MHH-RCUG/scrnaseq_app/network)
[![GitHub stars](https://img.shields.io/github/stars/MHH-RCUG/scrnaseq_app)](https://github.com/MHH-RCUG/scrnaseq_app/stargazers)

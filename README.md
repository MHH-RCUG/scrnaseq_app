<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://www.mhh.de/">
    <img src="shiny/www/Logo_engl_schwarz.png" height="25%" width="25%" alt="MHH">
  </a>
  
  <a href="https://www.mhh.de/genomics">
    <img src="shiny/www/RCUG_Logo.png" height="15%" width="15%" alt="Logo">
  </a>
</p>

# scrnaseq_app

This is an R/Shiny app for visualisation of single cell RNA-seq (scRNASeq) data. Based on a provided Seurat object, one or more genes can be selected through a user interface, and several plots are generated, including feature plots, ridge plots, violin plots and dot plots. This project is based on [scrnaseq](https://github.com/ktrns/scrnaseq).
The scrnaseq workflow was and is being developed by [Katrin Sameith](https://github.com/ktrns) and [Andreas Petzold](https://github.com/andpet0101) at the [DRESDEN-concept Genome Center (Dresden, Germany)](https://genomecenter.tu-dresden.de/about-us), in collaboration with the [Research Core Unit Genomics (Hannover, Germany)](https://www.mhh.de/genomics).

[:bug: Report Bug](https://github.com/MHH-RCUG/scrnaseq_app/issues) 
·
[:question: Request Feature](https://github.com/MHH-RCUG/scrnaseq_app/issues)

<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
  * [Built With](#built-with)
* [Usage](#usage)
  * [Upload .rds file](#upload-rds-file)
  * [Select genes](#select-genes)
* [Getting Started](#getting-started)
* [Contributing](#contributing)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)

<!-- ABOUT THE PROJECT -->
## About The Project

![screenshot](https://github.com/MHH-RCUG/scrnaseq_app/blob/master/shiny/www/scrnaseq_app_upload.png)
This app is dependent on data generated by Seurat, which must be run first. The app reads an .rds file which contains the Seurat object and all associated data. 
The user can then select genes through a list or an Excel file (.xlsx). The expression of these genes is then displayed in a) UMAP plots b) Ridge plots c) Violin plots d) Dot plots or e) a global Heatmap. The generated plots can then be downloaded in .pdf and .png format.

### Built With

* [R](https://www.r-project.org/)
* [Shiny from RStudio](https://shiny.rstudio.com/)
* [Seurat](https://satijalab.org/seurat/)

<!-- USAGE EXAMPLES -->
## Usage

To get the app running, you need to follow these simple steps.

### Upload .rds file

First, you need to upload an .rds file containing your Seurat object. If you do not have this file yet but would still like to try the app, you can use this [example](https://owncloud.gwdg.de/index.php/s/rRawkhIOVe1T5qi).

### Select genes

#### Select genes through the user interface

Once the .rds file is uploaded and processed, you can click on the field below `Genes` and select your genes by clicking or searching. To delete genes again, you have to click on the gene and hit the backspace key.

#### Select genes through an Excel file with header

You can select genes by uploading an .xlsx file ([example](https://owncloud.gwdg.de/index.php/s/ZwY0iVPji6uBVKO)).  
The spreadsheet has to look as follows:

| Ensembl_IDs | Gene_Names |
| :--- | :--- |
| First ID  | First Name |
| Second ID | Second Name |

The app will only read the first column including the Ensembl IDs. The second column `Genes` is optional.

#### Select genes through an Excel file without header

If your .xlsx file contains no header, please uncheck the checkbox. Otherwise, the first Ensembl ID will not show up in the selection.

### Aspect Ratio

In order to change the aspect ratio of the plots you can change the pixels of the x and y-axis. By clicking on `Default settings for axes` the default settings will be restored: `X = 1024px`and `Y = 576px`

### Download

Once the plots are generated, you can download them. The download is an archive which contains a **PDF and PNG** version of the selected plots.
If you want to rename your archive before downloading, you can do so in the given field.

### Troubleshooting

Should the app hang or stop responding, you can reset it by reloading the data (RDS file) into the app.

<!-- GETTING STARTED -->
## Getting Started

### Installation

#### Singularity

If you are within the network of the MHH, use this [link](http://172.24.148.210:3838/) if the singularity container is running.

#### Local

To run a copy of this app locally, follow these steps:

1. Go to [Releases](https://github.com/MHH-RCUG/scrnaseq_app/releases) and download the latest release
2. Unzip the project und run scrnaseq_app.R with [RStudio](https://rstudio.com/)
3. Install required libraries

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

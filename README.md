



<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/MHH-RCUG/scrnaseq_app">
    <img src="www/MHH.png" alt="Logo">
  </a>

  <h2 align="center">scrnaseq_app</h2>

  <p align="center">
    A shiny app to generate and download plots
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
Once the .rds file is uploaded and processed, you can click on the field below "Genes" and select your genes by clicking or searching. To delete genes you have to click on the gene and hit the backspace key.

#### Select genes through an excel file with header
You can select genes by uploading an .xlsx file ([example](https://owncloud.gwdg.de/index.php/s/ZwY0iVPji6uBVKO)).  
The spreadsheet has to look something like this:

| EnemblIDs | Genes       |
|-----------|-------------|
| First ID  | First Gene  |
| Second ID | Second Gene |

The app will only read the first column where it expects EnsemblIDs, whereas rhe column "Genes" is optional.

#### Select genes through an excel file without header
If your .xlsx file contains no header please uncheck the checkbox. Otherwise the first row/ EnsemblID will not show up in the selection.

<!-- GETTING STARTED -->
## Getting Started

### Installation

#### Singularity
If you are in the network of the MHH, you can use this [link](http://172.24.148.210:3838/)

#### Local
This is an example of how you may give instructions on setting up your project locally. To get a local copy up and running follow these simple example steps.

1. Clone the repo
```sh
git clone https://github.com/MHH-RCUG/scrnaseq_app
```

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

# Help

**Please visit the GitHub-Page for any further information!**

[<img src="www/GitHub_logo.png" width="30%" height="30%">](https://github.com/MHH-RCUG/scrnaseq_app/blob/master/README.md)

## 1. Upload .rds file

First of all, you have to upload the .rds file which contains the seurat object. Once the upload has completed, the app will process the data. Please wait until everything has finished.

## 2. Select genes

After the data from the .rds file have been processed, you are able to select genes through two methods:

### 2.1 Select Genes Manually

If you select the genes manually, you will get a list of all genes to select from.
In order to find certain genes, you can simply just type in the name and/ or ensemblID.

### 2.2 Select through Gene List

You can chose a predefined list and all associated genes will be loaded, after clicking  "Select List".

### 2.3 Select through Excel File

Create an excel file (.xlsx) which contains the ensemblIDs in the first column. It is optional to include column headings.  
Template:

[<img src="www/excel_file_template.png" width="20%" height="20%">]()

Upload this file to the app and wait for the data to be processed.
Once processed, the selected genes will be displayed and you can correct your selection.

## 3. Download

By default the download will contain every plot. If you only want certain plots, you can uncheck the ones you don't want. You can rename the archive and hit download.
This can take a moment, because the .png and .pdf files have to be created.

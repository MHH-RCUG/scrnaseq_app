# singularity_scrna_app
A test environment for the R Shiny scrnaseq app by Marius Rueve

Singularity app by Colin Davenport



## How to build using Singularity: build process command:

```
### takes ~10 minutes (uses up to 24 cores on a multicore machine!)
### build main version
sudo singularity build ../scrnaseq_app.simg singularity_app_recipe.txt
### build dev version
sudo singularity build ../scrnaseq_app.simg singularity_app_recipe_dev.txt


```

## Run app - exits without msg
```
singularity run ../scrnaseq_app.simg

singularity exec ../scrnaseq_app.simg ls  /data/scrnaseq/

singularity shell ../scrnaseq_app.simg 
```

### Actual app and repo to run:
```
https://github.com/MHH-RCUG/scrnaseq_app
```

### Test app in container. Currently missing pbmc_2020-06-08.rds', probable reason 'No such file or directory'
```
singularity exec ../scrnaseq_app.simg     R -e "options('shiny.port'=3838,shiny.host='0.0.0.0'); shiny::runApp('/data/scrnaseq_app/scrnaseq_app.R')"
```


```
Tested/built on:

singularity version
-3.5.3
-Ubuntu 1604
```

### Helpful links - singularity setup and SLURM

https://njstem.wordpress.com/2018/08/02/r-script-seurat-with-a-singularity-container-using-slurm/


### Dependencies : downloaded in script
```
git clone https://github.com/MHH-RCUG/scrnaseq_app
wget https://owncloud.gwdg.de/index.php/s/rRawkhIOVe1T5qi/download
```

# singularity_scrna_app
An environment for the R Shiny scrnaseq app by Marius Rueve

Singularity app by Colin Davenport



## How to build using Singularity: build process command:

```
### takes ~10 minutes (uses up to 24 cores on a multicore machine!)
### build main version
sudo singularity build ../../scrnaseq_app.simg singularity_app_recipe.txt
### build dev version
sudo singularity build ../../scrnaseq_app.simg singularity_app_recipe_dev.txt


```

## To check if the build process did well - exits without msg
```
singularity run ../../scrnaseq_app.simg
singularity exec ../../scrnaseq_app.simg ls  /data/scrnaseq/
singularity shell ../../scrnaseq_app.simg 
```

### Actual app and repo to run:
```
https://github.com/MHH-RCUG/scrnaseq_app
```

### Run app in container and make available for any browser that has contact to the server (e.g. inside the campus network).
### Currently missing pbmc_2020-06-08.rds', probable reason 'No such file or directory'
```
# kill previous container, if exists
ps -aux | grep R
kill <PID> 

# Start container
singularity exec  ../../scrnaseq_app.simg     R -e "options('shiny.port'=3838,shiny.host='0.0.0.0'); shiny::runApp('/data/scrnaseq_app/runAppp.R')"
```


```
Tested/built on:

singularity version
-3.5.3
-Ubuntu 1604
-Ubuntu 2004
```

### Helpful links - singularity setup and SLURM

https://njstem.wordpress.com/2018/08/02/r-script-seurat-with-a-singularity-container-using-slurm/


### Writable containers
-Just for sandbox directories for testing
-Instead, write images to directories which are bound into the container, eg /home/$USER, or /tmp

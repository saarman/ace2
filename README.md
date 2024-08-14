# ace2
Culex pipiens species complex PCR assay results, Ace2 gene, samples from VecTech collaboration

# Logging into the CHPC server:  
One way to get to your CHPC user interface and have access to your home directory and group storage directory:
https://ondemand.chpc.utah.edu/pun/sys/dashboard/files/fs/uufs/chpc.utah.edu/common/home/u6036559  
...replacing uXXXXXXX with your U of U user ID.  

# Connecting to RStudio server on CHPC:
https://ondemand.chpc.utah.edu/pun/sys/dashboard/batch_connect/sys/rstudio_server_app/session_contexts/new       
**R version:** R 4.0.3 Geospatial packages      
**Cluster:** notchpeak     
**Account:Partition** saarman-np:saarman-shared-np (very important!!! this allows multiple jobs to run simultaneously, if you don't use this partition you will block entry for everyone else!)     
**Number of cores:** 4 (there are 32 total)  
**Number of hours:** 100 (336 is the max, definitely put more than 24)  
**Memory per job:** 128G (1000 GB total is the max across all of the 32 cores + all users, but good to stay well below half of the limit) 

**You can access a running session by going to https://ondemand.chpc.utah.edu/pun/sys/dashboard/batch_connect/sessions 

#  Slurm
**Slurm** to connect and allocate a job on the correct node. For example to initiate an “interactive” job, run commands such as:

```
salloc --time=1:00:00 --ntasks 1 --mem=100G --account=saarman-np --partition=saarman-shared-np
```

# Secure Shell (ssh) from Terminal (in bash):
**ssh** directions: https://www.chpc.utah.edu/documentation/software/ssh.php
ssh <uNID>@kingspeak.chpc.utah.edu #syntax  
```
ssh u6036559@notchpeak.chpc.utah.edu #example   
```

Allocate a slurm interactive job
```
salloc --time=1:00:00 --ntasks 1 --mem=100G --account=saarman-np --partition=saarman-shared-np
```
   
# Modules can be used with FastX3 or Secure Shell
Modules are used to access installed software:  https://www.chpc.utah.edu/documentation/software/modules.php   

```
module spider plink
module load plink/2.0
```      

# Connect github with Rstudio 
1) create a personal access token  
2) Two ways to do this  
** In browser: https://docs.github.com/en/enterprise-server@3.4/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token
** In command line: https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token#using-a-token-on-the-command-line  
3) enable git in Rstudio if not already available, more or less following these steps https://hansenjohnson.org/post/sync-github-repository-with-existing-r-project/
While doing this process, you will be asked to start a new RStudio project, for this it is best to use a temp Rproject in a "new directory"... once you have activated Github this folder can be deleted.

4) Once Git is enabled, create a new Project through the File --> New Project --> Version Control --> Git

You will be asked for your username and token if you are not already signed into github on the Rstudio server

5) Then use the Git menu in your RStudio Server top right window to change to the remote repository, commit, pull, push, etc.

# Remember to update Permissions
For example:
```
chmod -R g+w /uufs/chpc.utah.edu/common/home/saarman-group1/raster_files  
```





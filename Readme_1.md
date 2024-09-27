# HPC Summer School 2024: Foundation in Computational Biomolecular and Biosystem Research

This guide demonstrates how to analyze molecular dynamics data using the MAHAMERU BRIN HPC. In this workshop, all necessary software is pre-installed as modules, including ORCA (for quantum chemistry calculations), while Yasara Structure (licensed) is already set up on local computers. Therefore, no additional environment setup is required. To run the examples, simply click the top-right corner of each code box to copy the provided snippets, and paste them directly into your terminal.

## Login to MAHAMERU
Open your terminal (Linux, MacOS) or PowerShell (Windows) and enter the following commands
```
ssh username@login2.hpc.brin.go.id
```
You will be connected to the `trembesi02` node, which serves as a **login** node.

## Creating a New Working Directory
To keep your files organized, create a directory named `workshop_yasara` under your home directory by typing the following command:
```
mkdir workshop_yasara
```
To navigate to the newly created `workshop_yasara` directory, type:
```
cd workshop_yasara
```
Press Enter.

## Sending Your Preparation Job from Your Laptop
To send files from your laptop to the remote directory using SCP (Secure Copy Protocol), type the following command:
```
scp /path/to/your/file username@remote_host:/home/username/workshop_yasara/
```
Here:

`/path/to/your/file` is the location of the file on your laptop.
`username` is your remote machine's username.
`remote_host` is the IP address or hostname of the remote machine.

For example, if your file is located at `/Documents/Yasara/testfile.txt` and your remote server’s username is user123, and the host is 192.168.1.5, the command would look like:
```
scp /Documents/Yasara/testfile.txt user123@192.168.1.5:/home/user123/workshop_yasara/
```
Alternatively, if you want to copy multiple files:
```
scp -rp * fari025@login2.hpc.brin.go.id:/mgpfs/home/fari025/workshop_yasara
```
After running the command, Once authenticated, the files will be transferred to the `workshop_yasara` directory.

## Verifying File Upload
Back on the terminal connected to the HPC, verify that the files have been uploaded successfully by listing the contents of the current directory, type:
```
ls
```
You should see a list of your uploaded files.

![Screenshot 2024-09-02 at xx](XX.jpg)

## Submitting the Job
To submit a job for execution on the HPC (High-Performance Computing) cluster, use the sbatch command followed by the script name. In this case, the script is named `Yasara_MD.sh`. This script contains the instructions for the HPC to run the molecular dynamics simulation.
```
sbatch Yasara_MD.sh
```
When you run this command, the job is submitted to the scheduler (SLURM), which queues your task and assigns the necessary computational resources (CPUs, memory, etc.) to execute it.

## Checking Job Status
After submitting your job, you can monitor its progress using the `squeue` command. This command allows you to view the status of jobs currently in the queue.
```
squeue -u whoami
```
Explanation:

`squeue` is used to query the scheduler for jobs that are either waiting or running.
`-u` specifies the user, and `whoami` automatically fetches your username, so it shows the jobs associated with your account.

## Interpreting the Job Status
Once you run the squeue command, you'll see a list of jobs along with their details. The important column to note is the "ST" column, which represents the status of your job. Here are common job statuses:
**R**: The job is Running.
**PD**: The job is Pending, meaning it's in the queue and waiting for resources.
**CG**: The job is Completing, meaning it is finishing up.

If everything is working correctly, the status (**ST**) should display as '**R**', indicating that your job is actively running on the cluster.

## Analysing MD results
After the simulation is completed, you can analyze the molecular dynamics results. To do so, modify your Yasara_MD.sh script to use md_analyze.mcr as the input for the analysis.

Step-by-Step Guide to Edit Using vi
1. Open the File In the terminal, type:
```
vi Yasara_MD.sh
```
2. Enter Insert Mode Press `i` to start editing the file.
```
#!/bin/bash
  
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32GB
#SBATCH --partition=medium-large

#SBATCH --output=md_runmembrane.out
#SBATCH --error=md_runmembrane.err

yasara_exec=/mgpfs/home/lala002/apps/yasara/yasara

FILE_INPUT=md_runmembrane.mcr

nice -n 20 ${yasara_exec} -txt ${FILE_INPUT}

```

Edit the File Add the line:
`md_runmembrane`
into
`md_analyze`

4. Save and Exit
Press `Esc` to exit editing mode.
Type `:wq` and press Enter to save and quit vi.

## Closing the Terminal
Even though the job may take some time to complete, you can safely close the terminal. The job will continue running in the background on the HPC system until it's finished. You can always reconnect later and check the status using the `squeue` command again.

By using `sbatch`, the job is handled by the scheduler independently of your terminal session, so there's no need to keep the terminal open while the job runs.



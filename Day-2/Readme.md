![Screenshot 2024-09-02 at xx](https://github.com/lala002-brin/BRIN_ComChem_workshop/blob/main/attachment/header.jpg) 

# HPC Summer School 2024: Foundation in Computational Biomolecular and Biosystem Research


## Login to MAHAMERU
To start, ensure you have your private key file (`keybrinapctp0i`) handy. This file is essential for authenticating your login to the HPC system.

Steps:
1. Open Terminal:
    - Linux/MacOS: Press `Ctrl + Alt + T` to open your terminal or search for "Terminal" in your applications.
    - Windows: Open PowerShell by typing "PowerShell" in the Start Menu.
      
2. Locate Your Private Key:
   Find the directory where you saved your private key file (`keybrinapctp0i`).

   For example:

    - If the key file is saved in the `Downloads` directory on Linux/MacOS, the path might look like
    
      `/home/your-username/Downloads/keybrinapctp0i`.
   
    - On Windows, it might be something like
    
      `C:\Users\your-username\Downloads\keybrinapctp0i`.

    **For Linux/Mac users**: 
    To ensure your private key has the correct permissions, you need to change its permission level to `600` to restrict access. This can be        done with the `chmod` command, which allows only the owner to read and write the key file. Use the following command:
    ```
    chmod 600 ~/Downloads/keybrinapctp03
    ``` 
5. Connect to the HPC: Run the `SSH` command with the correct path to your private key:
   
    ```
    ssh wsbrinapctp0i@login2.hpc.brin.go.id -i /full/path/to/keybrinapctp0i
    ``` 
    Example
    If the `id_ws000i` file is stored in `/home/username/Downloads/keybrinapctp0i` on a Linux/MacOS system, the command will look like:
    ```
    ssh wsbrinapctp0i@login2.hpc.brin.go.id -i /home/username/Downloads/keybrinapctp0i
    ```
    For a Windows system, if the file is in `C:\Users\username\Downloads\keybrinapctp0i`, the command will be:
    ```
    ssh wsbrinapctp0i@login2.hpc.brin.go.id -i C:\Users\username\Downloads\keybrinapctp0i
    ```

6. Authenticate:
    If this is your first time connecting to `login2.hpc.brin.go.id`, you'll be asked to confirm the host identity. Type `yes` to continue.
    If your private key (`id_ws000i`) is encrypted, you'll be prompted to enter the passphrase associated with the key.
    type:

    `brinapctp` as the passphrase.

    The system will not display characters as you type, but simply press `Enter` after typing.

    After completing these steps, you should be successfully connected to the `trembesi02` login node of the HPC BRIN system.

![Screenshot 2024-09-02 at xx](LoginHPC.jpeg)


## Creating a New Working Directory
Now, let's create a dedicated directory for your workshop files to stay organized:
1. **Create the Directory**: To keep your files organized, create a directory named `workshop_yasara` under your home directory by typing the following command:
     ```
    mkdir workshop_yasara
    ```

2. **Navigate to the Directory**:
    ```
    cd workshop_yasara
    ```

3. **Confirm Your Location**: To ensure you are in the correct directory, use the `pwd` (print working directory) command, which will display the full path of your current directory:
    ```
    pwd
    ```
    This will show the current location in the directory structure, helping you confirm that you are now inside the `workshop_yasara` folder.

    For example, if everything is correct, you will see output similar to this:

    `/mgpfs/home/wsbrinapctp0i/workshop_yasara`

 Now, you're ready to start working within your new workshop_yasara folder.
 
## Sending Files to the HPC
To upload files from your laptop to the HPC system, use SCP (Secure Copy Protocol). This example shows how to upload your private key file:

1. **Send a File**: To send files from your laptop to the remote directory using SCP, type the following command (in your laptop):
    ```
    scp -i ~/Downloads/id_ws000i file_input wsbrinapctp0i@login2.hpc.brin.go.id:/mgpfs/home/wsbrinapctp0i/workshop_yasara/

    ```
    Here:

    - `~/Downloads/id_ws000i` is the location of your private key file in the Downloads directory.
    - `wsbrinapctp0i` is your remote username.
    - `login2.hpc.brin.go.id` is the remote host.
    - `/home/wsbrinapctp0i/workshop_yasara/` is the destination folder on the remote server.

or, alternatively:
    
2. **Upload a Folder**: if you want to copy a folder, make sure that the "workshop_yasara" directory exists on the HPC. If it hasn’t been created yet, you can create the directory while copying the folder using this command:

    ```
    scp -rp workshop_yasara wsbrinapctp0i@login2.hpc.brin.go.id:/mgpfs/home/wsbrinapctp0i/
    ```

Back on the terminal connected to the HPC, verify that the files have been uploaded successfully by listing the contents of the current directory, type:
```
ls
```
You should see a list of your uploaded files.

![Screenshot 2024-09-02 at xx](YASARA.jpg)


## Submitting the Job
To submit a job for execution on the HPC cluster, use the sbatch command followed by the script name. In this example, the script is called Yasara_MD.sh, which contains the necessary instructions for running the molecular dynamics simulation.

#### Yasara_MD.sh Script Contents
```
#!/bin/bash

#SBATCH --nodes=1              # Request 1 node
#SBATCH --ntasks=4             # Request 4 CPU cores
#SBATCH --mem=16GB             # Request 16 GB of memory
#SBATCH --partition=medium-large  # Specify the partition to use

#SBATCH --output=md_analyze.out  # Save standard output to md_analyze.out
#SBATCH --error=md_analyze.err   # Save error messages to md_analyze.err

FILE_INPUT=md_analyze.mcr  # Define the input file for the simulation

yasara -txt ${FILE_INPUT}  # Run Yasara using the specified input file
```

Explanation of the Script:
- `#!/bin/bash`: This line specifies the script will be executed in a Bash shell.
- `#SBATCH --nodes=1`: Allocates 1 compute node for the job.
- `#SBATCH --ntasks=4`: Requests 4 tasks (or CPU cores) for parallel execution.
- `#SBATCH --mem=16GB`: Allocates 16 GB of memory to the job.
- `#SBATCH --partition=medium-large`: Specifies the partition (queue) to submit the job to, in this case, a medium-large queue suitable for moderately resource-intensive tasks.
- `#SBATCH --output=md_analyze.out`: Redirects the standard output to the md_analyze.out file.
- `#SBATCH --error=md_analyze.err`: Redirects error messages to the md_analyze.err file for debugging.
- `FILE_INPUT=md_analyze.mcr`: Defines the input file for the molecular dynamics simulation.
- `yasara -txt ${FILE_INPUT}`: Runs the Yasara program in text mode using the input file md_analyze.mcr.


Then, to submit this job to the HPC scheduler (SLURM), run the following command in your terminal:
```
sbatch Yasara_MD.sh
```
#### What Happens When You Submit:
- The `sbatch` command submits the job to SLURM, the job scheduling system used by the HPC.
- The scheduler then queues your job and assigns it the required resources (CPU, memory, etc.) once available.
- The simulation will run, with its output saved in `md_analyze.out` and any errors captured in `md_analyze.err`.
- This process ensures that your molecular dynamics simulation runs efficiently on the HPC cluster.


## Checking Job Status
After submitting your job, you can monitor its progress using the `squeue` command. This command allows you to view the status of jobs currently in the queue.
```
squeue -u wsbrinapctp0i
```
Explanation:

- `squeue` is used to query the scheduler for jobs that are either waiting or running.
- `-u` specifies the user, and 
- `wsbrinapctp0i` automatically fetches your username, so it shows the jobs associated with your account.

## Interpreting the Job Status
Once you run the squeue command, you'll see a list of jobs along with their details. The important column to note is the "ST" column, which represents the status of your job. 
Here are common job statuses:
- `R`: The job is Running.
- `PD`: The job is Pending, meaning it's in the queue and waiting for resources.
- `CG`: The job is Completing, meaning it is finishing up.

If everything is working correctly, the status (`ST`) should display as `R`, indicating that your job is actively running on the cluster.

## Monitoring Output in Real-Time
To check the real-time progress of your job, you can use the `tail -f` command to monitor your output file. For example, if your job is running a molecular dynamics simulation and writing output to a file called `md_runmembrane.mcr`, you can view the live updates with this command:
```
tail -f md_runmembrane.mcr
```
Explanation:
- `tail -f` continuously displays the last few lines of the file, and it updates in real-time as new lines are added.
- `md_runmembrane.mcr` is the name of the output file where the progress of your job is being logged.

This allows you to track the progress of your job without waiting for it to finish.

## Analysing MD results
After the simulation is completed, you can analyze the molecular dynamics results. To do so, use the `Yasara_analysis.sh` script, which is specifically designed to process and analyze the output from your simulation.

#### Yasara_analysis.sh Script Contents
```
#!/bin/bash

#SBATCH --nodes=1               # Request 1 node
#SBATCH --ntasks=4              # Request 4 CPU cores
#SBATCH --mem=8GB               # Request 8 GB of memory
#SBATCH --partition=medium       # Use a medium partition for analysis

#SBATCH --output=analysis.out    # Save standard output to analysis.out
#SBATCH --error=analysis.err     # Save error messages to analysis.err

FILE_ANALYSIS=md_analyze.mcr     # Define the input file for analysis

yasara -txt ${FILE_ANALYSIS}     # Run Yasara in text mode for analysis

```

To execute the `Yasara_analysis.sh` script and analyze the results, use the following command:
```
sbatch Yasara_analysis.sh
```
This will submit the analysis job to the HPC, and the system will allocate the necessary resources. The results of the analysis will be saved in `analysis.out`, and any errors will be logged in `analysis.err`.

To confirm that the analysis has completed successfully and view the results, you can list the files in your directory:
```
ls
```

You should see a list of files:

![Screenshot 2024-09-02 at xx](analysis.png)


## Closing the Terminal
Even though the job may take some time to complete, you can safely close the terminal. The job will continue running in the background on the HPC system until it's finished. You can always reconnect later and check the status using the `squeue` command again.

By using `sbatch`, the job is handled by the scheduler independently of your terminal session, so there's no need to keep the terminal open while the job runs.

To download the result files from the HPC to your laptop, you can use the `scp` command. Here's how to do it:

**Download the Result Files to Your Laptop**: 

Use the following command on your laptop to download the 4hjo_001_report.html files from the HPC to your local machine:
```
scp wsbrinapctp0i@login2.hpc.brin.go.id:/mgpfs/home/wsbrinapctp0i/workshop_yasara/4hjo_001_report.html ~/Downloads/.
```

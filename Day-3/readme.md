![Screenshot 2024-09-02 at xx](https://github.com/lala002-brin/BRIN_ComChem_workshop/blob/main/attachment/header.jpg) 

# HPC Summer School 2024: Foundation in Computational Biomolecular and Biosystem Research

### Software:
- Avogadro: Avogadro adalah alat penting untuk pengeditan dan visualisasi molekul. Pastikan sudah terinstal di laptop. Unduh di sini: [Avogadro Download.](https://sourceforge.net/projects/avogadro/files/latest/download).

- Multiwfn: Multiwfn adalah alat analisis fungsi gelombang yang sangat berguna dalam kimia komputasi. Pastikan untuk menginstalnya sebelum mulai. Unduh di sini: [Multiwfn Download ]([https://sobereva.com/multiwfn/](https://drive.google.com/drive/folders/1N65AWDBRvWIPrqWztZKNEbuylwdOAg9M)).

- Orca: Orca sudah terinstal di HPC, jadi tidak perlu mengunduhnya sendiri. Pastikan akses ke HPC siap untuk digunakan.

- Gabedit: Gabedit adalah antarmuka grafis penting untuk alat kimia komputasi. Pastikan sudah terinstal. Unduh di sini: [Gabedit Download](https://gabedit.sourceforge.net/).

# Tutorial Perhitungan Muatan ESP Menggunakan Orca dan Multiwfn

## Optimasi Geometri Molekul Menggunakan Program Orca

1. Buat folder dengan nama `opt` atau apapun yang akan dijadikan tempat untuk melakukan optimasi geometri
   ```bash
   mkdir opt
   ```

   ```bash
   cd opt
   ```

2. Buat molekul menggunakan program `Avogadro`, kemudian simpan dengan nama `geom.xyz`  di dalam folder `opt`. 

3. Buat input file dengan nama `orca.in` yang berisikan:
   ```bash
   !B3LYP def2-svp opt
   %pal 
    nprocs 1 
   end
   %geom
     maxiter 9999
   end
        
   *xyzfile 0 1 geom.xyz
   ```

   Dalam hal ini, `nprocs` adalah jumlah processor yang digunakan. Dalam kasus ini akan digunakan 8 CPU core. `maxiter` dalam blok `%geom` adalah jumlah iterasi optimasi geometri. Angka `0` dan `1` pada baris `*xyzfile 0 1` masing-masing bermakna muatan dan multiplisitas spin. 

4. Buat *submission script* dengan nama `run.sh` yang berisikan:
   ```bash
   #!/bin/bash
   #SBATCH --nodes=1
       #SBATCH --ntasks=1
       #SBATCH --cpus-per-task=12
       #SBATCH --time=168:0:0
       export LD_LIBRARY_PATH=/opt/openmpi411/lib:$LD_LIBRARY_PATH
       export PATH=/opt/openmpi411/bin:$PATH
       export OMP_NUM_THREADS=1
       cd $PWD
       orca orca.in > orca.out --oversubscribe
   ```

   Dalam *script* di atas, pastikan lokasi `LD_LIBRARY_PATH` dan `PATH` untuk `openmpi` sudah benar. Cara memastikannya adalah dengan melakukan `ls` terhadap lokasi yang memungkinkan terdapat library `openmpi`, misalnya, `ls /opt/openmpi`. Jika perhitungan dilakukan di komputer sendiri yang tidak memiliki sistem antrian `slurm`, tidak perlu membuat file ini dan lewatkan saja langkah (4) dan (5).

5. Submit perhitungan dengan mengetik:
   ```bash
   sbatch run.sh
   ```

6. Lakukan perhitungan optimasi geometri dengan mengetik:
   ```bash
   orca orca.in
   ```

7. Setelah perhitungan selesai (dipastikan dengan mengetik `squeue`), cek bahwa terdapat file bernama `orca.gbw`. Bergantung nama input file yang digunakan, asalkan file tersebut memiliki ekstensi `.gbw` maka artinya perhitungan optimasi geometri berhasil dilakukan. 



## Konversi File Fungsi Gelombang Orca ke Dalam Format MOLDEN

1. Lakukan konversi file `orca.gbw` ke dalam format molden dengan mengetik:
   ```bash
   orca_2mkl orca -molden
   ```

   Dalam beberapa kasus, bisa saja perintah `orca_2mkl` tidak terdeteksi. Apabila hal ini terjadi, pastikan dahulu di mana lokasi folder utama `orca` berada. Binary `orca_2mkl` berada di lokasi yang sama dengan binary program `orca`. 

2. Jika perintah di atas berhasil, maka akan terdapat keterangan seperti berikut:
   ```bash
   Reading the input file orca.gbw               ... done
   Writing the output file orca.molden.input     ... done
   ```

3. File `orca.molden.input` berisikan informasi koefisien orbital, energi orbital, dan beragam informasi lainnya terkait dengan fungsi gelombang (basis set) yang digunakan dalam perhitungan DFT. File ini siap digunakan untuk berbagai perhitungan menggunakan program `multiwfn`. 



## Perhitungan Muatan ESP Menggunakan Program Multiwfn

1. Buat *script* dengan nama `calc_esp.sh` (nama bisa disesuaikan bergantung keinginan).
   ``` bash
   vi calc_esp.sh
   ```

   ```bash
   multiwfn $1 << EOF
   7 
   12 
   1
    
   "y"
   EOF
   ```

   Pastikan bahwa perintah `multiwfn` dapat diakses di terminal. Jika tidak, ubah `multiwfn` dalam script tersebut menjadi path di mana program`multiwfn` sudah terinstall. Simpan file tersebut dengan mengetik `:wq` seperti biasanya. 

2. Ubah perizinan file menjadi *executable* dengan mengetik:
   ```bash
   chmod +x calc_esp.sh
   ```

3. Jalankan perhitungan dengan mengetik:
   ```bash
   ./calc_esp.sh orca.molden.input
   ```

   Sesuaikan nama `orca.molden.input` dengan nama file `molden` yang diperoleh sebelumnya. File berisikan muatan ESP atom-atom dalam molekul yang digunakan akan bernama `orca.molden.chg`. Sebagai contoh, untuk molekul air sebagai berikut:
   ```bash
   O     0.946312   -0.021880    0.067820  -0.7070948892
   H     1.911733   -0.055997    0.109553   0.3535569670
   H     0.675325   -0.609363    0.786407   0.3535379223
   ```

   Kolom terakhir berisi informasi muatan ESP yang diinginkan.  

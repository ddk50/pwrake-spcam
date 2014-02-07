# Pwrake workflow demo : Subaru SprimeCam

* [Subaru Data Reduction website](http://subarutelescope.org/Observing/DataReduction/index.html)

## Reference
* Yagi et al., 2002, AJ, 123, 66
* Ouchi et al., 2004, ApJ, 611, 660

## Usage
### Compile SExtractor

    $ tar xzf sextractor-2.8.6.tar.gz
    $ cd sextractor-2.8.6
    $ ./configure
    $ make
    $ cd ..

### Compile SDFRED

    $ cd sdfred20130924_mf2
    $ ./configure
    $ make
    $ cd ..

### Download SprimeCam sample image data
(not included in this package)

    $ wget http://www.naoj.org/Observing/DataReduction/mtk/subaru_red/SPCAM/data/spcam_training_data.tar.gz
      (file size : ~800MB)
    $ tar xvzf spcam_training_data.tar.gz
    $ mv spcam_training_data/* .

### Execute Workflow

Non parallel execution

    $ rake
     or
    $ pwrake

Parallel execution

* Setup SSH for your cluster
* Store this directory to NFS or Gfarm
* Make nodefile

        $ cat nodefile
        host1.ne.jp 2
        host2.ne.jp 2
        ...

* Run Pwrake:

        $ pwrake --hostfile=nodefile

Target file: all.fits

# Setting up python enviroments with Conda
Andy@SantaCruzIntegration.com

Goal create reproducable python environments with the miniumal set of dependecies

## ref:
- [conda cheat sheet](conda-cheatsheet.pdf)
- [conda users guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

- conda basics
  - [getting started](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html)
  - [installing and using packages](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/packages.html#)
  - [conda enviroments](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html)
  

## Cheat Sheet
- environments
  * [create an env](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments)
  * list envs
    ```
    $ conda info --envs
    ```
  * use activate and deactivate to switch envs
    ```
    $  conda create --name BME-230B-hw2 python=3.7.3
    ...
    (base) $ conda info --envs
    # conda environments:
    #
    base                  *  /Users/andrewdavidson/anaconda3
    BME-230B-hw2             /Users/andrewdavidson/anaconda3/envs/BME-230B-hw2
    (base) $ conda activate BME-230B-hw2
    (BME-230B-hw2) $
    ```
  * managing packages
    + use conda install
    + use conda list
    

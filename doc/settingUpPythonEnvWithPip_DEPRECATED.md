# Setting up Python Development Environments
Andy@SantaCruzIntegration.com

<span style="color:red">This file is Deprecated</span>
It describes how to use pip to set up virtual python environments. <span style="color:red">You are better off using conda.</span> pip only knows about python. Conda is a more modern packaging system. It can handle packages with dependencies that are not written in pythong. For example libraries writting in C. see [settingUpPythonEnvWithConda.md](./settingUpPythonEnvWithConda.md)

<span style="color:red">Becareful there are several verson of python installed on our mac book pro</span>

```
$ which pip
/Library/Frameworks/Python.framework/Versions/3.4/bin/pip
$ which pip3
/Library/Frameworks/Python.framework/Versions/3.6/bin/pip3
$ which python
/usr/bin/python
$ python --version
Python 2.7.10
$ which python3
/Library/Frameworks/Python.framework/Versions/3.6/bin/python3
$ python3 --version
Python 3.6.1
$ 
```

## 1. pip cheet sheet

list: shows all all the packages currently installed. The output in the example bellow has been trunckated
```
$ pip3 list 
alabaster (0.7.10)
appnope (0.1.0)
Babel (2.5.1)
```

install or upgrade if needed
```
$ pip3 install packageName --upgrade
```

## 2. Virtual Environments

ref: [Python Tutorial: virtualenv and why you should use virtual environments](https://www.youtube.com/watch?v=N5vscPTWKOk)

Challenge: 
* How can we work on different python projects each with different dependencies?
* How can I quickly set up a new development env on a remote machine wiht the correct dependencies for my project?

Installing python global site packages is a bad idea.

Python virtualenv is similar to maven or gradle builds for Java

### 2.1. Install virtualenv
```
$ pip3 install virtualenv --upgrade
```

### 2.2 Create a directory to store virtual envs
<span style="color:red">Do not store our actual project python files in this directory.</span>
```
$ mkdir ~/workSpace/pythonEnv
```

### 2.3 Create a project specfic virtual env
use virtualenv "my project name". It will create a directory
```
$cd  ~/workSpace/pythonEnv
$ virtualenv project1_env
Using base prefix '/Library/Frameworks/Python.framework/Versions/3.6'
New python executable in /Users/andrewdavidson/workSpace/pythonEnv/project1_env/bin/python3.6
Also creating executable in /Users/andrewdavidson/workSpace/pythonEnv/project1_env/bin/python
Installing setuptools, pip, wheel...done.
$ 
```

### 2.4 How to start and stop  a virtual env
source path to project/bin/activate. Notice the prompt has changed

```
$ pwd
/Users/andrewdavidson/workSpace/pythonEnv
$ ls
project1_env/
$ source project1_env/bin/activate
(project1_env) $
```

use deactivate to stop a virtual environment. Notice the prompt and default python have been reverted to their original values.
```
(project1_env) $ which python
/Users/andrewdavidson/workSpace/pythonEnv/project1_env/bin/python
(project1_env) $ deactivate
$ which python
/usr/bin/python
$ 
```

### 2.5 how to test what environment you are in?
shell prompt should be the name of the env. You can use which to check which python and pip you are using
```
(project1_env) $ which python
/Users/andrewdavidson/workSpace/pythonEnv/project1_env/bin/python
(project1_env) $ 
```

### 2.6 How to install packages into your new env
```
project1_env) $ which pip
/Users/andrewdavidson/workSpace/pythonEnv/project1_env/bin/pip
(project1_env) $ pip --version
pip 9.0.1 from /Users/andrewdavidson/workSpace/pythonEnv/project1_env/lib/python3.6/site-packages (python 3.6)
(project1_env) $ pip install numpy
```

### 2.7 How to export all the packages and their version from our local virtual env
use pip freeze. The '--local' prevents global site packages from being included. Note we previous used pip install to add numpy and pytz. 
```
(project1_env) $ pip list
numpy (1.14.1)
pip (9.0.1)
pytz (2018.3)
setuptools (38.5.2)
wheel (0.30.0)

(project1_env) $ pip freeze --local | tee requirements.txt
numpy==1.14.1
pytz==2018.3
(project1_env) $ 
```

### 2.8 How to delete virtual environment
Its a good idea to delete environments you are not using to save disk space. its probably a good idea to use step 2.7 to create a requirements.txt file so you can restore it if needed. To delete and enviroment 
a. make sure you deactivate your virtual environment. see 2.4
b. use 'rm -rf' path to directory you created in 2.3

### 2.9 How to create a new enviroment specifying the python version and a requirements.txt file

A. create the new env specifiying the version of python you want
```
$ pwd
/Users/andrewdavidson/workSpace/pythonEnv
$ virtualenv -p /usr/bin/python2.7 py27_env
Running virtualenv with interpreter /usr/bin/python2.7
New python executable in /Users/andrewdavidson/workSpace/pythonEnv/py26_env/bin/python
Installing setuptools, pip, wheel...done.
$ 
```

B. Start the new virtual env.
```
$ source py27_env/bin/activate
(py26_env) $ python --version
Python 2.7.10
(py27_env) $ 
```

C. install the required packages
```
(py27_env) $ cat requirements.txt 
numpy==1.14.1
pytz==2018.3

(py27_env) $ pip install -r requirements.txt 
```

### 2.9 using Virtual environments in eclipse
ref: [see section "Using PyDev with virtualenv"](https://www.caktusgroup.com/blog/2011/08/31/getting-started-using-python-eclipse/)

## 3.0 Creating proejects in eclipse
ref: [https://www.rose-hulman.edu/class/csse/resources/Eclipse/eclipse-python-configuration.htm](https://www.rose-hulman.edu/class/csse/resources/Eclipse/eclipse-python-configuration.htm)


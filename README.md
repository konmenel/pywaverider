# Introduction
My code of a waverider design tool and results of the paper (still in progress).

The tool computes the optimal waverider shape for a given mission. This is done using the *Nealder-Mead Simplex* optimization algorithm in order to find the geometry with the best *L/D* ratio for the specific mission. The geometry of the waverider is computed using a method called *osculating cones*, which an extension of the *cone-derived waveriders*.

When I started to write the code, I had very little experience with python so some parts of the code are not considered pythonic. Also, there is very little optimization to the code.

# To do
- [ ] Code refactoring
- [ ] Code optimazations
- [ ] Find a reliable way to calculate the initial inputs
- [ ] Find more reliable optimization algorithms to escape local minima
- [ ] Code it in **C++**

# Dependencies
- numpy
- scipy
- fluid
- matplotlib (only for 2D plot of waverider)
- plotly (only for 3D plot of waverider)
- tringle (only for contour plots and 3D plot of waverider)
- xlsxwriter (only to output the points in .xlsx file)

# Setting up the virtual environment
The process described below is for a windows machine but it is similar for a linux machine.

You can create a virual environment using either conda or venv. For venv open Command Prompt or Powershell in the directory and type: 
```
python -m venv ./venv
```

To install the dependencies first you need to activate the environment:
```
venv\Scripts\Activate
```
After venv is activated type:
```
python -m pip install -r dependencies.txt
```
Alternatively, you can install the dependancies listed above one by one using:
```
python -m pip install <nameofdependency>
```
# Running the code
To run the code on Windows you can simple run the `missionRun.bat` file. Alternatively, you can opening a terminal on the directory and typing the following command:
```
python Code/main.py
```
A window will be displayed for you to choose the mission inputs

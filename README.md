# Introduction
My code of a waverider design tool and results of the paper (still in progress).

The tool computes the optimal waverider shape for a given mission. This is done using the Nealder-Mead Simplex optimization algorithm in order to find the geometry of the mission inputs.

When I started to write the code, I had very little experience with python so some parts of the code are not considered pythonic. Also, there is very little optimization to the code.

# Dependencies
* numpy
* scipy
* fluid
* matplotlib (only for 2D plot of waverider)
* plotly (only for 3D plot of waverider)
* tringle (only for 3D plot of waverider)
* xlsxwriter (only to output the points in .xlsx file)

# Setting up the virtual environment
## Windows
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

# Running the code
To run the code you can simple double-click the missionRun.bat file (if on Windows) or by opening a terminal on the directory and typing the following command:
```
python Code/main.py
```
A window will be displayed for you to choose the mission inputs

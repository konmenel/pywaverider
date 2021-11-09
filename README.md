# Introduction
My code of the waverider design tool and results for the paper.

When I started to write the code, I had very little experience with python so some parts of the code are not considered pythonic. Also, there is very little optimization to the code.

# Dependencies
numpy
scipy
fluid
matplotlib (only for 2D plot of waverider)
plotly (only for 3D plot of waverider)
tringle (only for 3D plot of waverider)
xlsxwriter (for output the points in .xlsx file)

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

# PyWaverider
***PyWaverider*** is a waverider design tool based on the osculating cones method using Bezier curves as the input. This is the tool used in the paper ***TBD***. 

The tool computes the optimal waverider shape for a given mission. This is done using the *Nealder-Mead Simplex* optimization algorithm in order to find the geometry with the best *L/D* ratio for the specific mission. The geometry of the waverider is computed using a method called *osculating cones*, which an extension of the *cone-derived waveriders*.

When I started to write the code, I had very little experience with python so some parts of the code are not *"pythonic"*. Also, it is very probable that there are many code optimizations to be done.

# To do
- [ ] Make it a command-line program 
- [ ] Fix all docstrings
- [ ] Code refactoring
- [ ] Code optimizations
- [ ] Find more reliable optimization algorithms to escape local minima
- [ ] Add user defined inputs (Generalize the tool)

# Dependencies
- python 3.9
- numpy
- scipy
- fluid
- matplotlib (only for 2D plot of waverider)
- plotly (only for 3D plot of waverider)
- triangle (only for contour plots and 3D plot of waverider)
- xlsxwriter (only to output the points in .xlsx file)

# Setting up the virtual environment
The process described below is for a windows machine but it is similar for a linux machine.

## Using Conda
Open Terminal, Command Prompt or Powershell in the directory and type: 
```
conda env create -f environment.yml
```
Activate the environment by typing:
```
conda activate waverider
```

## Using venv
Open Command Prompt or Powershell in the directory and type: 
```
python -m venv ./venv
```
To install the dependencies first you need to activate the environment:
```
venv\Scripts\Activate
```
After venv is activated type:
```
python -m pip install -r requirements.txt
```
Alternatively, you can install the dependancies listed above one by one using:
```
python -m pip install <nameofdependency>
```
# Running the code
To run the code just open the terminal in the root directory of the code and typing the following command:
```
python src/main.py
```
A window will be displayed for you to choose the mission inputs

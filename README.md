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
```console
$ conda env create -f environment.yml
```
Activate the environment by typing:
```console
$ conda activate waverider
```

## Using venv
Open Command Prompt or Powershell in the directory and type: 
```console
$ python -m venv ./venv
```
To install the dependencies first you need to activate the environment:
```console
$ venv\Scripts\Activate
```
After venv is activated type:
```console
$ python -m pip install -r requirements.txt
```
Alternatively, you can install the dependancies listed above one by one using:
```console
$ python -m pip install <nameofdependency>
```
# Running the code
To run the code just open the terminal in the root directory of the code and typing the following command:
```console
python src/main.py
```
A window will be displayed for you to choose the mission inputs

# Displaying results
After the optimization is finished you can see the results using:
```console
$ python src/results_run.py
```
This will lunch a window that will let you choose the output file you want to see. 

You will be asked to if you want to create a 3D plot as well along with the 2D plots (Base, Top, and Local Radii 2D plots). Additioally, you will be asked if you want to create the plot using localhost or save it as an html file. The html files are saved in the directory `plots` (will be created if doesn't exists) in the format 
`YYYY.MM.DD-hh.mm-<VISCOUS_TYPE>.html`.

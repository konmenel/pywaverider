"""Module containing the necessary functions to perform a Nelder-Mead optimazation
for a waverider object.
"""
from typing import List
import time
import sys
import numpy as np
import matplotlib.pyplot as plt
import waverider as wrdr
import config as cfg


def optimization_function(wr_inputs) -> float:
    """The function that will be used for the Nelder-Mead optimization.
    Return the L/D of the waverider if the contraints are met, otherwise
    returns 1000.
    """
    b, s, l, per_l, *le = wr_inputs
    yle, zle = le[:3], le[3:]

    wr_in = wrdr.WRInputs(b, s, l, per_l, yle, zle)

    if wr_in.check():
        return 1000.

    wr = wrdr.Waverider(wr_in)

    if wr.check():
        return 1000.

    return -wr.L_D


def NelderMead(guess: wrdr.Waverider, optVar: str, minMax=-1, parms='ANMS', maxIt=10000, Tol=1e-8) -> dict:
    """Function that performs the optimization using the Nelder-Mead Simplex Method
        Arguments:
            guess: The initial guess needed for the optimization to start (instance of waverider)
            optVar: The name of the attribute that is the objective function (string)
            minMax: Integer that defines if the problem is maximazation (-1) or minimazation (+1) (defualt=-1)
            parms: String that defines the method used to calculate the Nelder-Mead parameters a,b,c,d. 
                   Either 'ANMS' (Adaptive Nelder-Mead Simplex) or 'SNMS' (Standard Nelder-Mead Simplex) (defualt='ANMS')
            maxIt: Maximum number of iterations (defualt=10000)
            Tol: Tolerance, value of standard deviation that terminates the optimization (defualt=1e-8)

        Outputs:
            results: A dictionary with the following fields:
                'best': a list with the best L/D ratio of each iteration.
                'avg':  a list with the average L/D ratio of each iteration.
                'totalTime': the total time needed for convergence.
                'st_dev': a list with the standard deviation L/D ratio of each iteration.
                'nfeval': the number of function evaluations.
                'nRef': the number of reflection steps.
                'nExp': the number of reflection steps.
                'nInCon': the number of inside contraction steps.
                'nOutCon': the number of outside contraction steps.
                'nShr': the number of shrinkage steps.
                'nConv': the number of iterations needed for convergence.
                'WR': a dictionary with the inputs of the best Waverider.
    """
    # Nested Function for Nelder-Mead steps
    def centroid(simplex: List[wrdr.WRInputs]) -> wrdr.WRInputs:
        """Computes the centroid of the simplex"""
        cntr = 0
        for v in simplex[:-1]:
            cntr += v
        cntr /= nVar
        return cntr


    def shrinkage(simplex: List[wrdr.WRInputs]) -> List[wrdr.WRInputs]:
        """Performs shrinkage on the simplex"""
        v1 = simplex[0]
        shr = [simplex[0]]
        for vrtx in simplex[1:]:
            shr.append(v1 + d*(vrtx - v1))
        return shr


    if minMax not in (-1, 1):
        raise Exception('minMax should be either -1 or 1')

    nVar = len(guess.inputs())
    a = 1.
    if parms == 'ANMS':
        b = 1. + 2/nVar
        c = 0.75 - 1/(2*nVar)
        d = 1. - 1/nVar
    elif parms == 'SNMS':
        b = 2
        c = 0.5
        d = 0.5
    else:
        raise Exception('parms should be either "ANMS" (adaptive) or "SNMS" (standard)')

    avg = []    # List of average value for each step
    best = []   # List of best value for each step
    st_dev = [] # List of standard diviation for each step

    nfeval = 1  # Number of function evaluation
    nRef = 0    # Number of reflection steps
    nExp = 0    # Number of expantion steps
    nInCon = 0  # Number of inside contraction steps
    nOutCon = 0 # Number of outsice contraction steps
    nShr = 0    # Number of shrickage steps

    # Initialization (First Simplex)
    print('Initialization:')
    tic = time.time()
    WR = [guess]
    for i in range(nVar):
        in_ = guess.inputs()
        if in_[i] != 0:
            in_[i] *= 1.05
        else:
            in_[i] += 0.1
        in_ = in_.reInit()

        if in_.check():
            print(f'Waverider {i+1} First Check Error')
            sys.exit()
        else:
            WR.append(wrdr.Waverider(in_))

        if WR[i+1].check():
            print(f'Waverider {i+1} Second Check Error')
            sys.exit()
    WR = np.array(WR)
    print('Completed!\n')

    # Optimization
    print('Optimization:')
    for i in range(maxIt):
        # Sorting
        L_D = np.array([getattr(i, optVar) for i in WR])
        idx = np.argsort(minMax*L_D)
        L_D = L_D[idx]
        WR = WR[idx]
        in_ = [i.inputs() for i in WR]
        
        # Termination
        best.append(L_D[0])
        avg.append(np.mean(L_D))
        st_dev.append(np.sqrt(np.sum((L_D - avg[i])**2) / nVar))

        print(f'Iteration: {i}')
        print(f'Best {optVar} = {best[i]:0.6}, Average {optVar} = {avg[i]:0.6}')
        print(f'Standard Deviation = {st_dev[i]:0.6}')

        if st_dev[i] <= Tol:
            break
        
        # Centroid
        m = centroid(in_)

        # Reflection
        r = m + a*(m - in_[-1])

        if r.check():
            fr = 1e10
        else:
            wr_r = wrdr.Waverider(r)
            nfeval += 1
            if wr_r.check():
                fr = 1e10
            else:
                fr = minMax * getattr(wr_r, optVar)

        if minMax*L_D[0] <= fr < minMax*L_D[-2]:
            WR[-1] = wr_r
            nRef += 1
            print('Reflection\n')
            continue

        # Expansion
        elif minMax*L_D[0] > fr:
            e = m + b*(r - m)

            if e.check():
                fe = 1e10
            else:
                wr_e = wrdr.Waverider(e)
                nfeval += 1
                if wr_e.check():
                    fe = 1e10
                else:
                    fe = minMax*getattr(wr_e, optVar)
            
            if fe < fr:
                WR[-1] = wr_e
                nExp += 1
                print('Expansion\n')
                continue
            else:
                WR[-1] = wr_r
                nRef += 1
                print('Reflection\n')
                continue
        
        # Outside Contraction
        elif minMax*L_D[-2] <= fr < minMax*L_D[-1]:
            oc = m + c*(r - m)

            if oc.check():
                foc = 1e10
            else:
                wr_oc = wrdr.Waverider(oc)
                nfeval += 1
                if wr_oc.check():
                    foc = 1e10
                else:
                    foc = minMax*getattr(wr_oc, optVar)

            if foc <= fr:
                WR[-1] = wr_oc
                nOutCon += 1
                print('Outside Contraction\n')
                continue
            else:
    
        # Shrinkage
                in_ = shrinkage(in_)
                for j, v in enumerate(in_[1:], 1):
                    WR[j] = wrdr.Waverider(v)
                    nfeval += 1
                nShr += 1
                print('Shrinkage\n')
                continue
    
        # Inside Contraction
        elif fr >= minMax*L_D[-1]:
            ic = m - c*(r - m)

            if ic.check():
                foc = 1e10
            else:
                wr_ic = wrdr.Waverider(ic)
                nfeval += 1
                if wr_ic.check():
                    fic = 1e10
                else:
                    fic = minMax*getattr(wr_ic, optVar)

            if fic < minMax*L_D[-1]:
                WR[-1] = wr_ic
                nInCon += 1
                print('Inside Contraction\n')
                continue
            else:
    
        # Shrinkage
                in_ = shrinkage(in_)
                for j, v in enumerate(in_[1:], 1):
                    WR[j] = wrdr.Waverider(v)
                    nfeval += 1
                nShr += 1
                print('Shrinkage\n')    

    # Post optimization
    nConv = len(best)
    toc = time.time()
    totalTime = (toc - tic)/60/60 # time in hours

    results = {
        'best': best, 'avg': avg, 'totalTime': totalTime, 'st_dev': st_dev,
        'nfeval': nfeval, 'nRef': nRef, 'nExp': nExp, 'nInCon': nInCon,
        'nOutCon': nOutCon, 'nShr': nShr, 'nConv': nConv, 'WR': WR[0].todict()
        }
    return results


def initalGuesses():
    """
    A function to find Initial guesses within the constrains randomly.
    The guesses' inputs are stored from best to worst in a json file.
    """
    import random as rng
    import json
    from config_modify import config_create_window

    config_create_window()
    check = True
    count = 0
    WR = []
    while len(WR) < 100:
        try:
            count += 1
            b = rng.uniform(10, 30)
            l = rng.uniform(50, 80)
            s = rng.uniform(0.2, 0.5)*l
            per_l = rng.uniform(0.2, 1)
            yle = []
            zle = []
            for i in range(3):
                if i == 0:
                    yle.append(rng.uniform(0.1, 1))
                else:
                    yle.append(rng.uniform(0.1, 1))
                    zle.append(rng.uniform(0, 1.1))

            yle.sort()
            yle = np.array(yle)*s
            zle = np.array(zle)*l*np.tan(b*np.pi/180)
            wr_in = wrdr.WRInputs(b, s, l, per_l, yle, zle)
            check = wr_in.check()
            if not(check):
                wr = wrdr.Waverider(wr_in)
                check = wr.check()
                if not(check):
                    WR.append(wr)
        except KeyboardInterrupt:
            print(count)
            print(len(WR))
            sys.exit()
    
    WR = np.array(WR)
    L_D = np.array([i.L_D for i in WR])
    idx = np.argsort(-L_D)
    L_D = L_D[idx]
    WR = WR[idx]

    json_dict = {'L_D': L_D.tolist(), 'guesses': [i.todict() for i in WR]}

    filename = "WaveriderGuesses.json"
    with open(filename, 'w') as output:  # Overwrites any existing file.
        json.dump(json_dict, output, indent=4)


def results_to_screen(results: dict, plot_3D: bool=True, htmlExport: bool=False):
    """A function that plot and prints the results to the screen.
    Arguments:
        results: The dictionary of the results
        plot_3D: A boolean to create a 3D plot if True. By default True.
        htmlExport: A boolean to export an offline html file of the 3D plot instead of
                    opening it in the browser. By default False.
    """
    WR = wrdr.Waverider(wrdr.WRInputs(**results['WR']))
    iterations = [i for i in range(results['nConv'])]

    # Print Results
    convTime = f'{int(np.floor(results["totalTime"]))} hours and {int((results["totalTime"] - np.floor(results["totalTime"]))*60)} minutes' 
    print(f'Best L_D: {WR.L_D:0.5}')
    print(f'Method for Viscous Effects: {cfg.CONFIG["Viscous Method"]}')
    print(f'Time needed for convergence: {convTime}')
    print(f'Number of Function Evaluations: {results["nfeval"]}')
    print('Number of Iterations:')
    print(f'\tReflection: {results["nRef"]}')
    print(f'\tExpansion: {results["nExp"]} ')
    print(f'\tOutside Contraction: {results["nOutCon"]}')
    print(f'\tInstide Contraction: {results["nInCon"]}')
    print(f'\tShrinkage: {results["nShr"]}')
    print(f'\tTotal: {results["nConv"]}')
    print('----------------------------------------')
    print('Best Waverider:')
    print(WR)

    # Poria Sigklisi
    _, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8.5, 6))
    mananger = plt.get_current_fig_manager()
    mananger.set_window_title(cfg.CONFIG['Viscous Method'])
    plt.suptitle('Convergence')

    ax1.plot(iterations, results['best'], label='Best ')
    ax1.plot(iterations, results['avg'], label='Average')
    ax1.set_ylabel(cfg.optVar)
    ax1.set_xlim(left=0)
    ax1.legend(loc='lower right')
    ax1.grid()

    ax2.semilogy(iterations, results['st_dev'])
    ax2.set_xlabel('Interations')
    ax2.set_ylabel('Standard Deviation')
    ax2.grid()

    # Waverider Plots
    WR.plotAll(plot_3D, htmlExport)

    plt.show()



if __name__ == "__main__":
    def testing() -> None:
        pass
    
    testing()

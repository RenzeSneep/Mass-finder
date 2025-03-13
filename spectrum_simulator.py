import itertools
import math
import operator as op
from functools import reduce
import all_elements
import numpy as np
import matplotlib.pyplot as plt
import mass_calculator

def spectrum_simulator(molecular_formula, charge):
    '''
    Simulates a mass spectrum based on a molecular formula and charge.

    Args:
        molecular_formula (str): The molecular formula to get the ion mass from.
        charge (str): Either + or -.

    Returns:
        peaks (list[float]): The mass peaks in the spectrum.
        heights (list[float]): The heigths of corresponding peaks.
    '''
    significant_atoms, other_atoms, rest_masses, rest_abundances, elements=get_significant_atoms(molecular_formula)
    possible_isotope_combinations=[]
    abundances_possible_combinations=[]

    for i in range(len(elements)):
        masses=itertools.combinations_with_replacement(elements[i].natural_isotopes, significant_atoms[i])
        masses=[list(j) for j in masses]
        possible_isotope_combinations.append(masses)
        abun=itertools.combinations_with_replacement(elements[i].natural_abundance, significant_atoms[i])
        abun=[list(k) for k in abun]
        abundances_possible_combinations.append(abun)
    for i in range(len(possible_isotope_combinations)):
        for j in range(len(possible_isotope_combinations[i])):
            possible_isotope_combinations[i][j]=possible_isotope_combinations[i][j]+rest_masses[i]
    for i in range(len(abundances_possible_combinations)):
        for j in range(len(abundances_possible_combinations[i])):
            abundances_possible_combinations[i][j]=abundances_possible_combinations[i][j]+rest_abundances[i]
    possible_combinations_per_set=[]
    for i in range(len(possible_isotope_combinations)):
        y=[]
        for j in possible_isotope_combinations[i]:
            isotopes = []
            for k in j:
                if k not in isotopes:
                    isotopes.append(k)
            x=nPr(len(j), len(j))
            for i in isotopes:
                x = int(x/math.factorial(j.count(i)))
            y.append(x)
        possible_combinations_per_set.append(y)
            
    possible_combinations_entire_molecule=combinations_overall(possible_isotope_combinations)
    number_of_possible_combinations_entire_molecule=amounts_overall(possible_combinations_per_set)
    abundance=combinations_overall(abundances_possible_combinations)
    peaks=sum_combinations(possible_combinations_entire_molecule, charge)
    heights=prod_combinations(abundance, number_of_possible_combinations_entire_molecule)
    return peaks, heights

def get_significant_atoms(molecular_formula):
    elements_as_string,amount_of_element,molecular_mass=mass_calculator.mass_and_element_calculator(molecular_formula, "+")
    elements=[]
    for i in elements_as_string:
        for j in all_elements.elements_list:
            if i == j.symbol:
                elements.append(j)
    
    significant_atoms=[]
    other_atoms=[]
    for i in range(len(amount_of_element)):
        x=amount_of_element[i]-int(round(((1-max(elements[i].natural_abundance))*500)))
        if x < 0:
            x=0
        y=amount_of_element[i]-x
        significant_atoms.append(y)
        other_atoms.append(x)
    rest_masses=[]
    for i in range(len(other_atoms)):
        rest_mass=[]
        for j in range(other_atoms[i]):
            rest_mass.append(elements[i].natural_isotopes[0])
        rest_masses.append(rest_mass)
    rest_abundances=[]
    for i in range(len(other_atoms)):
        rest_abundance=[]
        for j in range(other_atoms[i]):
            rest_abundance.append(elements[i].natural_abundance[0])
        rest_abundances.append(rest_abundance)
    return significant_atoms, other_atoms, rest_masses, rest_abundances, elements

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom

def nPr(n, r):
    return reduce(op.mul, range(n, n-r, -1), 1)

def combinations_overall(possible_isotope_combinations):
    combs=itertools.product(*possible_isotope_combinations)
    combs=[list(i) for i in combs]
    possible_combinations_entire_molecule=[]
    for i in combs:
        possible_combinations_entire_molecule.append(list(itertools.chain.from_iterable(i)))
    return possible_combinations_entire_molecule

def amounts_overall(amount):  
    ams=itertools.product(*amount)
    ams=[list(i) for i in ams]
    number_of_possible_combinations_entire_molecule=[]
    for i in ams:
        x=np.prod(i)
        number_of_possible_combinations_entire_molecule.append(x)
    return number_of_possible_combinations_entire_molecule

def sum_combinations(values, charge):
    y=[]
    for i in values:
        if charge == "+":
            x=sum(i)- mass_calculator.electron_mass
        elif charge == "-":
            x=sum(i) + mass_calculator.electron_mass
        else:
            print("no charge")
        y.append(x)
    return y

def prod_combinations(values, number_of_possible_combinations_entire_molecule):
    y=[]
    for i in range(len(values)):
        x=np.prod(values[i])*number_of_possible_combinations_entire_molecule[i]
        y.append(x)
    return y


def plot(peak, height, molecular_formula, charge):
    heightnorm=[]
    heightmax=max(height)
    
    for i in height:
        sig=(i/heightmax)*100
        heightnorm.append(sig)
        
    mass=mass_calculator.mass_and_element_calculator(molecular_formula, charge)[2]
    low=mass
    for i in range(len(heightnorm)):
        if heightnorm[i]>0.2:
            if peak[i] < low:
                low=peak[i]  
    low-=0.3  
    up=0
    for i in range(len(heightnorm)):
        if heightnorm[i]>0.2:
            if peak[i] > up:
                up=peak[i]
    up+=0.3
    plt.xlim(low,up)
    plt.ylim(0,110)
    markerline, stemlines, baseline = plt.stem(peak,heightnorm,markerfmt=" ", basefmt=" ", linefmt = 'k')
    plt.setp(stemlines, linewidth=0.3)
    plt.show()
    

if __name__ == "__main__":
    molecular_formula = "C71H102N2O9S2Na"
    charge = "+"
    peak, height = spectrum_simulator(molecular_formula, charge)
    plot(peak, height, molecular_formula, charge)
    significant_atoms, other_atoms, rest_masses, rest_abundances, elements = get_significant_atoms(molecular_formula)

    

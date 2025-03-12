import all_elements

electron_mass = 0.000548

def mass_and_element_calculator(molecular_formula, charge):  
    '''
    Separates a chemical formula into elements and their coefficients. Also calculates the molecular mass based on the charge given.
    
    Args:
        molecular_formula (str): A string of elements and their coefficients.
        charge (str): Either 0, + or - for neutral, positive or negative.
    
    Returns:
        elements_in_formula (list[str]): A list of the elements present in the molecule.
        amount_of_element (list[int]): The amounts of the corresponding element.
        molecular_mass (float): The mass of the molecule, given its charge.
    '''
    
    fragment=[]
    elements_in_formula=[]
    amount_of_element=[]
    
    position=0
    for i in range(len(molecular_formula)):
        if molecular_formula[i].isupper():
            a=molecular_formula[position:i]
            position=i
            fragment.append(a)
    a=molecular_formula[position:]
    fragment.append(a)
    fragment.pop(0)
    
    for i in fragment:
        position2=0
        dig=0
        for j in range(len(i)):
            if i[j].isdigit():
                b=i[position2:j]
                position2=j
                elements_in_formula.append(b)
                dig=1
                break
        c=i[position2:]
        if dig==0:
            elements_in_formula.append(i)
            c=1
        amount_of_element.append(int(c))
        
    molecular_mass=0
    for i in range(len(elements_in_formula)):
        molecular_mass+=all_elements.atoms[elements_in_formula[i]]*int(amount_of_element[i])
    if charge == 0:
        pass
    elif charge == "+":
        molecular_mass -= electron_mass
    elif charge == "-":
        molecular_mass += electron_mass
    return elements_in_formula,amount_of_element,molecular_mass

def minus_H(molecular_formula):
    '''
    Subtracts one H atom from a molecular formula.
    
    Args:
        molecular_formula (str): A molecular formula.
        
    Returns:
        molecular_formula_minus_H (str): The molecular formula minus one H.
    
    '''
    fragment=[]
    elements_in_formula=[]
    number=[]
    
    pos=0
    for i in range(len(molecular_formula)):
        if molecular_formula[i].isupper():
            a=molecular_formula[pos:i]
            pos=i
            fragment.append(a)
    a=molecular_formula[pos:]
    fragment.append(a)
    fragment.pop(0)
    for i in fragment:
        pos2=0
        dig=0
        for j in range(len(i)):
            if i[j].isdigit():
                b=i[pos2:j]
                pos2=j
                elements_in_formula.append(b)
                dig=1
                break
        c=i[pos2:]
        if dig==0:
            elements_in_formula.append(i)
            c=1
        number.append(c)
    molecular_formula_minus_H = ""
    for j in range(len(number)):
        if number[j] == 1:
            number[j] = ""
    for i in range(len(elements_in_formula)):
        if elements_in_formula[i] == "H":
            number[i] = int(number[i])-1
        molecular_formula_minus_H += elements_in_formula[i] + str(number[i])

    return molecular_formula_minus_H

if __name__ == "__main__":
    molecular_formula = "C8H10N4O2"
    charge = "+"
    elements_in_formula,amount_of_element,molecular_mass = mass_and_element_calculator(molecular_formula, charge)

    print(f"Elements:\n{elements_in_formula}\nAmounts:\n{amount_of_element}\nMass:\n{molecular_mass}")
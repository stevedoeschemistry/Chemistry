# -*- coding: utf-8 -*-

"""
Created on Thu Apr  4 11:10:00 2024

@author: SteveAponte
"""

import csv
import os
import pandas as pd

#Creates a dictionary for the periodic table from a csv file
with open('ptable.csv', 'r') as periodic_table_file:
    csv_reader = csv.DictReader(periodic_table_file)
    periodic_table = [row for row in csv_reader]
df = pd.read_csv('ptable.csv')
df['Atom'] = 'H','He', 'Li','Be','B','C','N','O','F','Ne','Na', 'Mg','Al','Si', 'P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og','Uuu'

#takes a molecule and returns a dictionary of the atoms as keys and the subscript as a value, if no subscript value is set to 1
def sep_mol(molecule):
    mol_dict = {}
    current_atom = ""
    atom_list = []
    sub_list = []
    for x in molecule:
        if x.isupper() == True:
            if current_atom != "":
                atom_list.append(current_atom)
                sub_list.append(1)
                current_atom = ""
                current_atom = current_atom + x
            else:
                current_atom = current_atom + x
        elif x.islower() == True:
            current_atom = current_atom + x
        else:
            atom_list.append(current_atom)
            sub_list.append(int(x))
            current_atom = ""
    if current_atom == "":
        pass
    else:
        atom_list.append(current_atom)
        sub_list.append(1)
    for i in range(len(atom_list)):
            mol_dict[atom_list[i]] = sub_list[i]    
    return mol_dict, atom_list
 

#returns the molecular weight of a given molecule
def molar_mass(molecule):
    x,y = sep_mol(molecule)
    atom_sum_multiplied = 0
    current_key = []
    for atom in y:
        atom_sum = df.loc[df['Atom'] == atom, 'atomic_mass'].iloc[0]
        atom_sum_multiplied += atom_sum * x[atom]

    return atom_sum_multiplied
    """
    #returns the rows where the Atoms appear
    for atom in molecule:
        #.loc locates the specific atomic_mass and iloc seperates out the row number and spaces 
        atom_mass = df.loc[df['Atom'] == atom, 'atomic_mass'].iloc[0]
        atom_sum += atom_mass
        
    return(float(atom_sum))
    """
def stoichiometry (start_chem,start_mol,end_chem,end_mol,start_mass=None,start_vol = None):
    #Take an entered value in grams or L and give you the conversion through a chemical reaction.
    molecule_start = molar_mass(start_chem)
    molecule_end = molar_mass(end_chem)
    
    if start_vol == None:
        return str((molecule_start*start_mass)*(start_mol/end_mol)/molecule_end) + "g"
    else: return str(start_vol*(start_mol/end_mol)) + "L"

#This function lets you calculate either p,n,T,or V using the ideal gas law. The function calculates for the variable not enetered, Pressure (atm) = p;  volume (L) = v; temperature (K) = temp; n = moles    
def ideal_gas_law(p = None,v = None,n = None, temp=None):
    r = 0.08206
    if p == None :
        return str((n * r * temp)/v) + "atm"
    elif n == None:
        return str((p*v)/(r*temp)) + "moles"
    elif temp == None:
        return str((p*v)/(n*r)) + "K"
    elif v == None:
        return str((n * r * temp)/p) + "L"

           
print(df.head)


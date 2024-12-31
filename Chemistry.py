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

#returns the molecular weight of a given molecule
def molar_mass(atom_list):
    atom_sum = 0
    #returns the rows where the Atoms appear
    for atom in atom_list:
        #.loc locates the specific atomic_mass and iloc seperates out the row number and spaces 
        atom_mass = df.loc[df['Atom'] == atom, 'atomic_mass'].iloc[0]
        atom_sum += atom_mass
        
    return(float(atom_sum))

def stoichiometry (start_chem,start_mol,end_chem,end_mol,start_mass=None,start_vol = None):
    #Take an entered value in grams or L and give you the conversion through a chemical reaction.
    molecule_start = molar_mass(start_chem)
    mol_start = start_mol
    molecule_end = molar_mass(end_chem)
    mol_end = end_mol
    starting_mass = start_mass
    starting_volume = start_vol
    
    if starting_volume == None:
        return str((molecule_start*starting_mass)*(mol_start/mol_end)/molecule_end) + "g"
    else: return str(start_vol*(mol_start/mol_end)) + "L"

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

print(molar_mass(["K", "O", "H"]))
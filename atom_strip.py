# -*- coding: utf-8 -*-
"""
Created on Fri Aug 04 17:16:26 2017

@author: av.sofronov
"""
from pdbremix import pdbatoms, util
from datetime import date 
import os, glob, string

def diviser(original_pdb, soup):
  util.goto_dir(str(date.today()))      
  util.goto_dir(original_pdb)
  soup.write_pdb(original_pdb)
  vals = [0, 10, 20, 30, 40]
  for v in vals:
    soup = pdbatoms.Soup(original_pdb)
    from_val = v
    to_val = from_val + 10
    for res in soup.residues():
      for atom in res.atoms():
        if atom.bfactor <= from_val or atom.bfactor > to_val:
          res.erase_atom(atom.type)
          
    soup.write_pdb('('+str(from_val)+','+str(to_val)+']'+original_pdb)  
  
  soup = pdbatoms.Soup(original_pdb)  
  for res in soup.residues():
    for atom in res.atoms():
      if atom.bfactor != 0:
        res.erase_atom(atom.type)
  soup.write_pdb('0_'+original_pdb)      
  
  
parent_dir = os.getcwd()
    
for file in glob.glob('asa_*'):
  original_pdb = string.replace(file, 'asa_', '')
  soup = pdbatoms.Soup(file)
  diviser(file, soup)
  print (original_pdb+' done')
  os.chdir(parent_dir)

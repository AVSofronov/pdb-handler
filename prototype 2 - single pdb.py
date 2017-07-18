from pdbremix import pdbatoms, util, asa, data, pdbtext
#from pdbremix import fetch
from datetime import datetime, date 
import matplotlib.pyplot as plt
import collections 

def counter(atom, borders, distribution):
  elem = atom.element
  if atom.bfactor == 0:
    distribution[elem, '=0'] += 1
    distribution['t', '=0'] += 1
  for b in borders:
    if atom.bfactor > b:
      distribution[elem, b] += 1
      distribution['t', b] += 1
                  
#asa calculation                  
def calculate_asa(soup):            
  atoms = []
  for res in soup.residues():
    if res.type not in data.solvent_res_types:
      atoms.extend(res.atoms())
  pdbatoms.add_radii(atoms)
  areas = asa.calculate_asa_optimized(atoms, 1.4)
  for a, area in zip(atoms, areas):
    a.asa = area
    a.bfactor = area
  for a in soup.atoms():
    if not hasattr(a, 'asa'):
      a.bfactor = 0.0
  soup.write_pdb('asa_'+original_pdb+'.pdb')                  

                                   
#function to make distribution of number of atoms by ASA values
# for each element and in total
def elements_asa_distr(soup, borders, mod):

  H=0
  C=0
  O=0
  N=0
  S=0
  P=0
  for a in soup.atoms():
    if a.element == 'H':
      H += 1
    elif a.element == 'C':
      C += 1
    elif a.element == 'O':
      O += 1
    elif a.element == 'N':
      N += 1
    elif a.element == 'S':
      S += 1
    elif a.element == 'P':
      P += 1
    
  util.goto_dir(protein_name+' '+method+' '+original_pdb)

  #creating of dictionary  
  considered_types = ['O', 'N', 'C', 't']
  
  distribution = collections.OrderedDict()
  reldistribution = collections.OrderedDict()
  types = list(considered_types)
  types.remove('t')
  
  for t in considered_types:
    distribution[t,'=0'] = 0 
    for b in borders:
      distribution[t,b] = 0
  
  #calculating the distribution    var H_N  
  if mod == 'H_N':
    for a in soup.atoms():
      valid_atom_types = [
#              "O",                  
              "CG2",
              "ND2",
              "NH1",
              "NH2",
              "NZ",
              "OD1",
              "OE1",
              "OE2",
              "OG",
              "OG1",
              "OH",
              "SG",
              "OXT"
              ]
      """near_valid_atom_types = [
              "CG1",
              "CB",
              "CD1",
              "CD2",
              "CE",
              "NE2",
              ]"""     
      res = a.res_type
      atom_type = a.type
      if a.element not in ['S','P']:
        if atom_type in valid_atom_types:      
          counter(a, borders, distribution)
        elif atom_type == 'CB':
          if res == 'ALA':
            counter(a, borders, distribution)
        elif atom_type == 'CD1':
          if res in ['LEU','ILE']:
            counter(a, borders, distribution)
        elif atom_type == 'CD2':
          if res == 'LEU':
           counter(a, borders, distribution)
        elif atom_type == 'CE':
          if res == 'MET':
            counter(a, borders, distribution)
        elif atom_type == 'NE2':
          if res == 'GLN':
            counter(a, borders, distribution)
        elif atom_type == 'CG1':
          if res == 'VAL':
            counter(a, borders, distribution)
            
  elif mod == 'ALL':
    for a in soup.atoms():      
      if a.element != 'S':
        counter(a, borders, distribution)
  else:
    return 'nothing'
                      
  #relative distribution
  for t in considered_types:    
    for b in borders:
      try:
        reldistribution[t,b] = round(float(distribution[t,b])/float(distribution['t',b] ),2)
        reldistribution[t,'=0'] = round(float(distribution[t,'=0'])/float(distribution['t','=0']),2)
      except ZeroDivisionError:
        print ('ZeroDivisionError with '+t+str(b))
        
  #result representation
  total = H + C + O + N + S + P
  tC = distribution['C','=0']+distribution['C',0]
  tN = distribution['N','=0']+distribution['N',0]
  tO = distribution['O','=0']+distribution['O',0]
  t = distribution['t','=0']+distribution['t',0]

  result = open('result.txt','w')  
  for a in distribution:    
    c = (a,distribution[a],reldistribution[a])
    print c
    b = str(c)+'\n'
    result.write(b)
  result.write('\natoms counted: '+str(total)+' of them with 1 H_N: '+str(t))
  result.write('\n'+'H: '+str(H)+'\n')
  result.write('C: '+str(C)+' of them with 1 H_N: '+str(tC)+'\n')
  result.write('O: '+str(O)+' of them with 1 H_N: '+str(tO)+'\n')
  result.write('N: '+str(N)+' of them with 1 H_N: '+str(tN)+'\n')
  result.write('S: '+str(S)+'\n')
  result.write('P: '+str(P)+'\n')    
  result.write('\nProtein name: '+protein_name)
  result.write('\nMethod: '+method)
  result.write('\nPDB name: '+original_pdb)
  result.close()
  

  xticks = ['=0']
  for b in borders:
    xticks.append('>'+str(b))

  #barchart for elements
  tr = collections.OrderedDict()
  for t in types:    
    tr['=0'] = reldistribution[t, '=0']     
    for b in borders:
      tr[b] = reldistribution[t,b]
    x = range(len(tr.keys()))
    y = tr.values()
    plt.bar(x,y,align = 'center')    
    plt.xticks(range(len(tr)), xticks)
    plt.suptitle(protein_name+' '+method+' '+original_pdb)
    plt.title('distribution of '+t+' atoms by ASA')
    for j,i in zip(x,y):
      plt.annotate(i, xy=(j,i), xytext=(-10,2), textcoords='offset points')
    plt.xlabel('ASA')
    plt.ylabel('%')
    plt.savefig(t+'.jpeg')
    plt.show()    
    
  #barchart total 
  total = collections.OrderedDict()
  total['=0'] = distribution['t', '=0']
  for b in borders:
    total[b] = distribution['t',b]
  #print ('Distribution of total atoms by ASA in ' +original_pdb)
  #print tr
  x = range(len(total.keys()))
  y = total.values() 
  plt.bar(x,y,align = 'center')
  plt.xticks(range(len(total)), xticks)
  plt.suptitle(protein_name+' '+method+' '+original_pdb)
  plt.title('distribution of total atoms by ASA')
  for j,i in zip(x,y):
    plt.annotate(i, xy=(j,i), xytext=(-8,1), textcoords='offset points')
  plt.xlabel('ASA')
  plt.ylabel('#atoms')
  plt.savefig('total.jpeg')
  plt.show()
  plt.close()


start = datetime.now()

util.goto_dir('results')
#fetch.get_pdbs_with_http('1btv') #To download pdb directly from web


# >>>>>>>>>>>>ENTER PARAMETERS HERE<<<<<<<<<<<<<<    
                    
mode = 'H_N'    
borders = [0, 10, 20, 30, 40, 50]   
protein_name = 'test'
method = 'test'
original_pdb = 'calcineurin nmr 2jog'

#>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<


#    DISABLE/ENABLE BELOW FUNCTIONS AS YOU NEED

#pdbtext.clean_pdb(original_pdb+'.pdb', 'clean_'+original_pdb+'.pdb') 

#our_soup = pdbatoms.Soup('clean_'+original_pdb+'.pdb')
#calculate_asa(our_soup) 

asa_soup = pdbatoms.Soup('asa_'+original_pdb+'.pdb') 
util.goto_dir(str(date.today()))
elements_asa_distr(asa_soup, borders, mode)

finish = datetime.now()
print finish - start

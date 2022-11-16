#Script 1

#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#This notebook houses the PubchemPy Script
#install PubchemPy if needed then import
#import sys
#!{sys.executable} -m pip install pubchempy


# In[19]:


from pprint import pprint
import pubchempy as pcp


# In[2]:


#Proof of concept tutorial examples
#results = pcp.get_compounds('Glucose', 'name')
#print (results)


# In[3]:


#uncomment below to see the top 5 lines of the read in file
with open('InChI_Identifiers.txt') as f:
    inchicontents=[line.strip() for line in f.readlines()]
    #print(inchicontents[:5])


# In[14]:


#uncomment below to see the top 5 lines of the read in file
with open('SMILES_Identifiers_ind.txt') as f:
    SMILEScontents=[line.strip() for line in f.readlines()]
    print(SMILEScontents[:5])


# In[5]:


#list comprehension version using tuples to link inchiIDs and CIDs
def InChIKeyConverter(inchikey_identifier_list):
    return [(inchiID, pcp.get_compounds(inchiID, 'inchikey')) for inchiID in inchikey_identifier_list]


# In[15]:


#list comprehension version using tuples to link SMILES and CIDs
def SMILESConverter(SMILES_identifier_list):
    return [(SMILESid, pcp.get_compounds(SMILESid, 'smiles')) for SMILESid in SMILES_identifier_list]


# In[7]:


#InChIKeyConverter proof of concept on the first 50 elements
pprint(InChIKeyConverter(inchicontents[:50]))


# In[20]:


#Convert the full InChI list
pprint(InChIKeyConverter(inchicontents))


# In[21]:


#SMILESConverter proof of concept on the first 5 elements
print(SMILESConverter(SMILEScontents[:5]))


# In[31]:


#Convert the full SMILES list
#bugged
pprint(SMILESConverter(SMILEScontents))


# In[ ]:





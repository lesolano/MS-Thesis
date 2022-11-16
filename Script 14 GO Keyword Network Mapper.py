#Script 14

#!/usr/bin/env python
# coding: utf-8

# In[1]:


#This notebook is intended to house the fully commented and ready to run functions

#import libraries
import networkx
import obonet
import csv
import pandas as pd
import numpy as np
import functools as ft
from pprint import pprint


# In[2]:


#Save the OBO formatted GO network from a url to a obonet graph
url = 'http://purl.obolibrary.org/obo/go.obo'
graph = obonet.read_obo(url)


# In[3]:


#KWNameQuery takes a user input keyword and searches all GO Terms for that keyword.
#Hits are stored in a list called keyword_nodes
def KWNameQuery(userinput):
    name_nodes = []

    for node in graph.nodes:
        if (userinput in graph.nodes[node]['name']):
            name_nodes.append(node)
    return(name_nodes)


# In[4]:


#KWDefQuery takes a user input keyword and searches all GO Terms definitions for that keyword.
#Hits are stored in a list called def_nodes
def KWDefQuery(userinput):
    def_nodes = []

    for node in graph.nodes:
        if (userinput in graph.nodes[node]['def']):
            def_nodes.append(node)
    return(def_nodes)


# In[5]:


#ConcatKeywordHits takes in the userinput, calls KWNameQuery and KWDefQuery seperately with the same input
#then uses set to return only unique hits
def ConcatKeywordHits(userinput):
    return(set(KWNameQuery(userinput) + KWDefQuery(userinput)))


# In[6]:


#ID to Name and Name to ID fetch either ID or Name, respectively from the other. Exact matches only
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}


# In[7]:


#Function outputs parent terms to an input term
def Parentfinder(term):
    node = term
    for child, parent, key in graph.out_edges(node, keys=True):
        return(f'• {id_to_name[child]} ⟶ {key} ⟶ {id_to_name[parent]}')


# In[8]:


#Function outputs parent terms to an input term
def Childfinder(term):
    node = term
    for parent, child, key in graph.in_edges(node, keys=True):
        return(f'• {id_to_name[child]} ⟵ {key} ⟵ {id_to_name[parent]}')


# In[9]:


#Demonstrating use of the function with lipid and saving lipid nodes for future work
lipid_nodes = ConcatKeywordHits('lipid')


# In[10]:


lipid_nodes


# In[11]:


#Identifies superterm relationships to the input GOTERM
def SupertermIdentifier(GOTERM):
    return(sorted(id_to_name[superterm] for superterm in networkx.descendants(graph, str(GOTERM))))


# In[12]:


def SubtermIdentifier(GOTERM):
    return(sorted(id_to_name[subterm] for subterm in networkx.ancestors(graph, str(GOTERM))))


# In[13]:


def AllPathsToRoot (GOTERM, ROOT_GOTERM):
    paths = networkx.all_simple_paths(
    graph,
    source=GOTERM,
    target=ROOT_GOTERM
    )
    for path in paths:
        return('•', ' ⟶ '.join(id_to_name[node] for node in path))


# In[14]:


#One function to house all the relationships per GOTERM
#print(NetworkMapper("GO:0001727"))
def NetworkMapper(GOTERM):
    relationships={
        'GOTERM':GOTERM,
        'GONAME':id_to_name[GOTERM],
        'Parents': [],
        'Children': [],
        'Superterms': [],
        'Subterms': [],
        'PathtoMF': [],
        'PathtoBP': [],
        'PathtoCC': []
    }
    relationships['Parents']= Parentfinder(GOTERM)
    relationships['Children']= Childfinder(GOTERM)
    relationships['Superterms']= SupertermIdentifier(GOTERM)
    relationships['Subterms']= SubtermIdentifier(GOTERM)
    relationships['PathtoMF']= AllPathsToRoot(GOTERM, name_to_id['molecular_function'])
    relationships['PathtoBP']= AllPathsToRoot(GOTERM, name_to_id['biological_process'])
    relationships['PathtoCC']= AllPathsToRoot(GOTERM, name_to_id['cellular_component'])
    return(relationships)


# In[15]:


#create a list to store the dictionaries. Term by term, append the result of NetworkMapper.
#This must be output to a tsv or csv table
LipidRelationships=[]
for term in lipid_nodes:
    LipidRelationships.append(NetworkMapper(term))

#print(LipidRelationships)


# In[16]:


#Import DEG outputs from DESEQ2 as dataframes
HeLa0RvsHeLaControl_DEG = pd.read_csv('VERIFYHeLa0RvsHeLaControl.csv')
HeLa8RvsHeLaControl_DEG = pd.read_csv('VERIFYHeLa8RvsHeLaControl.csv')
HeLa8RvsHeLa0R_DEG = pd.read_csv('VERIFYHeLa8RvsHeLa0Rv2.csv')


# In[17]:


#Designate columns of interest. baseMean columns are equivalent, so we only need 1.
HeLa0RvsHeLaControl_DEGstats =HeLa0RvsHeLaControl_DEG[['Gene', 'baseMean', 'log2FoldChange', 'padj']]
HeLa8RvsHeLaControl_DEGstats =HeLa8RvsHeLaControl_DEG[['Gene', 'log2FoldChange', 'padj']]
HeLa8RvsHeLa0R_DEGstats =HeLa8RvsHeLa0R_DEG[['Gene', 'log2FoldChange', 'padj']]


# In[18]:


#create a list of imported dataframes and use a lambda function to left and right join into 1 df
dfs =[HeLa0RvsHeLaControl_DEGstats,HeLa8RvsHeLaControl_DEGstats,HeLa8RvsHeLa0R_DEGstats]
df_final = ft.reduce(lambda left, right: pd.merge(left, right, on='Gene'), dfs)


# In[19]:


#Read in tsv files from QuickGO
df1= pd.read_csv('QuickGO-annotations1-600.tsv', sep='\t', 
                dtype={"GENE PRODUCT DB":'str',
                       "GENE PRODUCT ID":'str',
                       "SYMBOL":'str',
                       "QUALIFIER":'str',
                       "GO TERM":'str',
                       "GO NAME":'str',
                       "ECO ID":'str',
                       "GO EVIDENCE CODE":'str',
                       "REFERENCE":'str',
                       "WITH/FROM":'str',
                       "TAXON ID":'int',
                       "ASSIGNED BY":'str',
                       "ANNOTATION EXTENSION":'str',
                       "GENE_PRODUCT_TYPE":'str',
                       "GO ASPECT":'str'})
df2= pd.read_csv('QuickGO-annotations601-921.tsv', sep='\t', 
                dtype={"GENE PRODUCT DB":'str',
                       "GENE PRODUCT ID":'str',
                       "SYMBOL":'str',
                       "QUALIFIER":'str',
                       "GO TERM":'str',
                       "GO NAME":'str',
                       "ECO ID":'str',
                       "GO EVIDENCE CODE":'str',
                       "REFERENCE":'str',
                       "WITH/FROM":'str',
                       "TAXON ID":'int',
                       "ASSIGNED BY":'str',
                       "ANNOTATION EXTENSION":'str',
                       "GENE_PRODUCT_TYPE":'str',
                       "GO ASPECT":'str'})


# In[20]:


GOTermGeneAssociator ={
'GO:0000030':[],
 'GO:0000139':[],
 'GO:0000170':[],
 'GO:0000329':[],
 'GO:0000421':[],
 'GO:0000740':[],
 'GO:0000741':[],
 'GO:0000742':[],
 'GO:0001575':[],
 'GO:0001727':[],
 'GO:0001765':[],
 'GO:0001786':[],
 'GO:0002450':[],
 'GO:0002457':[],
 'GO:0002468':[],
 'GO:0002469':[],
 'GO:0002470':[],
 'GO:0002471':[],
 'GO:0002472':[],
 'GO:0002473':[],
 'GO:0002475':[],
 'GO:0002493':[],
 'GO:0002494':[],
 'GO:0002598':[],
 'GO:0002599':[],
 'GO:0002600':[],
 'GO:0002933':[],
 'GO:0003373':[],
 'GO:0003374':[],
 'GO:0004308':[],
 'GO:0004376':[],
 'GO:0004377':[],
 'GO:0004378':[],
 'GO:0004481':[],
 'GO:0004576':[],
 'GO:0004583':[],
 'GO:0004620':[],
 'GO:0004623':[],
 'GO:0004629':[],
 'GO:0004859':[],
 'GO:0004879':[],
 'GO:0005227':[],
 'GO:0005253':[],
 'GO:0005261':[],
 'GO:0005319':[],
 'GO:0005354':[],
 'GO:0005543':[],
 'GO:0005544':[],
 'GO:0005545':[],
 'GO:0005548':[],
 'GO:0005633':[],
 'GO:0005635':[],
 'GO:0005637':[],
 'GO:0005640':[],
 'GO:0005641':[],
 'GO:0005740':[],
 'GO:0005741':[],
 'GO:0005743':[],
 'GO:0005758':[],
 'GO:0005765':[],
 'GO:0005774':[],
 'GO:0005777':[],
 'GO:0005778':[],
 'GO:0005789':[],
 'GO:0005790':[],
 'GO:0005811':[],
 'GO:0005886':[],
 'GO:0006012':[],
 'GO:0006050':[],
 'GO:0006071':[],
 'GO:0006084':[],
 'GO:0006114':[],
 'GO:0006490':[],
 'GO:0006497':[],
 'GO:0006498':[],
 'GO:0006501':[],
 'GO:0006505':[],
 'GO:0006506':[],
 'GO:0006580':[],
 'GO:0006629':[],
 'GO:0006638':[],
 'GO:0006643':[],
 'GO:0006644':[],
 'GO:0006646':[],
 'GO:0006649':[],
 'GO:0006650':[],
 'GO:0006655':[],
 'GO:0006656':[],
 'GO:0006658':[],
 'GO:0006659':[],
 'GO:0006660':[],
 'GO:0006661':[],
 'GO:0006664':[],
 'GO:0006665':[],
 'GO:0006670':[],
 'GO:0006671':[],
 'GO:0006673':[],
 'GO:0006675':[],
 'GO:0006676':[],
 'GO:0006678':[],
 'GO:0006684':[],
 'GO:0006687':[],
 'GO:0006688':[],
 'GO:0006743':[],
 'GO:0006744':[],
 'GO:0006869':[],
 'GO:0006910':[],
 'GO:0006982':[],
 'GO:0007006':[],
 'GO:0007086':[],
 'GO:0007107':[],
 'GO:0008203':[],
 'GO:0008250':[],
 'GO:0008289':[],
 'GO:0008373':[],
 'GO:0008378':[],
 'GO:0008392':[],
 'GO:0008404':[],
 'GO:0008405':[],
 'GO:0008417':[],
 'GO:0008429':[],
 'GO:0008521':[],
 'GO:0008525':[],
 'GO:0008526':[],
 'GO:0008610':[],
 'GO:0008611':[],
 'GO:0008653':[],
 'GO:0008654':[],
 'GO:0008713':[],
 'GO:0008759':[],
 'GO:0008779':[],
 'GO:0008825':[],
 'GO:0008915':[],
 'GO:0008951':[],
 'GO:0009029':[],
 'GO:0009245':[],
 'GO:0009247':[],
 'GO:0009279':[],
 'GO:0009395':[],
 'GO:0009526':[],
 'GO:0009527':[],
 'GO:0009528':[],
 'GO:0009529':[],
 'GO:0009568':[],
 'GO:0009569':[],
 'GO:0009668':[],
 'GO:0009705':[],
 'GO:0009706':[],
 'GO:0009707':[],
 'GO:0009941':[],
 'GO:0010008':[],
 'GO:0010236':[],
 'GO:0010287':[],
 'GO:0010344':[],
 'GO:0010517':[],
 'GO:0010518':[],
 'GO:0010519':[],
 'GO:0010742':[],
 'GO:0010743':[],
 'GO:0010744':[],
 'GO:0010745':[],
 'GO:0010795':[],
 'GO:0010872':[],
 'GO:0010873':[],
 'GO:0010876':[],
 'GO:0010877':[],
 'GO:0010883':[],
 'GO:0010884':[],
 'GO:0010888':[],
 'GO:0010899':[],
 'GO:0010900':[],
 'GO:0010901':[],
 'GO:0010902':[],
 'GO:0010903':[],
 'GO:0012506':[],
 'GO:0012507':[],
 'GO:0012508':[],
 'GO:0012509':[],
 'GO:0012510':[],
 'GO:0012511':[],
 'GO:0014045':[],
 'GO:0014803':[],
 'GO:0014804':[],
 'GO:0015161':[],
 'GO:0015168':[],
 'GO:0015220':[],
 'GO:0015221':[],
 'GO:0015247':[],
 'GO:0015269':[],
 'GO:0015648':[],
 'GO:0015713':[],
 'GO:0015726':[],
 'GO:0015733':[],
 'GO:0015735':[],
 'GO:0015736':[],
 'GO:0015737':[],
 'GO:0015738':[],
 'GO:0015745':[],
 'GO:0015749':[],
 'GO:0015750':[],
 'GO:0015751':[],
 'GO:0015752':[],
 'GO:0015753':[],
 'GO:0015754':[],
 'GO:0015756':[],
 'GO:0015757':[],
 'GO:0015761':[],
 'GO:0015762':[],
 'GO:0015792':[],
 'GO:0015793':[],
 'GO:0015810':[],
 'GO:0015836':[],
 'GO:0015871':[],
 'GO:0015876':[],
 'GO:0015882':[],
 'GO:0015905':[],
 'GO:0015914':[],
 'GO:0015917':[],
 'GO:0015920':[],
 'GO:0016004':[],
 'GO:0016020':[],
 'GO:0016042':[],
 'GO:0016254':[],
 'GO:0016298':[],
 'GO:0016320':[],
 'GO:0017089':[],
 'GO:0017121':[],
 'GO:0017128':[],
 'GO:0018281':[],
 'GO:0019031':[],
 'GO:0019196':[],
 'GO:0019216':[],
 'GO:0019374':[],
 'GO:0019375':[],
 'GO:0019376':[],
 'GO:0019377':[],
 'GO:0019563':[],
 'GO:0019695':[],
 'GO:0019817':[],
 'GO:0019866':[],
 'GO:0019882':[],
 'GO:0019883':[],
 'GO:0019884':[],
 'GO:0019915':[],
 'GO:0019930':[],
 'GO:0020005':[],
 'GO:0022010':[],
 'GO:0022011':[],
 'GO:0030114':[],
 'GO:0030148':[],
 'GO:0030149':[],
 'GO:0030228':[],
 'GO:0030258':[],
 'GO:0030259':[],
 'GO:0030290':[],
 'GO:0030658':[],
 'GO:0030659':[],
 'GO:0030660':[],
 'GO:0030661':[],
 'GO:0030662':[],
 'GO:0030663':[],
 'GO:0030665':[],
 'GO:0030666':[],
 'GO:0030667':[],
 'GO:0030668':[],
 'GO:0030669':[],
 'GO:0030670':[],
 'GO:0030671':[],
 'GO:0030672':[],
 'GO:0030673':[],
 'GO:0030674':[],
 'GO:0030867':[],
 'GO:0030868':[],
 'GO:0030882':[],
 'GO:0030883':[],
 'GO:0030884':[],
 'GO:0030936':[],
 'GO:0031088':[],
 'GO:0031090':[],
 'GO:0031092':[],
 'GO:0031095':[],
 'GO:0031161':[],
 'GO:0031164':[],
 'GO:0031210':[],
 'GO:0031225':[],
 'GO:0031362':[],
 'GO:0031579':[],
 'GO:0031887':[],
 'GO:0031898':[],
 'GO:0031899':[],
 'GO:0031900':[],
 'GO:0031901':[],
 'GO:0031902':[],
 'GO:0031903':[],
 'GO:0031965':[],
 'GO:0031966':[],
 'GO:0031967':[],
 'GO:0031968':[],
 'GO:0031969':[],
 'GO:0031970':[],
 'GO:0031972':[],
 'GO:0031973':[],
 'GO:0031974':[],
 'GO:0031975':[],
 'GO:0032093':[],
 'GO:0032127':[],
 'GO:0032220':[],
 'GO:0032322':[],
 'GO:0032365':[],
 'GO:0032368':[],
 'GO:0032369':[],
 'GO:0032370':[],
 'GO:0032377':[],
 'GO:0032378':[],
 'GO:0032379':[],
 'GO:0032398':[],
 'GO:0032578':[],
 'GO:0032580':[],
 'GO:0032585':[],
 'GO:0032586':[],
 'GO:0032588':[],
 'GO:0032594':[],
 'GO:0032595':[],
 'GO:0032596':[],
 'GO:0032599':[],
 'GO:0032810':[],
 'GO:0032865':[],
 'GO:0032978':[],
 'GO:0032987':[],
 'GO:0032991':[],
 'GO:0032994':[],
 'GO:0033016':[],
 'GO:0033017':[],
 'GO:0033096':[],
 'GO:0033097':[],
 'GO:0033098':[],
 'GO:0033102':[],
 'GO:0033105':[],
 'GO:0033106':[],
 'GO:0033110':[],
 'GO:0033111':[],
 'GO:0033112':[],
 'GO:0033113':[],
 'GO:0033115':[],
 'GO:0033116':[],
 'GO:0033118':[],
 'GO:0033162':[],
 'GO:0033163':[],
 'GO:0033164':[],
 'GO:0033545':[],
 'GO:0033548':[],
 'GO:0033556':[],
 'GO:0033606':[],
 'GO:0033644':[],
 'GO:0033648':[],
 'GO:0033691':[],
 'GO:0033700':[],
 'GO:0033964':[],
 'GO:0033993':[],
 'GO:0034040':[],
 'GO:0034202':[],
 'GO:0034203':[],
 'GO:0034204':[],
 'GO:0034228':[],
 'GO:0034229':[],
 'GO:0034235':[],
 'GO:0034358':[],
 'GO:0034359':[],
 'GO:0034363':[],
 'GO:0034364':[],
 'GO:0034365':[],
 'GO:0034368':[],
 'GO:0034369':[],
 'GO:0034370':[],
 'GO:0034371':[],
 'GO:0034372':[],
 'GO:0034373':[],
 'GO:0034374':[],
 'GO:0034375':[],
 'GO:0034376':[],
 'GO:0034377':[],
 'GO:0034378':[],
 'GO:0034379':[],
 'GO:0034380':[],
 'GO:0034385':[],
 'GO:0034389':[],
 'GO:0034425':[],
 'GO:0034426':[],
 'GO:0034430':[],
 'GO:0034433':[],
 'GO:0034434':[],
 'GO:0034435':[],
 'GO:0034439':[],
 'GO:0034440':[],
 'GO:0034441':[],
 'GO:0034478':[],
 'GO:0034638':[],
 'GO:0034646':[],
 'GO:0034702':[],
 'GO:0035014':[],
 'GO:0035091':[],
 'GO:0035103':[],
 'GO:0035339':[],
 'GO:0035348':[],
 'GO:0035478':[],
 'GO:0035577':[],
 'GO:0035579':[],
 'GO:0035621':[],
 'GO:0035627':[],
 'GO:0035727':[],
 'GO:0036012':[],
 'GO:0036013':[],
 'GO:0036014':[],
 'GO:0036020':[],
 'GO:0036103':[],
 'GO:0036104':[],
 'GO:0036186':[],
 'GO:0036313':[],
 'GO:0036338':[],
 'GO:0036362':[],
 'GO:0036405':[],
 'GO:0036407':[],
 'GO:0036504':[],
 'GO:0038036':[],
 'GO:0039641':[],
 'GO:0039661':[],
 'GO:0039662':[],
 'GO:0042081':[],
 'GO:0042082':[],
 'GO:0042157':[],
 'GO:0042158':[],
 'GO:0042159':[],
 'GO:0042160':[],
 'GO:0042161':[],
 'GO:0042170':[],
 'GO:0042281':[],
 'GO:0042283':[],
 'GO:0042284':[],
 'GO:0042285':[],
 'GO:0042425':[],
 'GO:0042426':[],
 'GO:0042577':[],
 'GO:0042584':[],
 'GO:0042589':[],
 'GO:0042599':[],
 'GO:0042611':[],
 'GO:0042627':[],
 'GO:0042717':[],
 'GO:0042869':[],
 'GO:0042870':[],
 'GO:0042873':[],
 'GO:0042874':[],
 'GO:0042875':[],
 'GO:0042882':[],
 'GO:0042892':[],
 'GO:0042899':[],
 'GO:0042920':[],
 'GO:0042953':[],
 'GO:0043036':[],
 'GO:0043208':[],
 'GO:0043227':[],
 'GO:0043228':[],
 'GO:0043231':[],
 'GO:0043232':[],
 'GO:0043233':[],
 'GO:0043264':[],
 'GO:0043495':[],
 'GO:0043548':[],
 'GO:0043550':[],
 'GO:0043551':[],
 'GO:0043592':[],
 'GO:0043654':[],
 'GO:0043808':[],
 'GO:0043810':[],
 'GO:0043838':[],
 'GO:0043839':[],
 'GO:0043842':[],
 'GO:0044162':[],
 'GO:0044167':[],
 'GO:0044169':[],
 'GO:0044171':[],
 'GO:0044173':[],
 'GO:0044175':[],
 'GO:0044178':[],
 'GO:0044185':[],
 'GO:0044186':[],
 'GO:0044188':[],
 'GO:0044190':[],
 'GO:0044191':[],
 'GO:0044192':[],
 'GO:0044193':[],
 'GO:0044199':[],
 'GO:0044200':[],
 'GO:0044201':[],
 'GO:0044202':[],
 'GO:0044232':[],
 'GO:0044233':[],
 'GO:0044241':[],
 'GO:0044242':[],
 'GO:0044255':[],
 'GO:0044258':[],
 'GO:0044385':[],
 'GO:0044386':[],
 'GO:0045017':[],
 'GO:0045026':[],
 'GO:0045052':[],
 'GO:0045054':[],
 'GO:0045117':[],
 'GO:0045121':[],
 'GO:0045125':[],
 'GO:0045332':[],
 'GO:0045833':[],
 'GO:0045834':[],
 'GO:0046027':[],
 'GO:0046335':[],
 'GO:0046336':[],
 'GO:0046337':[],
 'GO:0046338':[],
 'GO:0046341':[],
 'GO:0046346':[],
 'GO:0046347':[],
 'GO:0046411':[],
 'GO:0046433':[],
 'GO:0046460':[],
 'GO:0046461':[],
 'GO:0046465':[],
 'GO:0046466':[],
 'GO:0046467':[],
 'GO:0046468':[],
 'GO:0046470':[],
 'GO:0046471':[],
 'GO:0046474':[],
 'GO:0046475':[],
 'GO:0046479':[],
 'GO:0046480':[],
 'GO:0046485':[],
 'GO:0046486':[],
 'GO:0046488':[],
 'GO:0046493':[],
 'GO:0046503':[],
 'GO:0046505':[],
 'GO:0046506':[],
 'GO:0046512':[],
 'GO:0046527':[],
 'GO:0046623':[],
 'GO:0046624':[],
 'GO:0046625':[],
 'GO:0046658':[],
 'GO:0046834':[],
 'GO:0046836':[],
 'GO:0046839':[],
 'GO:0046858':[],
 'GO:0046859':[],
 'GO:0046860':[],
 'GO:0046861':[],
 'GO:0046862':[],
 'GO:0046889':[],
 'GO:0046890':[],
 'GO:0046981':[],
 'GO:0046989':[],
 'GO:0047066':[],
 'GO:0047157':[],
 'GO:0047177':[],
 'GO:0047178':[],
 'GO:0047179':[],
 'GO:0047600':[],
 'GO:0047909':[],
 'GO:0048003':[],
 'GO:0048006':[],
 'GO:0048007':[],
 'GO:0048017':[],
 'GO:0048210':[],
 'GO:0048269':[],
 'GO:0048279':[],
 'GO:0048280':[],
 'GO:0048288':[],
 'GO:0048475':[],
 'GO:0050494':[],
 'GO:0050738':[],
 'GO:0050746':[],
 'GO:0050747':[],
 'GO:0050748':[],
 'GO:0050872':[],
 'GO:0050994':[],
 'GO:0050995':[],
 'GO:0050996':[],
 'GO:0051055':[],
 'GO:0051267':[],
 'GO:0051377':[],
 'GO:0051469':[],
 'GO:0051665':[],
 'GO:0051697':[],
 'GO:0051861':[],
 'GO:0051872':[],
 'GO:0051977':[],
 'GO:0051978':[],
 'GO:0051999':[],
 'GO:0052631':[],
 'GO:0052712':[],
 'GO:0055035':[],
 'GO:0055036':[],
 'GO:0055037':[],
 'GO:0055038':[],
 'GO:0055085':[],
 'GO:0055088':[],
 'GO:0055091':[],
 'GO:0055102':[],
 'GO:0060105':[],
 'GO:0060106':[],
 'GO:0060191':[],
 'GO:0060192':[],
 'GO:0060193':[],
 'GO:0060198':[],
 'GO:0060199':[],
 'GO:0060200':[],
 'GO:0060201':[],
 'GO:0060203':[],
 'GO:0060229':[],
 'GO:0060230':[],
 'GO:0060510':[],
 'GO:0060587':[],
 'GO:0060588':[],
 'GO:0060642':[],
 'GO:0060696':[],
 'GO:0060697':[],
 'GO:0060856':[],
 'GO:0060857':[],
 'GO:0060987':[],
 'GO:0060988':[],
 'GO:0060989':[],
 'GO:0060990':[],
 'GO:0061007':[],
 'GO:0061008':[],
 'GO:0061024':[],
 'GO:0061025':[],
 'GO:0061091':[],
 'GO:0061092':[],
 'GO:0061093':[],
 'GO:0061200':[],
 'GO:0061202':[],
 'GO:0061474':[],
 'GO:0061588':[],
 'GO:0061592':[],
 'GO:0061701':[],
 'GO:0061724':[],
 'GO:0061725':[],
 'GO:0061739':[],
 'GO:0061796':[],
 'GO:0062040':[],
 'GO:0062242':[],
 'GO:0062243':[],
 'GO:0062245':[],
 'GO:0065005':[],
 'GO:0065010':[],
 'GO:0070075':[],
 'GO:0070081':[],
 'GO:0070083':[],
 'GO:0070088':[],
 'GO:0070089':[],
 'GO:0070112':[],
 'GO:0070113':[],
 'GO:0070114':[],
 'GO:0070115':[],
 'GO:0070118':[],
 'GO:0070258':[],
 'GO:0070268':[],
 'GO:0070381':[],
 'GO:0070391':[],
 'GO:0070392':[],
 'GO:0070394':[],
 'GO:0070395':[],
 'GO:0070396':[],
 'GO:0070505':[],
 'GO:0070782':[],
 'GO:0070821':[],
 'GO:0070915':[],
 'GO:0071071':[],
 'GO:0071072':[],
 'GO:0071073':[],
 'GO:0071129':[],
 'GO:0071130':[],
 'GO:0071210':[],
 'GO:0071223':[],
 'GO:0071396':[],
 'GO:0071449':[],
 'GO:0071617':[],
 'GO:0071723':[],
 'GO:0071813':[],
 'GO:0071814':[],
 'GO:0071825':[],
 'GO:0071827':[],
 'GO:0071967':[],
 'GO:0071968':[],
 'GO:0072492':[],
 'GO:0072562':[],
 'GO:0072564':[],
 'GO:0075513':[],
 'GO:0080064':[],
 'GO:0080065':[],
 'GO:0080177':[],
 'GO:0085017':[],
 'GO:0085019':[],
 'GO:0085026':[],
 'GO:0090077':[],
 'GO:0090078':[],
 'GO:0090107':[],
 'GO:0090108':[],
 'GO:0090153':[],
 'GO:0090154':[],
 'GO:0090155':[],
 'GO:0090156':[],
 'GO:0090159':[],
 'GO:0090174':[],
 'GO:0090218':[],
 'GO:0090219':[],
 'GO:0090318':[],
 'GO:0090319':[],
 'GO:0090520':[],
 'GO:0090522':[],
 'GO:0097001':[],
 'GO:0097003':[],
 'GO:0097004':[],
 'GO:0097005':[],
 'GO:0097035':[],
 'GO:0097040':[],
 'GO:0097045':[],
 'GO:0097092':[],
 'GO:0097093':[],
 'GO:0097209':[],
 'GO:0097212':[],
 'GO:0097232':[],
 'GO:0097233':[],
 'GO:0097234':[],
 'GO:0097302':[],
 'GO:0097303':[],
 'GO:0097304':[],
 'GO:0097348':[],
 'GO:0097381':[],
 'GO:0097384':[],
 'GO:0097478':[],
 'GO:0097488':[],
 'GO:0097598':[],
 'GO:0097691':[],
 'GO:0097707':[],
 'GO:0098585':[],
 'GO:0098588':[],
 'GO:0098852':[],
 'GO:0098856':[],
 'GO:0098857':[],
 'GO:0098895':[],
 'GO:0098896':[],
 'GO:0098897':[],
 'GO:0098920':[],
 'GO:0098944':[],
 'GO:0098954':[],
 'GO:0098993':[],
 'GO:0099012':[],
 'GO:0099022':[],
 'GO:0099025':[],
 'GO:0099026':[],
 'GO:0099027':[],
 'GO:0099028':[],
 'GO:0099029':[],
 'GO:0099030':[],
 'GO:0099031':[],
 'GO:0099032':[],
 'GO:0099033':[],
 'GO:0099034':[],
 'GO:0099035':[],
 'GO:0099036':[],
 'GO:0099037':[],
 'GO:0099039':[],
 'GO:0099041':[],
 'GO:0099044':[],
 'GO:0099144':[],
 'GO:0099169':[],
 'GO:0099501':[],
 'GO:0099541':[],
 'GO:0099552':[],
 'GO:0101003':[],
 'GO:0101004':[],
 'GO:0102003':[],
 'GO:0102070':[],
 'GO:0102431':[],
 'GO:0102657':[],
 'GO:0102771':[],
 'GO:0102772':[],
 'GO:0102834':[],
 'GO:0102865':[],
 'GO:0102993':[],
 'GO:0103003':[],
 'GO:0103015':[],
 'GO:0103028':[],
 'GO:0106006':[],
 'GO:0106073':[],
 'GO:0106125':[],
 'GO:0106126':[],
 'GO:0106175':[],
 'GO:0106236':[],
 'GO:0106254':[],
 'GO:0106402':[],
 'GO:0110112':[],
 'GO:0110113':[],
 'GO:0110114':[],
 'GO:0110146':[],
 'GO:0120009':[],
 'GO:0120010':[],
 'GO:0120011':[],
 'GO:0120012':[],
 'GO:0120013':[],
 'GO:0120014':[],
 'GO:0120015':[],
 'GO:0120016':[],
 'GO:0120017':[],
 'GO:0120019':[],
 'GO:0120020':[],
 'GO:0120021':[],
 'GO:0120146':[],
 'GO:0120149':[],
 'GO:0120201':[],
 'GO:0120202':[],
 'GO:0120281':[],
 'GO:0120322':[],
 'GO:0120323':[],
 'GO:0140025':[],
 'GO:0140042':[],
 'GO:0140043':[],
 'GO:0140282':[],
 'GO:0140303':[],
 'GO:0140326':[],
 'GO:0140327':[],
 'GO:0140328':[],
 'GO:0140329':[],
 'GO:0140331':[],
 'GO:0140333':[],
 'GO:0140353':[],
 'GO:0140354':[],
 'GO:0140366':[],
 'GO:0140443':[],
 'GO:0140444':[],
 'GO:0140474':[],
 'GO:0140504':[],
 'GO:0140505':[],
 'GO:0140513':[],
 'GO:0140522':[],
 'GO:0140523':[],
 'GO:0140567':[],
 'GO:0140622':[],
 'GO:0140735':[],
 'GO:1900130':[],
 'GO:1900131':[],
 'GO:1900132':[],
 'GO:1900161':[],
 'GO:1900162':[],
 'GO:1900163':[],
 'GO:1901373':[],
 'GO:1901759':[],
 'GO:1901760':[],
 'GO:1902068':[],
 'GO:1902069':[],
 'GO:1902070':[],
 'GO:1902300':[],
 'GO:1902388':[],
 'GO:1902994':[],
 'GO:1902995':[],
 'GO:1902999':[],
 'GO:1903000':[],
 'GO:1903001':[],
 'GO:1903002':[],
 'GO:1903059':[],
 'GO:1903060':[],
 'GO:1903061':[],
 'GO:1903725':[],
 'GO:1903726':[],
 'GO:1903727':[],
 'GO:1904121':[],
 'GO:1904729':[],
 'GO:1904730':[],
 'GO:1904731':[],
 'GO:1905038':[],
 'GO:1905329':[],
 'GO:1905691':[],
 'GO:1905952':[],
 'GO:1905953':[],
 'GO:1905954':[],
 'GO:1990044':[],
 'GO:1990050':[],
 'GO:1990052':[],
 'GO:1990199':[],
 'GO:1990379':[],
 'GO:1990455':[],
 'GO:1990458':[],
 'GO:1990482':[],
 'GO:1990530':[],
 'GO:1990531':[],
 'GO:1990668':[],
 'GO:1990669':[],
 'GO:1990670':[],
 'GO:1990671':[],
 'GO:1990672':[],
 'GO:1990674':[],
 'GO:1990675':[],
 'GO:1990676':[],
 'GO:1990684':[],
 'GO:1990685':[],
 'GO:1990686':[],
 'GO:1990687':[],
 'GO:1990688':[],
 'GO:1990689':[],
 'GO:1990690':[],
 'GO:1990691':[],
 'GO:1990692':[],
 'GO:1990725':[],
 'GO:1990777':[],
 'GO:1990816':[],
 'GO:1990836':[],
 'GO:2001138':[],
 'GO:2001139':[],
 'GO:2001140':[],
 'GO:2001289':[]
}


# In[21]:


#Create the function with def and have it take 2 inputs(dataframe with Symbols and GO terms and list of GO Terms)
#iterate through the list of possible GO Terms (921 iterations)
#inside the above for loop, iterate through the dataframe column "GOTERM" by index (iterations equal number of df rows)
#if the current GO term matches the index iterated GO term in the df column, execute the code below
#using the current GO term as a key, add the index matched row of the symbol column to the dictionary
#This function will return the filled in dictionary when called
def GeneCollector(dfinput, golist):
    for term in golist:
        for i in range(len(dfinput["GO TERM"])):
            if term == dfinput["GO TERM"][i]:
                GOTermGeneAssociator[term].append(dfinput["SYMBOL"][i])  
    return(GOTermGeneAssociator)


# In[22]:


First600= GeneCollector(df1, lipid_nodes)


# In[23]:


Full921 =GeneCollector(df2, lipid_nodes)


# In[24]:


Full921


# In[25]:


#create an empty dictionary to hold unique gene values
UniqueGenestoGOTerms={}

#loop through keys in filled in dictionary 
#run set() on the values to acquire unique genes
#store current key and unique genes in the empty dictionary above as a key:value pair
for key in GOTermGeneAssociator:
    UniqueGenestoGOTerms[key]=set(GOTermGeneAssociator[key])

#print the final dictionary
print(UniqueGenestoGOTerms)


# In[26]:


#Creates an empty list called Ultimate list. If UniqueGenestoGOTerms has any contents
#for a given GO term, the if statement executes. A dictionary is created then converted to house the GO Term ID and Name.
#For a given GO Term, the associated genes and the specific respective DGEstats are stored in a df called MappedGenes.
#output concatenates TERM IDs/Names/DGEstats dataframes and is added to the ultimatelist.
#Result concatenates the contents to make a single final comprehensive df.

Ultimate_listv3=[]

for GOT in list(lipid_nodes):
    if len(list(UniqueGenestoGOTerms[GOT]))>0:
        #create a new df to append the df_final df onto
        GO_dict = {'GO TERM': [GOT], 'GO NAME': [id_to_name[GOT]]}
        GO_TERMNAME = pd.DataFrame(data = GO_dict)
        MappedGenes = df_final[df_final['Gene'].isin(list(UniqueGenestoGOTerms[GOT]))]
        output = pd.concat([GO_TERMNAME, MappedGenes], ignore_index=True)
        Ultimate_listv3.append(output)

Result= pd.concat(Ultimate_listv3)


# In[27]:


#writes Result to a csv dropping the index column
Result.to_csv("Full920MappedDGEstats.csv",encoding='utf-8', index=False)


# In[28]:


#Hierarchy list holds the known relationships of GOTs of interest
#Calls NetworkMapper on each GOT of interest.
Hierarchy_list=[]

for GOT in list(lipid_nodes):                    
        Hierarchy_list.append(NetworkMapper(GOT))

print(Hierarchy_list)


# In[29]:


dfcolumnnames= [str(GOT) for GOT in lipid_nodes]


# In[30]:


Hierarchy_listdf= pd.DataFrame(Hierarchy_list)
Hierarchy_listdf= (np.transpose(Hierarchy_listdf))


# In[31]:


Hierarchy_listdf


# In[32]:


#writes Result to a csv dropping the index column
#Hierarchy_listdf.to_csv("Full920Relationships.csv",encoding='utf-8', index=True)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[33]:


#def GetMaxLenV(dfToBe):
#    return max((len(v), k) for k,v in dfToBe.items())


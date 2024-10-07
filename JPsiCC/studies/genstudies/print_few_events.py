'''
Print few events.
'''

from kytools import jsonreader
import os
import ROOT
import pandas as pd
from particle import Particle as P

# Disable multithreading for now
ROOT.DisableImplicitMT()

# Get TChain of generator-level events in the NANOAOD format
genevents = jsonreader.get_chain_from_json(anpath='JPsiCC',
                                           jsonname='NANOAOD.json',
                                           keys=['mariadlf', '20240321', 'GEN-signal'],
                                           treename='Events')

# Create RDF from TChain
rdf = ROOT.RDataFrame(genevents)

# Make Numpy array of select columns
arr_pdgid  = rdf.AsNumpy(['GenPart_pdgId'])['GenPart_pdgId']
arr_idxmom = rdf.AsNumpy(['GenPart_genPartIdxMother'])['GenPart_genPartIdxMother']

# Make Pandas DataFrame
pddf = pd.DataFrame(data=(arr_pdgid, arr_idxmom))

# Print few events
for row in range(0, 3):
    pdgid = arr_pdgid[row]
    idxmom = arr_idxmom[row]

    for i in range(0, len(idxmom)):
        if idxmom[i] < 0: continue
        if i <= 1:
            #continue
            print (P.from_pdgid(pdgid[i]))
        elif idxmom[i] <= 1:
            print (P.from_pdgid(pdgid[i]), P.from_pdgid(pdgid[idxmom[i]]))
            print(i, idxmom[i])
        elif idxmom[idxmom[i]] <= 1:
            print (P.from_pdgid(pdgid[i]), P.from_pdgid(pdgid[idxmom[i]]), P.from_pdgid(pdgid[idxmom[idxmom[i]]]))
            print(i, idxmom[i], idxmom[idxmom[i]])
        elif idxmom[idxmom[idxmom[i]]] <= 1:
            print (P.from_pdgid(pdgid[i]), P.from_pdgid(pdgid[idxmom[i]]), P.from_pdgid(pdgid[idxmom[idxmom[i]]]), P.from_pdgid(pdgid[idxmom[idxmom[idxmom[i]]]]))
            print(i, idxmom[i], idxmom[idxmom[i]], idxmom[idxmom[idxmom[i]]])
        else:
            print (P.from_pdgid(pdgid[i]), P.from_pdgid(pdgid[idxmom[i]]), P.from_pdgid(pdgid[idxmom[idxmom[i]]]), P.from_pdgid(pdgid[idxmom[idxmom[idxmom[i]]]]), P.from_pdgid(pdgid[idxmom[idxmom[idxmom[idxmom[i]]]]]))
            print(i, idxmom[i], idxmom[idxmom[i]], idxmom[idxmom[idxmom[i]]], idxmom[idxmom[idxmom[idxmom[i]]]])
    print ('--------------------------------')
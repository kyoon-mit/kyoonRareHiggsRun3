import json
import os
from ROOT import TChain

def get_chain_from_json(anpath, jsonname, keys, treename='events'):
    """Opens ROOT files from a JSON file and returns a TChain.

    Args:
        anpath (str): Path to the analysis folder, relative to HRARE_DIR.
            e.g. 'JPsiCC'
        jsonname (str): Name of the JSON file.
        keys (list(str)): List of the keys to the JSON object, in sequential order.
        treename (str): 
            Defaults to 'events'.

    Returns:
        chain (ROOT.TChain): TChain of the ROOT files.

    Raises:
        TypeError: If keys object is not list(str).
        TypeError: If json object is not list(str).
    """
    jsonfile = json.load(open(os.path.join(os.environ['HRARE_DIR'], anpath, 'json', jsonname)))
    jsonobject = jsonfile
    for i in range(len(keys)):
        jsonobject = jsonobject[keys[i]]
    if not all(isinstance(item, str) for item in keys):
        raise TypeError('Argument provided in keys is not list(str).')
    if not all(isinstance(item, str) for item in jsonobject):
        raise TypeError('JSON object is not list(str).')
    chain = TChain(treename)
    for item in jsonobject:
        chain.Add(item)
    return chain


def get_chain_from_json_xrtd(anpath, jsonname, keys, treename='events', xrtd_proxy='root://xrootd.cmsaf.mit.edu//',):
    """Opens ROOT files from a JSON file and returns a TChain.

    Uses XRootD to open the files.

    Args:
        anpath (str): Path to the analysis folder, relative to HRARE_DIR.
            e.g. 'JPsiCC'
        jsonname (str): Name of the JSON file.
        keys (list(str)): List of the keys to the JSON object, in sequential order.
        treename (str): 
            Defaults to 'events'.
        xrtd_proxy (str):
            Name of the XRootD proxy. Defaults to root://xrootd.cmsaf.mit.edu'.

    Returns:
        chain (ROOT.TChain): TChain of the ROOT files.

    Raises:
        TypeError: If keys object is not list(str).
        TypeError: If json object is not list(str).
    """
    jsonfile = json.load(open(os.path.join(os.environ['HRARE_DIR'], anpath, 'json', jsonname)))
    jsonobject = jsonfile
    for i in range(len(keys)):
        jsonobject = jsonobject[keys[i]]
    if not all(isinstance(item, str) for item in keys):
        raise TypeError('Argument provided in keys is not list(str).')
    if not all(isinstance(item, str) for item in jsonobject):
        raise TypeError('JSON object is not list(str).')
    chain = TChain(treename)
    print(f'{'\033[1;36m'}kytools: Opening files from json:{'\033[0m'}')
    for item in jsonobject:
        print('    ' + os.path.join(xrtd_proxy, *item.split('/')))
        chain.Add(os.path.join(xrtd_proxy, *item.split('/')))
    return chain
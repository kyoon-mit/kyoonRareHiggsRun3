import json
import os
from ROOT import TChain, RDF

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
        xrtd_proxy (str): Name of the XRootD proxy.
            Defaults to root://xrootd.cmsaf.mit.edu'.

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
    print('{}kytools: Opening files from json...{}'.format('\033[1;36m', '\033[0m'))
    for item in jsonobject:
        print('    ' + os.path.join(xrtd_proxy, *item.split('/')))
        chain.Add(os.path.join(xrtd_proxy, *item.split('/')))
    return chain

def get_object_from_json(anpath, jsonname, keys):
    """Retrieve an object from a JSON file.

    Args:
        anpath (str): Path to the analysis folder, relative to HRARE_DIR.
            e.g. 'JPsiCC'
        jsonname (str): Name of the JSON file.
        keys (list(str)): List of the keys to the JSON object, in sequential order.

    Returns:
        jsonobject (list or dict): Retrieved object.

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
    return jsonobject

def get_rdf_from_json_spec(anpath, jsonname, calc_nfiles=False):
    """Retrieve an RDataFrame from JSON spec file.

    Args:
        anpath (str): Path to the analysis folder, relative to HRARE_DIR.
            e.g. 'JPsiCC'
        jsonname (str): Name of the JSON file.
        calc_nfiles (bool, optional): Calculate the number of files in the JSON spec file.
            Defaults to False.

    Returns:
        rdf (ROOT.RDataFrame): Retrieved ROOT.RDataFrame object.
    """
    jsonfile = os.path.join(os.environ['HRARE_DIR'], anpath, 'json', jsonname)
    if calc_nfiles:
        nfiles_total = 0
        samples = get_object_from_json(anpath, jsonname, ['samples'])
        print_detailed = '..... SAMPLE              | NUMBER OF FILES'
        for key_sample in samples.keys():
            nfiles_sample = len(samples[key_sample]['files'])
            nfiles_total += nfiles_sample
            print_detailed += f'\n..... {'_'.join(key_sample.split()[:2]):19.19} | {nfiles_sample}'
        print(f'INFO: the total number of files in {jsonname} is {nfiles_total}.')
        print(print_detailed)
    rdf = RDF.Experimental.FromSpec(jsonfile)
    return rdf
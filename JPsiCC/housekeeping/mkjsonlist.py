'''
This housekeeping script prints a string output of a list of the file names in a directory,
which can be copied and pasted onto a json file.
'''

from kytools import jsonreader
from subprocess import check_output
from glob import glob
import os, sys, re, json

# mydir = '/data/submit/cms/store/user/mariadlf/nano/GluGluH_HJPsiCC/NANOAOD_test3/'
mydir = '/store/user/paus/nanohr/D04/'

def sorted_alphanumeric(data):
    ''' Sorts a list so that [1, 10, 2] is properly sorted [1, 2, 10].
    https://stackoverflow.com/questions/4813061/non-alphanumeric-list-order-from-os-listdir
    '''
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)

def find_directory(mydir, use_xrootd=True, xrtd_proxy='root://xrootd.cmsaf.mit.edu', f_endswith='.root'):
    dirlist = []
    if use_xrootd is True:
        if '/data/submit/cms' in mydir:
            mydir = mydir.replace('/data/submit/cms', '')
        xrtd_files = check_output(['xrdfs', xrtd_proxy, 'ls', mydir]).decode(sys.stdout.encoding)
        for item in xrtd_files.split():
            if ('failed/' in item) or ('log/' in item) or ('.txt' in item): continue
            elif (item.endswith(f_endswith)):
                dirlist.append(item)
    else:
        dirlist = sorted_alphanumeric(glob(os.path.join(mydir, f'*{f_endswith}')))
    return dirlist

def find_all_directories(topdir, use_xrootd=True, xrtd_proxy='root://xrootd.cmsaf.mit.edu'):
    dirdict = {}
    dirlist = find_directory(topdir, use_xrootd, xrtd_proxy, f_endswith='')
    for subdir in dirlist:
        print(subdir)
        files = find_directory(subdir, use_xrootd, xrtd_proxy, f_endswith='.root')
        dirdict[os.path.basename(os.path.normpath(subdir))] = files
    return dirdict

def make_json(dirs):
    json_str = json.dumps(dirs, indent=4)
    return json_str

def make_json_spec(anpath, jsonname, keys, year, data_type):
    '''Make the json spec file that will be sued in RDataFrame.

    See the 'Creating an RDataFrame from a dataset specification file' section in
    https://root.cern/doc/master/classROOT_1_1RDataFrame.html#crash-course.

    Args:
        anpath (str): Path to the analysis folder, relative to HRARE_DIR.
            e.g. 'JPsiCC'
        jsonname (str): Name of the JSON file.
            It must be a file that includes the year, dataset name, cross section,
            and luminosity, e.g. 'MC_bkg_names.json'.
        keys (list(str)): List of the keys to the JSON object, in sequential order.
        year (int): Year of the dataset.
        data_type (str): Type of data, i.e. 'DATA', 'MC_BKG', or 'MC_SIG'.

    Returns:
        jsondict (dict): Dictionary of the JSON itmes.
    '''
    jsondict = {}
    jsondict['samples'] = {}
    dataset_meta = jsonreader.get_object_from_json(anpath, jsonname, keys)
    nanoaod_json = jsonreader.get_object_from_json(anpath, 'NANOAOD.json', keys) # where the actual file is stored
    for item in dataset_meta[year]:
        jsondict['samples'][item['dataset']] = {
            'trees': ['Events'],
            'files': nanoaod_json[item['dataset']],
            'metadata': {
                'xsec': item['xsec'],
                'xsec_sigma': item['xsec_sigma'],
                'lumi': item['lumi'],
                'sample_category': data_type,
                'year': year
            }
        }
    return jsondict

if __name__=='__main__':
    dirdict = make_json_spec('JPsiCC', 'DATA_bkg_names.json', ['kraken', '202405'], 2018, 'DATA_BKG')
    json_str = make_json(dirdict)
    print(json_str)
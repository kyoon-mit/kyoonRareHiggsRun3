'''
This housekeeping script prints a string output of a list of the file names in a directory,
which can be copied and pasted onto a json file.
'''

from subprocess import check_output
from glob import glob
import os, sys, re, json

# mydir = '/data/submit/cms/store/user/mariadlf/nano/GluGluH_HJPsiCC/NANOAOD_test3/'
mydir = '/store/user/paus/nanohr/D04/Charmonium+Run2016B-ver1_HIPM_UL2016_MiniAODv2-v1+MINIAOD'

def sorted_alphanumeric(data):
    ''' Sorts a list so that [1, 10, 2] is properly sorted [1, 2, 10].
    https://stackoverflow.com/questions/4813061/non-alphanumeric-list-order-from-os-listdir
    '''
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)

def find_directory(mydir, use_xrootd=False, xrtd_proxy='root://xrootd.cmsaf.mit.edu/'):
    dirlist = []
    if use_xrootd is True:
        if '/data/submit/cms' in mydir:
            mydir = mydir.replace('/data/submit/cms', '')
        xrtd_files = check_output(['xrdfs', xrtd_proxy, 'ls', mydir]).decode(sys.stdout.encoding)
        for item in xrtd_files.split():
            if ('failed/' in item) or ('log/' in item) or ('.txt' in item): continue
            elif (item.endswith('.root')):
                dirlist.append(os.path.join(xrtd_proxy, item))
    else:
        dirlist = sorted_alphanumeric(glob(os.path.join(mydir, '*.root')))
    return dirlist

def make_json_list(dirlist):
    json_str = json.dumps(dirlist, indent=4)
    return json_str

if __name__=='__main__':
    json_str = make_json_list(find_directory(mydir, True))
    print(json_str)
'''
This housekeeping script prints a string output of a list of the file names in a directory,
which can be copied and pasted onto a json file.
'''

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

if __name__=='__main__':
    dirdict = find_all_directories(mydir)
    print(dirdict.items())
    json_str = make_json(dirdict)
    print(json_str)
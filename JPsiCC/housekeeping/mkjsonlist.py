'''
This housekeeping script prints a string output of a list of the file names in a directory,
which can be copied and pasted onto a json file.
'''

from glob import glob
import os
import re
import json

mydir = '/data/submit/cms/store/user/mariadlf/nano/GluGluH_HJPsiCC/NANOAOD/'

def sorted_alphanumeric(data):
    ''' Sorts a list so that [1, 10, 2] is properly sorted [1, 2, 10].
    https://stackoverflow.com/questions/4813061/non-alphanumeric-list-order-from-os-listdir
    '''
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)

dirlist = sorted_alphanumeric(glob(os.path.join(mydir, '*')))

json_str = json.dumps(dirlist, indent=4)
print(json_str)
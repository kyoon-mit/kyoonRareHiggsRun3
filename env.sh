#!/bin/bash

cmsenv
export HRARE_DIR=$(dirname $(realpath -s $BASH_SOURCE[0]))
python3 -m pip install $HRARE_DIR/python-package/.
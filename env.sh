#!/bin/bash

conda activate hrare
export HRARE_DIR=$(dirname $(realpath -s $BASH_SOURCE[0]))
pip install $HRARE_DIR/python-package/.

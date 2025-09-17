#test

import numpy as np
import pyreadr
import os
os.chdir("C:\\Users\\dan94\\rips_project\\")
os.listdir()

sepsis = pyreadr.read_r("sepsis_rna_array.rds")
sepsis = sepsis[None]

from builtins import object

__all__ = ["WriteASCII", "compile_ascii", "check_compiled", "compiled_exists"]

import os
import numpy as np
import subprocess
from powderday.nebular_emission.cloudy_tools import grouper
#getting cfg.par accessible outside the definitions (ala pd_front_end.py)
import powderday.config as cfg
import sys
#script, pardir, parfile, modelfile = sys.argv
print(sys.argv)
sys.path.insert(0, sys.argv[1])
par = __import__(sys.argv[2])
model = __import__(sys.argv[3])
cfg.par = par  # re-write cfg.par for all modules that read this in now
cfg.model = model

"""
--------------------------------------------------------------------------------------
Based on cloudyfsps written by Nell Byler.
(Source https://github.com/nell-byler/cloudyfsps/blob/master/cloudyfsps/ASCIItools.py 
retrieved in October 2019)
--------------------------------------------------------------------------------------
"""


if (cfg.par.add_neb_emission) and not cfg.par.use_cloudy_tables:

    try:
        CLOUDY_EXE = os.environ['CLOUDY_EXE']
    except KeyError:
        raise KeyError("[powderday/nebular_emission/ASCII_TOOLS:] You must have set a system environment CLOUDY_EXE")

    try:
        # taking the last in path, if more than one directory given
        CLOUDY_DATA_PATH = os.environ['CLOUDY_DATA_PATH'].split(':')[-1]
    except KeyError:
        print('Cloudy data path not set. Assuming standard cloudy structure')
        CLOUDY_DATA_PATH = '/'.join(CLOUDY_EXE.split('/')[:-2])+'/data'

class WriteASCII(object):
    """
    -------------------------------------------------------------------------
    Prints an input ascii files readable by CLOUDY
    Calling sequence:
        writeASCII('outfile.ascii', lam_arr, spec_arr, model_arr, **kwargs)
    Dictionary with header information - change any of these values by
    inputting them as kwargs.
    -------------------------------------------------------------------------
    """
    def __init__(self, outfile, lam, flu, **kwargs):
        self.nom_dict = {'nmod': 94, 'ndim': 1, 'npar': 1, 'nx': 1963,
                         'x': 'lambda', 'conv1': 1.0, 'peraa': False,
                         'conv2': 3.839e33, 'par1': 'age', 'par1_val': 1.e6}
        self.init_pars(**kwargs)
        self.file = open('/'.join([CLOUDY_DATA_PATH, outfile]), 'w')
        self.write_header()
        self.write_body(lam, flu)
        self.file.close()

    def init_pars(self, **kwargs):
        for key, value in list(kwargs.items()):
            self.nom_dict[key] = value
        if self.nom_dict['peraa']:
            self.nom_dict['f_type'] = 'F_lambda'
        else:
            self.nom_dict['f_type'] = 'F_nu'

    def write_header(self):
        """
        Header for cloudy ascii files
        """
        self.file.write("  20060612\n")
        self.file.write("  %i\n" % self.nom_dict['ndim'])
        self.file.write("  %i\n" % self.nom_dict['npar'])
        self.file.write("  %s\n" % self.nom_dict['par1'])        # first param
        self.file.write("  %i\n" % self.nom_dict['nmod'])        # total number of mods
        self.file.write("  %i\n" % self.nom_dict['nx'])          # number of lam pts
        self.file.write("  %s\n" % self.nom_dict['x'])           # lambda or freq
        self.file.write("  %.8e\n" % self.nom_dict['conv1'])     # AA or Hz
        self.file.write("  %s\n" % self.nom_dict['f_type'])      # F_lam or F_nu
        self.file.write("  %.8e\n" % self.nom_dict['conv2'])     # units
        self.file.write("  %.8e\n" % self.nom_dict['par1_val'])

    def write_data(self, array):
        """
        write array with 5 items per line in format 1.0000e+00
        """
        for chunk in grouper(5, array):
            self.file.write("  " + "  ".join("%1.7e" % x for x in chunk) + "\n")

    def write_body(self, lam, flu):
        self.write_data(lam)
        flu[(flu < 0.0)] = 0.0
        self.write_data(flu)


def compile_ascii(ascii_file):
    fname = format(ascii_file.split('.')[0])
    comp_file = CLOUDY_DATA_PATH + '/' + fname + '.in'
    f = open(comp_file, 'w')
    f.write('compile stars "{}"\n'.format(ascii_file))
    f.close()
    to_run = 'cd {} ; {} {}.in'.format(CLOUDY_DATA_PATH, CLOUDY_EXE, fname)
    proc = subprocess.Popen(to_run, shell=True)
    proc.communicate()


def check_compiled(ascii_file):
    """
    checks to make sure ascii_file.mod exists and that
    the words "Cloudy exited OK' are in compile.out
    """
    out_file = CLOUDY_DATA_PATH+'/' + format(ascii_file.split('.')[0]) + '.out'
    f = open(out_file, 'r')
    content = f.readlines()
    f.close()
    comp_mod = '{}/{}.mod'.format(CLOUDY_DATA_PATH, ascii_file.split('.')[0])
    check = np.all(['OK' in content[-1],
                    os.path.exists(comp_mod)])
    return check


def compiled_exists(filename):
    if filename.split('.')[-1] == 'mod':
        return os.path.exists('/'.join([CLOUDY_DATA_PATH, filename]))
    else:
        return os.path.exists('/'.join([CLOUDY_DATA_PATH, filename.split('.')[0]+'.mod']))

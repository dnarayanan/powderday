from powderday.nebular_emission.abund import getNebAbunds
from powderday.nebular_emission.ASCIItools import *
from astropy import constants
from powderday.nebular_emission.cloudy_tools import air_to_vac, calc_LogU
import logging
import numpy as np
import os
from scipy.interpolate import interp1d
import sys
import uuid

#logging.getLogger().setLevel(logging.INFO)
#logging.basicConfig(format='Adding Nebular Emission : %(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

"""
----------------------------------------------------------------------------------------------------------------
Based on cloudyfsps written by Nell Byler.
(Source https://github.com/nell-byler/cloudyfsps/tree/master/cloudyfsps retrieved in October 2019)
and FSPS written by Charlie Conroy
(Source https://github.com/cconroy20/fsps/blob/master/src/add_nebular.f90 retrieved in October 2019)
----------------------------------------------------------------------------------------------------------------
"""


def write_input_sed(wav, spec):
    """
    ---------------------------------------------------------------------
    Writes a input SED file for CLOUDY to use
    Assumes you have $CLOUDY_EXE and $CLOUDY_DATA_PATH set as sys vars.
    ---------------------------------------------------------------------
    """
    ascii_file = str(uuid.uuid4().hex) + ".ascii"
    while compiled_exists(ascii_file):
        ascii_file = str(uuid.uuid4().hex) + ".ascii"

    logging.info("Executing write ascii sequence...")
    logging.info("Writing.....")
    WriteASCII(ascii_file, wav, spec, nx=len(wav), nmod=1, par1_val=1.e6)

    logging.info("Compiling {} with Cloudy".format(ascii_file))
    compile_ascii(ascii_file)

    logging.info("Checking to see if compilation was successful...")
    if check_compiled(ascii_file):
        logging.info("Your model {} is ready to run.".format(ascii_file))
    else:
        sys.exit()

    return ascii_file


def write_cloudy_input(**kwargs):
    """
    ----------------------------------------
    Writes standard cloudy input file
    ----------------------------------------
    """
    pars = {"dir_": "/home/prerak/codes/modeling-DIG/cloudy/",
            "model_name": "temp_sed",
            "age": 1.0e6,      # age in years
            "logZ": -1.5,      # logZ/Zsol (-2.0 to 0.2)
            "gas_logZ": None,
            "logQ": 47.576,
            "logU": -3.0,
            "dens": 100.0,     # number density of hydrogen
            "r_inner": 1.e19,  # inner radius of cloud
            "r_in_pc": False,
            "use_Q": True,
            "dust": True,
            "efrac": -1.0,
            "geometry": "sphere",
            "abundance": "dopita"
            }
    for key, value in list(kwargs.items()):
        pars[key] = value

    file_name = os.path.join(pars["dir_"]+"/temp_files", pars["model_name"] + ".in")
    f = open(file_name, "w")

    def this_print(s, eol=True):
        if s is None:
            print('"None" parameter not printed')
        else:
            to_print = s.strip()
            if eol:
                to_print += "\n"
            f.write(to_print)

    # -----------------------------------------------------
    pars["gas_logZ"] = pars["logZ"]
    this_print('////////////////////////////////////')
    this_print('title {0}'.format(pars["model_name"].split('/')[-1]))
    this_print('////////////////////////////////////')
    this_print('set punch prefix "{0}"'.format(pars["model_name"]))
    this_print('print line precision 6')
    cloudy_mod = "{}.mod".format(pars["model_name"])
    this_print('table star "{0}" {1}={2:.2e}'.format(cloudy_mod, "age", pars['age']))
    if pars['use_Q']:
        this_print('Q(H) = {0:.3f} log'.format(pars['logQ']))
    else:
        this_print('ionization parameter = {0:.3f} log'.format(pars['logU']))

    if pars['dust']:
        this_print('grains {0:.2f} log'.format(pars['gas_logZ']))

    if pars['r_in_pc']:
        pc_to_cm = 3.08568e18
        r_out = np.log10(pars['r_inner'] * pc_to_cm)
    else:
        r_out = np.log10(pars['r_inner'])

    linefile = "cloudyLines.dat"
    abunds = getNebAbunds(pars["abundance"],
                          pars["logZ"],
                          dust=pars["dust"],
                          re_z=False)
    this_print(abunds.solarstr)
    for line in abunds.elem_strs:
        this_print(line)

    this_print('radius {0:.3f} log'.format(r_out))
    this_print('hden {0:.3f} log'.format(np.log10(pars['dens'])))
    this_print('{}'.format(pars['geometry']))
    this_print('cosmic ray background')
    this_print('iterate to convergence max=5')
    this_print('stop temperature 100.0')
    this_print('stop efrac {0:.2f}'.format(pars['efrac']))
    this_print('save last linelist ".lin" "{}" absolute column'.format(linefile))
    this_print('save last outward continuum ".outwcont" units Angstrom no title')
    this_print('save last incident continuum ".inicont" units Angstrom no title')
    print("Input written in {0}".format(file_name))
    f.close()


def get_output(model_name, dir_, qq):
    lsun = 3.839e33
    c = constants.c.cgs.value * 1.e8

    outcontfl = dir_ + "/temp_files/"+ model_name + ".outwcont"
    outlinefl = dir_ + "/temp_files/" + model_name + ".lin"
    refline_file = dir_ + "/data/refLines.dat"
    fsps_lam_file = dir_ + "/data/FSPSlam.dat"

    dat = np.genfromtxt(outlinefl, delimiter="\t", skip_header=2, dtype="S20,f8")
    datflu = np.array([d[1] for d in dat])
    wdat = np.genfromtxt(refline_file, delimiter=',')
    wl = np.array([dat[0] for dat in wdat])
    sinds = np.argsort(wl)
    wl = wl[sinds]

    datflu = datflu[sinds] / lsun / qq
    datflu = [float("{0:1.4e}".format(dat)) for dat in datflu]
    cont_data = np.genfromtxt(outcontfl, skip_header=1)
    fsps_lam = np.genfromtxt(fsps_lam_file)
    nu = c / fsps_lam
    ang_0, diffuse_0 = cont_data[:, 0], cont_data[:, 2]
    ang_v, diffuse_in = air_to_vac(ang_0[::-1]), diffuse_0[::-1]
    diffuse_y = interp1d(ang_v, diffuse_in, fill_value=0.0, bounds_error=False)(fsps_lam)
    diffuse_out = diffuse_y / nu / lsun / qq
    return wl, datflu, fsps_lam, diffuse_out


def clean_files(dir_, model_name, error=False):
    logging.info("Cleaning up temporary files")
    if not error:
        os.remove(os.path.join(dir_ + "/temp_files", model_name + ".out"))
        os.remove(os.path.join(dir_ + "/temp_files", model_name + ".in"))
        os.remove(os.path.join(dir_ + "/temp_files", model_name + ".inicont"))
        os.remove(os.path.join(dir_ + "/temp_files", model_name + ".lin"))
        os.remove(os.path.join(dir_ + "/temp_files", model_name + ".outwcont"))
        os.remove(os.path.join(os.environ['CLOUDY_DATA_PATH'], model_name + ".ascii"))

    os.remove(os.path.join(os.environ['CLOUDY_DATA_PATH'], model_name + ".mod"))
    os.remove(os.path.join(os.environ['CLOUDY_DATA_PATH'], model_name + ".in"))
    os.remove(os.path.join(os.environ['CLOUDY_DATA_PATH'], model_name + ".out"))


def get_nebular(spec_lambda, sspi, nh, logq, radius, logu, logz, logq_1, abund="dopita", useq=True, clean_up=True):
    nspec = len(spec_lambda)
    frac_obrun = 0.0
    clight = constants.c.cgs.value*1.e8
    nebular_smooth_init = 0

    logging.info("Writing the input SED file")
    filename = write_input_sed(spec_lambda, sspi)
    model_name = filename.split(".")[0]

    logging.info("Writing CLOUDY input file")
    dir_= os.getcwd() + "/powderday/nebular_emission"
    dir_base = os.getcwd()

    write_cloudy_input(dir_=dir_,
                       model_name=model_name,
                       r_inner=radius,
                       dens=nh,
                       use_Q=useq,
                       logQ=logq,
                       logU=logu,
                       logZ=logz,
                       abundance=abund,
                       dust=False)

    logging.info("Input SED file written")
    logging.info("Running CLOUDY")
    os.chdir(dir_ + "/temp_files/")
    os.system("$CLOUDY_EXE " + model_name + ".in")
    os.chdir(dir_base)
    f_out = open(dir_ + "/temp_files/"+model_name+".out", 'r')
    content = f_out.readlines()
    check = np.all(['OK' in content[-1]])
    if check:
        logging.info("CLOUDY run finished sucessfully")

    else:
        logging.info("CLOUDY run was unsucessful")
        if clean_up:
            clean_files(dir_, model_name, error=True)
        raise ValueError('CLOUDY run was unsucessful')

    logging.info("Getting output spectrum")

    nebem_line_pos, nebem_line, readlambneb, readcontneb = get_output(model_name, dir_, 10**logq)
    nebem_cont = interp1d(readlambneb, np.log10(readcontneb + 10 ** (-95.0)),
                          fill_value=0, bounds_error=False)(spec_lambda)

    logging.info("Adding nebular emission to input spectrum")
    nemline = sum(1 for line in open(dir_ + '/temp_files/' + 'cloudyLines.dat'))

    neb_res_min = np.zeros(nemline)
    for i in range(nemline):
        max_id = max(np.max(np.where(spec_lambda <= nebem_line_pos[i])), 1)
        j = min(max_id, nspec - 1)
        neb_res_min[i] = spec_lambda[j + 1] - spec_lambda[j]

    gaussnebarr = []
    for i in range(nemline):
        dlam = nebem_line_pos[i] * nebular_smooth_init / clight * 1e13
        dlam = max(dlam, neb_res_min[i] * 2)
        gaussnebarr.append(1 / np.sqrt(2 * np.pi) / dlam * np.exp(-(spec_lambda - nebem_line_pos[i]) ** 2 / 2 /
                                                                  dlam ** 2) / clight * nebem_line_pos[i] ** 2)
    gaussnebarr = np.array(gaussnebarr)

    sspo = np.copy(sspi)
    whlylim = np.max(np.where(spec_lambda <= 912.))

    sspo[: whlylim] = sspi[: whlylim] * max(min(frac_obrun, 1.0), 0.0)
    sspo = sspo + 10 ** nebem_cont * (10**logq_1)
    for i in range(nemline):
        sspo = sspo + nebem_line[i] * gaussnebarr[i] * (10**logq_1)

    if clean_up:
        clean_files(dir_, model_name,  error=False)

    return sspo

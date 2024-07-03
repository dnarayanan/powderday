from powderday.nebular_emission.abund import getNebAbunds
from powderday.nebular_emission.ASCIItools import *
from powderday.nebular_emission.cloudy_tools import air_to_vac, calc_LogQ, convert_metals
import powderday.config as cfg
from astropy import constants
import numpy as np
import os
from scipy.interpolate import interp1d
import sys
import uuid

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
      
    # Writing input SED file
    WriteASCII(ascii_file, wav, spec, nx=len(wav), nmod=1, par1_val=1.e6)

    # Compiling the input SED file with Cloudy
    compile_ascii(ascii_file)

    # Checking to see if compilation was successful
    if not check_compiled(ascii_file):
        ValueError('CLOUDY run was unsucessful')

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
            "dust": False,
            "efrac": 0.0,
            "abundance": "dopita",
            "id_val": 1,
            "metals": [-1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1.]
            }
    for key, value in list(kwargs.items()):
        pars[key] = value

    file_name = os.path.join(pars["dir_"]+"/temp_files", pars["model_name"] + ".in")
    f = open(file_name, "w+")

    def this_print(s, eol=True):
        if s is None:
            print('"None" parameter not printed')
        else:
            to_print = s.strip()
            if eol:
                to_print += "\n"
            f.write(to_print)

    # -----------------------------------------------------
    _id = int(pars["id_val"])

    pars["gas_logZ"] = pars["logZ"]
    this_print('////////////////////////////////////')
    this_print('title {0}'.format(pars["model_name"].split('/')[-1]))
    this_print('////////////////////////////////////')
    this_print('set punch prefix "{0}"'.format(pars["model_name"]))
    this_print('print line precision 6')
    
    cloudy_mod = "{}.mod".format(pars["model_name"])
    this_print('table star "{0}" {1}={2:.2e}'.format(cloudy_mod, "age", pars['age']))
    
    # If DIG (id = 3) then use the ionization parameter else use the rate of ionizing photons
    if _id == 3:
        this_print('ionization parameter = {0:.3f} log'.format(pars['logU']))
    else:
        this_print('Q(H) = {0:.3f} log'.format(pars['logQ']))
    
    if pars['r_in_pc']:
        pc_to_cm = 3.08568e18
        r_out = np.log10(pars['r_inner'] * pc_to_cm)
    else:
        r_out = np.log10(pars['r_inner'])

    
    linefile = "cloudyLines.dat"

    # If there was an error in getting metallicities from the simulation
    # then the code reverts back to using "dopita" abundances in place of "direct"
    if (any(q <= -1.0 for q in pars["metals"][1:]) and pars["abundance"] == "direct"):
        pars["abundance"] = "dopita"
        print ("Warning: Unable to get metallicities from the simulation. This can be because you are binning the stars.\n \
                Make sure that FORCE BINNED is set to False and set the max_age_direct appropriately to use this feature.\n \
                Reverting to using abundances from \"dopita et al. 2001\"")

    if pars["abundance"] == "direct":
        if pars['dust']:
            this_print('grains orion {0:.2f} log no qheat'.format(pars['gas_logZ']))

        abund_el = ['He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe']
        abund_metal = convert_metals(pars["metals"][1:])
        abund_str = "abundances "
    
        if cfg.par.FORCE_N_O_Pilyugin[_id]: # Adding N/O pilyugin relation
            if 12 + abund_metal[3] < 8.14:
                abund_metal[2] = -1.493 + abund_metal[3]
            else:
                abund_metal[2] =  1.489*(12 + abund_metal[3]) - 13.613 + abund_metal[3]
       
        if cfg.par.FORCE_N_O_ratio[_id]:
            abund_metal[2] = cfg.par.N_O_ratio[_id] + abund_metal[3]

        # Enhancing abundances for post-AGB stars
        if _id == 1:
            abund_metal[1] += cfg.par.PAGB_C_enhancement
            abund_metal[2] += cfg.par.PAGB_N_enhancement

        for e in range(len(abund_el)):
            el_str = str(abund_el[e]) + " " + str(abund_metal[e]) + " "
            abund_str = abund_str + el_str

        abund_extra = "init file=\"hii.ini\""
        this_print(abund_extra)
        this_print(abund_str)

    else:
        abunds = getNebAbunds(pars["abundance"],
                            pars["logZ"],
                            dust=pars["dust"],
                            re_z=False)

        abund_C = float(abunds.elem_strs[1].split(" ")[3])
        abund_N = float(abunds.elem_strs[2].split(" ")[3])
        abund_O = float(abunds.elem_strs[3].split(" ")[3])
        
        if cfg.par.FORCE_N_O_Pilyugin[_id]: # Forcing the Nitrogen abundances to follow the N/O relation from Pilyugin et al. 2012
            if 12 + abund_N < 8.14:
                abund_N = -1.493 + abund_O
            else:
                abund_N =  1.489*(12 + abund_O) - 13.613 + abund_O

        if cfg.par.FORCE_N_O_ratio[_id]:
            abund_N = cfg.par.N_O_ratio[_id] + abund_O

        # Enhancing abundances for post-AGB stars
        if _id == 1:
            abund_C += cfg.par.PAGB_C_enhancement
            abund_N += cfg.par.PAGB_N_enhancement

        abunds.elem_strs[1] = "element abundance carbon "+str(abund_C)+" log"
        abunds.elem_strs[2] = "element abundance nitrogen "+str(abund_N)+" log"
        abunds.elem_strs[3] = "element abundance oxygen "+str(abund_O)+" log"

        this_print(abunds.solarstr)
        for line in abunds.elem_strs:
            this_print(line)
    
    cf = 1 - pars['efrac']

    # For DIG, since we are making use of the ionization parameter we do need
    # to set the inner radius, the geometry is plane parallel
    # and the covering factor is set to 1. 
    if _id != 3:
        this_print('radius {0:.3f} log'.format(r_out))
        this_print('sphere')
        this_print('Covering factor {0:.3f}'.format(cf))
    
    this_print('hden {0:.3f} log'.format(np.log10(pars['dens'])))
    this_print('cosmic ray background')
    this_print('iterate to convergence max=5')
    this_print('stop temperature 100.0')
    this_print('stop efrac -1.0')
    this_print('save last linelist ".lin" "{}" absolute column'.format(linefile))
    this_print('save last outward continuum ".outwcont" units Angstrom no title')
    this_print('save last incident continuum ".inicont" units Angstrom no title')
    f.close()


def get_output(model_name, dir_, qq, fsps_lam, cell_width, id_val):
    lsun = 3.839e33
    c = constants.c.cgs.value * 1.e8

    outcontfl = dir_ + "/temp_files/"+ model_name + ".outwcont"
    outlinefl = dir_ + "/temp_files/" + model_name + ".lin"
    refline_file = dir_ + "/data/refLines.dat"

    dat = np.genfromtxt(outlinefl, delimiter="\t", skip_header=2, dtype="S20,f8")
    datflu = np.array([d[1] for d in dat])
    wdat = np.genfromtxt(refline_file, delimiter=',')
    wl = np.array([dat[0] for dat in wdat])
    sinds = np.argsort(wl)
    wl = wl[sinds]
     
    # For DIG we are making use of the ionization parameter (intensity case, see CLOUDY docs (Hazy 1))
    # We need to multiply the output by the area of the cloud -- 6 sides of the cell in this case
    if id_val == 3:
        datflu = (datflu[sinds]*6*cell_width**2) / lsun 

    else:
        datflu = datflu[sinds] / lsun / qq
    
    datflu = [float("{0:1.4e}".format(dat)) for dat in datflu]

    cont_data = np.genfromtxt(outcontfl, skip_header=1)
    nu = c / fsps_lam
    ang_0, diffuse_0 = cont_data[:, 0], cont_data[:, 2]
    ang_v, diffuse_in = air_to_vac(ang_0[::-1]), diffuse_0[::-1]
    diffuse_y = interp1d(ang_v, diffuse_in, fill_value=0.0, bounds_error=False)(fsps_lam)
    
    if id_val == 3:
        diffuse_out = (diffuse_y*6*cell_width**2) / nu / lsun 
    else:
        diffuse_out = diffuse_y / nu / lsun / qq
    
    return wl, datflu, fsps_lam, diffuse_out


def clean_files(dir_, model_name, id_val, error=False):
    # Cleaning up temporary files
    
    os.remove(os.path.join(dir_ + "/temp_files", model_name + ".out"))
    os.remove(os.path.join(dir_ + "/temp_files", model_name + ".in"))
    os.remove(os.path.join(dir_ + "/temp_files", model_name + ".inicont"))
    os.remove(os.path.join(dir_ + "/temp_files", model_name + ".lin"))
    os.remove(os.path.join(dir_ + "/temp_files", model_name + ".outwcont"))
    if id_val != 3:
        os.remove(os.path.join(os.environ['CLOUDY_DATA_PATH'], model_name + ".ascii"))
        
    if id_val != 3:
        os.remove(os.path.join(os.environ['CLOUDY_DATA_PATH'], model_name + ".mod"))
        os.remove(os.path.join(os.environ['CLOUDY_DATA_PATH'], model_name + ".in"))
        os.remove(os.path.join(os.environ['CLOUDY_DATA_PATH'], model_name + ".out"))


def get_nebular(spec_lambda, sspi, nh, Metals, logq = 0.0, radius = 1.e19, logu = 0.0, logz = 0.0, logq_1=0.0, 
        Cell_width=1, Dust=False, abund="dopita", clean_up=True, index=1, efrac=0.0):
    
    # Since we don't need the logQs for DIG as we are using 
    # the ionization parmeter instead. Set them to 0 just as a place holder 
    if index == 3:
        logq = 0.0
        logq_1 = 0.0
    
    nspec = len(spec_lambda)
    frac_obrun = efrac
    clight = constants.c.cgs.value*1.e8
    nebular_smooth_init = 0

    # Writing the input SED file
    filename = write_input_sed(spec_lambda, sspi)
    model_name = filename.split(".")[0]

    # Writing CLOUDY input file
    dir_= cfg.par.pd_source_dir + "powderday/nebular_emission"
    dir_base = os.getcwd()
    write_cloudy_input(dir_=dir_,
                       model_name = model_name,
                       r_inner = radius,
                       dens = nh,
                       logQ = logq,
                       logU = logu,
                       logZ = logz,
                       abundance = abund,
                       dust = Dust,
                       id_val = index,
                       efrac = frac_obrun,
                       metals = Metals)

    # Running CLOUDY
    os.chdir(dir_ + "/temp_files/")
    os.system("$CLOUDY_EXE " + model_name + ".in")
    os.chdir(dir_base)
    f_out = open(dir_ + "/temp_files/"+model_name+".out", 'r')
    content = f_out.readlines()
    check = np.all(['OK' in content[-1]])
    
    # Checking if CLOUDY run finished sucessfully
    if not check:
        print ("WARNING: The CLOUDY run was unsucessful. Please see the cloudy output file "+ model_name+".out"+
               " located in "+ dir_ + "/temp_files/"+" to figure out why this happened.")
        if clean_up:
            clean_files(dir_, model_name, index, error=True)
        raise ValueError('CLOUDY run was unsucessful')

    
    # Getting output spectrum
    nebem_line_pos, nebem_line, readlambneb, readcontneb = get_output(model_name, dir_, 10**logq, spec_lambda, Cell_width, index)
    nebem_cont = interp1d(readlambneb, np.log10(readcontneb + 10 ** (-95.0)),
                          fill_value=-95.0, bounds_error=False)(spec_lambda)

    # Adding nebular emission to input spectrum
    nemline = len(np.where(nebem_line_pos < spec_lambda[-1])[0])
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
        clean_files(dir_, model_name, index, error=False)

    return sspo, nebem_line_pos, np.array(nebem_line)*(10**logq_1)

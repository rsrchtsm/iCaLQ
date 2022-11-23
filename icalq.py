# Calculator Imports
import sys
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
import sympy as sym
from sympy.utilities.lambdify import lambdify
from sympy.utilities.iterables import flatten
from sympy import sympify
import math
from copy import deepcopy
import random
from functools import cmp_to_key
import scipy.optimize as optimize
# import readline

# Use this function as the comparator function while sorting lambdas
def compare_lambda(item1, item2):
    a1 = list(item1[0])
    a2 = list(item2[0])
    if a1[3] != a2[3]:
        return ord(a1[3]) - ord(a2[3])
    else:
        if a1[4] != a2[4]:
            if a1[4] == 'L':
                return -1
            return 1
        else:
            return ord(a1[2]) - ord(a2[2])
    return 0

# parsing string input to their appropriate format. Outputs to be used by home function.
def parse(mass_f, lambdas_f, ignore_f, margin_f, lam_values_f):
    mass = float(mass_f)
    margin = float(margin_f)
    original_lambdastring = lambdas_f.replace(',' , ' ').strip().split()
    ignorePairSingle = False
    original_lam_vals = []
    temp_lam_vals = []
    lambdastring = []
    if (ignore_f.strip().lower() == "yes" or ignore_f.strip().lower() == 'y'):
        ignorePairSingle = True
    for val in lam_values_f:
        val = val.replace(',', ' ')
        try:
            if len(val.strip().split()) != len(original_lambdastring):
                raise ValueError()
            for x in val.strip().split():
                float(x)
        except ValueError:
            print(original_lambdastring)
            sys.exit(f"[Query Error]: Query values for lambdas are either not {len(original_lambdastring)} (number of lambdas) in count or not convertible to float.")
        original_lam_vals.append(val.strip().split())
    for lam_val in original_lam_vals:
        combined_lambda = zip(original_lambdastring, lam_val)
        combined_lambda = sorted(combined_lambda, key=cmp_to_key(compare_lambda))
        combined_lambda = list(zip(*combined_lambda))
        lambdastring = list(combined_lambda[0])
        temp_lam_vals.append(list(combined_lambda[1]))
    lam_vals = temp_lam_vals
    return mass, lambdastring, original_lambdastring, ignorePairSingle, lam_vals, original_lam_vals, margin


# parsing queries while in interactive mode
def parse_lam(original_lambdastring, lam_val_f):
    temp_lam_vals = []
    lam_val_f = lam_val_f.replace(',', ' ')
    original_lam_vals = [lam_val_f.strip().split()]
    combined_lambda = zip(original_lambdastring, original_lam_vals[0])
    combined_lambda = sorted(combined_lambda, key=cmp_to_key(compare_lambda))
    combined_lambda = list(zip(*combined_lambda))
    temp_lam_vals.append(list(combined_lambda[1]))
    return temp_lam_vals, original_lam_vals

interpolation_type='slinear'
def interpolate_cs_func(df, ls):
    return lambda mass: [interp1d(data_mass_list, df[coupling][:5], kind=interpolation_type)([mass])[0] for coupling in ls]

# Interpolating cross-section of t-channel's cross terms
def interpolate_cs_ct_func(df):
    return lambda mass: [interp1d(data_mass_list, df[ij], kind=interpolation_type)([mass])[0] for ij in range(len(df)) ]



# Declaring variables which need not be reloaded every run
chi_sq_limits_1 = [1.00, 2.295748928898636, 3.5267403802617303, 4.719474460025883, 5.887595445915204, 7.038400923736641, 8.176236497856527, 9.30391276903717, 10.423363154355838, 11.535981713319316, 12.64281133339149, 13.744655587189282, 14.842148802786893, 15.935801892195538, 17.026033423371082, 18.11319133873574, 19.197568537049687]
chi_sq_limits_2 = [4.00, 6.180074306244173, 8.024881760266252, 9.715627154871333, 11.313855908361862, 12.848834791793395, 14.337110231671799, 15.789092974617745, 17.21182898078949, 18.610346565823498, 19.988381717650192, 21.348799569984315, 22.693854280452445, 24.025357063756637, 25.344789151124267, 26.653380234523553, 27.952164463248984]
chi_sq_limits = chi_sq_limits_2
data_mass_list = [1000, 1500, 2000, 2500, 3000]
luminosity_tau = 139 * 1000
luminosity_e_mu = 140 * 1000
standard_HHbT = pd.read_csv('./data/HEPdata/HHbT.csv',header=[0])
standard_HHbV = pd.read_csv('./data/HEPdata/HHbV.csv',header=[0])
standard_LHbT = pd.read_csv('./data/HEPdata/LHbT.csv',header=[0])
standard_LHbV = pd.read_csv('./data/HEPdata/LHbV.csv',header=[0])
standard_ee = pd.read_csv('./data/HEPdata/dielectron.csv',header=[0])
standard_mumu = pd.read_csv('./data/HEPdata/dimuon.csv',header=[0])
ND = [standard_HHbT['ND'].to_numpy(), standard_HHbV['ND'].to_numpy(),
      standard_LHbT['ND'].to_numpy(), standard_LHbV['ND'].to_numpy(),
      standard_ee['ND'].to_numpy(),standard_mumu['ND'].to_numpy()]
NSM = [standard_HHbT['Standard Model'].to_numpy(), standard_HHbV['Standard Model'].to_numpy(),
       standard_LHbT['Standard Model'].to_numpy(), standard_LHbV['Standard Model'].to_numpy(),
      standard_ee['Standard Model'].to_numpy(),standard_mumu['Standard Model'].to_numpy()]


# Path variables for cross section: 
cs_sc_path="./data/cross_section/"
df_pair = pd.read_csv(cs_sc_path + "pair.csv")
df_single = pd.read_csv(cs_sc_path + "single.csv")
df_interference = pd.read_csv(cs_sc_path + "interference.csv")
df_tchannel = pd.read_csv(cs_sc_path + "tchannel.csv")
df_pureqcd = pd.read_csv(cs_sc_path + "pureqcd.csv")
cross_terms_tchannel = "./data/cross_section/tchannel_doublecoupling.csv"
double_coupling_data_tchannel = pd.read_csv(cross_terms_tchannel, header=[0])


# Get cross sections from data files
def get_cs(mass, lambdastring, num_lam):
    cs_q = interpolate_cs_func(df_pureqcd, lambdastring)
    cs_p = interpolate_cs_func(df_pair, lambdastring)
    cs_s = interpolate_cs_func(df_single, lambdastring)
    cs_i = interpolate_cs_func(df_interference, lambdastring)
    cs_t = interpolate_cs_func(df_tchannel, lambdastring)
    cs_l = [cs_q(mass), cs_p(mass), cs_s(mass), cs_i(mass), cs_t(mass)]
    #
    ee_cs = []
    mumu_cs = []
    tautau_cs = []
    for process in cs_l:
        ee_temp = []
        mumu_temp = []
        tautau_temp = []
        for lamda,cs in zip(lambdastring,process):
            if lamda[3] == '1':
                ee_temp.append(cs)
            elif lamda[3] == '2':
                mumu_temp.append(cs)
            elif lamda[3] == '3':
                tautau_temp.append(cs)
        ee_cs.append(ee_temp)
        mumu_cs.append(mumu_temp)
        tautau_cs.append(tautau_temp)
    # cross terms: 
    ee_t_ct = [double_coupling_data_tchannel[lambdastring[i] + '_' + lambdastring[j]] for i in range(num_lam) for j in range(i+1, num_lam) if lambdastring[i][3] == lambdastring[j][3] and lambdastring[i][3] == '1']
    mumu_t_ct = [double_coupling_data_tchannel[lambdastring[i] + '_' + lambdastring[j]] for i in range(num_lam) for j in range(i+1, num_lam) if lambdastring[i][3] == lambdastring[j][3] and lambdastring[i][3] == '2']
    tautau_t_ct = [double_coupling_data_tchannel[lambdastring[i] + '_' + lambdastring[j]] for i in range(num_lam) for j in range(i+1, num_lam) if lambdastring[i][3] == lambdastring[j][3] and lambdastring[i][3] == '3']
    cs_ee_t_ct_func = interpolate_cs_ct_func(ee_t_ct)
    cs_ee_t_ct_temp = cs_ee_t_ct_func(mass)
    cs_mumu_t_ct_func = interpolate_cs_ct_func(mumu_t_ct)
    cs_mumu_t_ct_temp = cs_mumu_t_ct_func(mass)
    cs_tautau_t_ct_func = interpolate_cs_ct_func(tautau_t_ct)
    cs_tautau_t_ct_temp = cs_tautau_t_ct_func(mass)
    #
    ee_cntr = 0
    cs_ee_t_ct = cs_ee_t_ct_temp[:]
    mumu_cntr = 0
    cs_mumu_t_ct = cs_mumu_t_ct_temp[:]
    tautau_cntr = 0
    cs_tautau_t_ct = cs_tautau_t_ct_temp[:]
    for i in range(num_lam):
        for j in range(i+1, num_lam):
            if lambdastring[i][3] == lambdastring[j][3]:
                if lambdastring[i][3] == '1':
                    cs_ee_t_ct[ee_cntr] = cs_ee_t_ct_temp[ee_cntr] - cs_l[4][i] - cs_l[4][j]
                    ee_cntr += 1
                elif lambdastring[i][3] == '2':
                    cs_mumu_t_ct[mumu_cntr] = cs_mumu_t_ct_temp[mumu_cntr] - cs_l[4][i] - cs_l[4][j]
                    mumu_cntr += 1
                elif lambdastring[i][3] == '3':
                    cs_tautau_t_ct[tautau_cntr] = cs_tautau_t_ct_temp[tautau_cntr] - cs_l[4][i] - cs_l[4][j]
                    tautau_cntr += 1
    return [ee_cs, mumu_cs, tautau_cs, cs_ee_t_ct, cs_mumu_t_ct, cs_tautau_t_ct, cs_ee_t_ct_temp, cs_mumu_t_ct_temp, cs_tautau_t_ct_temp]

# Find the closest mass to the input mass
def get_closest_mass(mass):
    closest_mass = 0
    closest_diff = 10000
    for val in data_mass_list:
        if abs(mass-val) < closest_diff:
            closest_diff = abs(mass-val)
            closest_mass = val
    return closest_mass


efficiency_prefix = "./data/efficiency/"
t_ct_prefix = "./data/efficiency/t/"
tagnames = ["/HHbT.csv", "/HHbV.csv", "/LHbT.csv", "/LHbV.csv"]

# Load efficiencies from the data files
def get_efficiencies(closest_mass, lambdastring, num_lam, cs_list):
    [ee_cs, mumu_cs, tautau_cs, cs_ee_t_ct, cs_mumu_t_ct, cs_tautau_t_ct, cs_ee_t_ct_temp, cs_mumu_t_ct_temp, cs_tautau_t_ct_temp] = cs_list

    path_interference_ee = [efficiency_prefix + "i/" + str(coupling[2:]) for coupling in lambdastring if coupling[3]=='1']
    path_pair_ee = [efficiency_prefix + "p/" + str(coupling[2:]) for coupling in lambdastring if coupling[3]=='1']
    path_single_ee = [efficiency_prefix + "s/" + str(coupling[2:]) for coupling in lambdastring if coupling[3]=='1']
    path_tchannel_ee = [efficiency_prefix + "t/" + str(coupling[2:]) for coupling in lambdastring if coupling[3]=='1']
    path_pureqcd_ee = [efficiency_prefix + "q/" + str(coupling[2:]) for coupling in lambdastring if coupling[3]=='1']

    path_interference_mumu = [efficiency_prefix + "i/" + str(coupling[2:]) for coupling in lambdastring if coupling[3]=='2']
    path_pair_mumu = [efficiency_prefix + "p/" + str(coupling[2:]) for coupling in lambdastring if coupling[3]=='2']
    path_single_mumu = [efficiency_prefix + "s/" + str(coupling[2:]) for coupling in lambdastring if coupling[3]=='2']
    path_tchannel_mumu = [efficiency_prefix + "t/" + str(coupling[2:]) for coupling in lambdastring if coupling[3]=='2']
    path_pureqcd_mumu = [efficiency_prefix + "q/" + str(coupling[2:]) for coupling in lambdastring if coupling[3]=='2']

    path_interference_tautau = [efficiency_prefix + "i/" + str(coupling[2:]) + "/" + str(closest_mass) for coupling in lambdastring if coupling[3]=='3']
    path_pair_tautau = [efficiency_prefix + "p/" + str(coupling[2:]) + "/" + str(closest_mass) for coupling in lambdastring if coupling[3]=='3']
    path_single_tautau = [efficiency_prefix + "s/" + str(coupling[2:]) + "/" + str(closest_mass) for coupling in lambdastring if coupling[3]=='3']
    path_tchannel_tautau = [efficiency_prefix + "t/" + str(coupling[2:]) + "/" + str(closest_mass) for coupling in lambdastring if coupling[3]=='3']
    path_pureqcd_tautau = [efficiency_prefix + "q/" + str(coupling[2:]) + "/" + str(closest_mass) for coupling in lambdastring if coupling[3]=='3']
    
    ee_path_t_ct = []
    mumu_path_t_ct = []
    tautau_path_t_ct = []
    for i in range(num_lam):
        for j in range(i+1, num_lam):
            if lambdastring[i][3] == lambdastring[j][3]:
                if lambdastring[i][3] == '1':
                    ee_path_t_ct.append(t_ct_prefix + str(lambdastring[i][2:]) + "_" + str(lambdastring[j][2:]) )
                elif lambdastring[i][3] == '2':
                    mumu_path_t_ct.append(t_ct_prefix + str(lambdastring[i][2:]) + "_" + str(lambdastring[j][2:]) )
                elif lambdastring[i][3] == '3':
                    tautau_path_t_ct.append(t_ct_prefix + str(lambdastring[i][2:]) + "_" + str(lambdastring[j][2:]) + "/" + str(closest_mass))

    ee_eff_l = [[[pd.read_csv(path_pureqcd_ee[i] + "/" + str(int(closest_mass)) + ".csv",header=[0]).to_numpy()[:,2]] for i in range(len(path_pureqcd_ee))],
             [[pd.read_csv(path_pair_ee[i] + "/" + str(int(closest_mass)) + ".csv",header=[0]).to_numpy()[:,2]] for i in range(len(path_pureqcd_ee))],
             [[pd.read_csv(path_single_ee[i] + "/" + str(int(closest_mass)) + ".csv",header=[0]).to_numpy()[:,2]] for i in range(len(path_pureqcd_ee))],
             [[pd.read_csv(path_interference_ee[i] + "/" + str(int(closest_mass)) + ".csv",header=[0]).to_numpy()[:,2]] for i in range(len(path_pureqcd_ee))],
             [[pd.read_csv(path_tchannel_ee[i] + "/" + str(int(closest_mass)) + ".csv",header=[0]).to_numpy()[:,2]] for i in range(len(path_pureqcd_ee))]]

    mumu_eff_l = [[[pd.read_csv(path_pureqcd_mumu[i] + "/" + str(int(closest_mass)) + ".csv",header=[0]).to_numpy()[:,2]] for i in range(len(path_pureqcd_mumu))],
             [[pd.read_csv(path_pair_mumu[i] + "/" + str(int(closest_mass)) + ".csv",header=[0]).to_numpy()[:,2]] for i in range(len(path_pureqcd_mumu))],
             [[pd.read_csv(path_single_mumu[i] + "/" + str(int(closest_mass)) + ".csv",header=[0]).to_numpy()[:,2]] for i in range(len(path_pureqcd_mumu))],
             [[pd.read_csv(path_interference_mumu[i] + "/" + str(int(closest_mass)) + ".csv",header=[0]).to_numpy()[:,2]] for i in range(len(path_pureqcd_mumu))],
             [[pd.read_csv(path_tchannel_mumu[i] + "/" + str(int(closest_mass)) + ".csv",header=[0]).to_numpy()[:,2]] for i in range(len(path_pureqcd_mumu))]]


    tautau_eff_l = [[[pd.read_csv(path_pureqcd_tautau[i] + j,header=[0]).to_numpy()[:,2] for j in tagnames] for i in range(len(path_pureqcd_tautau))],
             [[pd.read_csv(path_pair_tautau[i] + j,header=[0]).to_numpy()[:,2] for j in tagnames] for i in range(len(path_pureqcd_tautau))],
             [[pd.read_csv(path_single_tautau[i] + j,header=[0]).to_numpy()[:,2] for j in tagnames] for i in range(len(path_pureqcd_tautau))],
             [[pd.read_csv(path_interference_tautau[i] + j,header=[0]).to_numpy()[:,2] for j in tagnames] for i in range(len(path_pureqcd_tautau))],
             [[pd.read_csv(path_tchannel_tautau[i] + j,header=[0]).to_numpy()[:,2] for j in tagnames] for i in range(len(path_pureqcd_tautau))]]

    ee_eff_t_ct_temp = [[pd.read_csv(ee_path_t_ct[i] + "/" + str(int(closest_mass)) + ".csv",header=[0]).to_numpy()[:,2]] for i in range(len(ee_path_t_ct))]
    mumu_eff_t_ct_temp = [[pd.read_csv(mumu_path_t_ct[i] + "/" + str(int(closest_mass)) + ".csv",header=[0]).to_numpy()[:,2]] for i in range(len(mumu_path_t_ct))]
    tautau_eff_t_ct_temp = [[pd.read_csv(tautau_path_t_ct[j] + i,header=[0]).to_numpy()[:,2] for i in tagnames]  for j in range(len(tautau_path_t_ct))]

    ee_eff_t_ct = deepcopy(ee_eff_t_ct_temp)
    mumu_eff_t_ct = deepcopy(mumu_eff_t_ct_temp)
    tautau_eff_t_ct = deepcopy(tautau_eff_t_ct_temp)

    ee_cntr = 0
    mumu_cntr = 0
    tautau_cntr = 0
    for i in range(len(path_interference_ee)):
        for j in range(i+1, len(path_interference_ee)):
            ee_eff_t_ct[ee_cntr][0] = (ee_eff_t_ct_temp[ee_cntr][0]*cs_ee_t_ct_temp[ee_cntr] - ee_eff_l[4][i][0]*ee_cs[4][i] - ee_eff_l[4][j][0]*ee_cs[4][j])/cs_ee_t_ct[ee_cntr]
            ee_cntr += 1
    for i in range(len(path_interference_mumu)):
        for j in range(i+1, len(path_interference_mumu)):
            mumu_eff_t_ct[mumu_cntr][0] = (mumu_eff_t_ct_temp[mumu_cntr][0]*cs_mumu_t_ct_temp[mumu_cntr] - mumu_eff_l[4][i][0]*mumu_cs[4][i] - mumu_eff_l[4][j][0]*mumu_cs[4][j])/cs_mumu_t_ct[mumu_cntr]
            mumu_cntr += 1
    for i in range(len(path_interference_tautau)):
        for j in range(i+1, len(path_interference_tautau)):
            for tag_num in range(4):
                tautau_eff_t_ct[tautau_cntr][tag_num] = (tautau_eff_t_ct_temp[tautau_cntr][tag_num]*cs_tautau_t_ct_temp[tautau_cntr] - tautau_eff_l[4][i][tag_num]*tautau_cs[4][i] - tautau_eff_l[4][j][tag_num]*tautau_cs[4][j])/cs_tautau_t_ct[tautau_cntr]
            tautau_cntr += 1
    ee_lambdas_len = len(path_pureqcd_ee)
    mumu_lambdas_len = len(path_pureqcd_mumu)
    tautau_lambdas_len = len(path_pureqcd_tautau)
    return [ee_eff_l, mumu_eff_l, tautau_eff_l, ee_eff_t_ct, mumu_eff_t_ct, tautau_eff_t_ct, ee_lambdas_len, mumu_lambdas_len, tautau_lambdas_len]

# Branching Fraction helper functions:

def momentum(Mlq, mq, ml):
    a = (Mlq + ml)**2 - mq**2
    b = (Mlq - ml)**2 - mq**2
    return math.sqrt(a*b)/(2*Mlq)

def abseffcoupl_massfactor(Mlq, mq, ml):
    return Mlq**2 - (ml**2 + mq**2) - (ml**2 - mq**2)**2/Mlq**2 - (6*ml*mq)

def decay_width_massfactor(Mlq, M):
    #M is a list with [Mlq, mq, ml]
    return momentum(Mlq,M[0], M[1])*abseffcoupl_massfactor(Mlq,M[0],M[1])/(8 * math.pi**2 * Mlq**2)


mev2gev = 0.001

mass_quarks = {'1': [2.3 , 4.8], '2': [1275, 95], '3': [173070, 4180]}
mass_leptons= {'1': [0.511, 0.0022], '2': [105.7, 0.17], '3': [1777, 15.5]}

for gen in mass_quarks:
    mass_quarks[gen] = [x*mev2gev for x in mass_quarks[gen]]

for gen in mass_leptons:
    mass_leptons[gen] = [x*mev2gev for x in mass_leptons[gen]]

def make_mass_dict(ls, num_lam):
    md = {}
    for i in range(num_lam):
        if ls[i][4] == 'L':
            md[ls[i]] = [[mass_quarks[ls[i][2]][1], mass_leptons[ls[i][3]][0]], [mass_quarks[ls[i][2]][0], mass_leptons[ls[i][3]][1]]]
        elif ls[i][4] == 'R':
            md[ls[i]] = [[mass_quarks[ls[i][2]][1], mass_leptons[ls[i][3]][0]]]
    return md

# Calculate branching fraction of dielectron, dimuon an ditau using decay_width
def branching_fraction(all_ls, all_lam, md, Mlq, width_const = 0):
    denom = width_const
    numer = [0,0,0]
    for i in range(len(all_ls)):
        for j in range(len(all_ls[i])):
            denom += all_lam[i][j]**2 * decay_width_massfactor(Mlq, md[all_ls[i][j]][0])
            numer[i] += all_lam[i][j]**2 * decay_width_massfactor(Mlq, md[all_ls[i][j]][0])
            if all_ls[i][j][4] == 'L':
                denom += all_lam[i][j]**2 * decay_width_massfactor(Mlq, md[all_ls[i][j]][1])
    bf = [nu/denom for nu in numer]
    # print(f"Branching fractions: {bf}")
    return bf

# separate sympy lambda symbols and lambda strings
def get_lam_separate(lam):
    ee_lam = []
    mumu_lam = []
    tautau_lam = []
    ee_ls = []
    mumu_ls = []
    tautau_ls = []

    for lamda in lam:
        temp_str_sym = str(lamda)
        if temp_str_sym[3] == '1':
            ee_lam.append(lamda)
            ee_ls.append(temp_str_sym)
        elif temp_str_sym[3] == '2':
            mumu_lam.append(lamda)
            mumu_ls.append(temp_str_sym)
        elif temp_str_sym[3] == '3':
            tautau_lam.append(lamda)
            tautau_ls.append(temp_str_sym)
    return [ee_lam, mumu_lam, tautau_lam], [ee_ls, mumu_ls, tautau_ls]


# Calculate a partial polynomial
def get_chisq_ind(tag, mass, all_lam, cs_list, eff_list, b_frac, ignorePairSingle, margin):
    [ee_lam, mumu_lam, tautau_lam] = all_lam
    [ee_eff_l, mumu_eff_l, tautau_eff_l, ee_eff_t_ct, mumu_eff_t_ct, tautau_eff_t_ct, ee_lambdas_len, mumu_lambdas_len, tautau_lambdas_len] = eff_list
    #[ee_cs, mumu_cs, tautau_cs, cs_ee_t_ct, cs_mumu_t_ct, cs_tautau_t_ct] = cs_list
    [ee_cs, mumu_cs, tautau_cs, cs_ee_t_ct, cs_mumu_t_ct, cs_tautau_t_ct, cs_ee_t_ct_temp, cs_mumu_t_ct_temp, cs_tautau_t_ct_temp] = cs_list
    if (tag<4 and len(tautau_eff_l[0])==0) or (tag==4 and len(ee_eff_l[0])==0) or (tag==5 and len(mumu_eff_l[0])==0):
        return 0
    if tag<4:
        num_bin = len(tautau_eff_l[0][0][tag])
    elif tag == 4:
        num_bin = len(ee_eff_l[0][0][0])
    elif tag == 5:
        num_bin = len(mumu_eff_l[0][0][0])
    nq = [0.0]*num_bin
    np = [0.0]*num_bin
    ns = [0.0]*num_bin
    ni = [0.0]*num_bin
    nt = [0.0]*num_bin
    ntc= [0.0]*num_bin
    nsm= NSM[tag]
    nd = ND[tag]
    denominator = [nd[bin_no] + margin*margin*nd[bin_no]**2 for bin_no in range(num_bin)]
    if tag<4:
        nq = [nq[bin_no] + tautau_cs[0][0]*tautau_eff_l[0][0][tag][bin_no] * b_frac**2 *luminosity_tau for bin_no in range(num_bin)]
        for i in range(tautau_lambdas_len):
            np = [np[bin_no] + tautau_cs[1][i]*tautau_eff_l[1][i][tag][bin_no]*tautau_lam[i]**4 *b_frac**2 *luminosity_tau for bin_no in range(num_bin)]
            ns = [ns[bin_no] + tautau_cs[2][i]*tautau_eff_l[2][i][tag][bin_no]*tautau_lam[i]**2 *b_frac    *luminosity_tau for bin_no in range(num_bin)]
            ni = [ni[bin_no] + tautau_cs[3][i]*tautau_eff_l[3][i][tag][bin_no]*tautau_lam[i]**2               *luminosity_tau for bin_no in range(num_bin)]
            nt = [nt[bin_no] + tautau_cs[4][i]*tautau_eff_l[4][i][tag][bin_no]*tautau_lam[i]**4               *luminosity_tau for bin_no in range(num_bin)]
        ntc_cntr = 0
        for i in range(tautau_lambdas_len):
            for j in range(i+1, tautau_lambdas_len):
                "use cross-terms"
                ntc = [ntc[bin_no] + cs_tautau_t_ct[ntc_cntr]*tautau_eff_t_ct[ntc_cntr][tag][bin_no]* tautau_lam[i]**2 * tautau_lam[j]**2            *luminosity_tau for bin_no in range(num_bin)]
                ntc_cntr += 1
    elif tag==4:
        nq = [nq[bin_no] + ee_cs[0][0]*ee_eff_l[0][0][0][bin_no] * b_frac**2 *luminosity_e_mu for bin_no in range(num_bin)]
        for i in range(ee_lambdas_len):
            np = [np[bin_no] + ee_cs[1][i]*ee_eff_l[1][i][0][bin_no]*ee_lam[i]**4 *b_frac**2 *luminosity_e_mu for bin_no in range(num_bin)]
            ns = [ns[bin_no] + ee_cs[2][i]*ee_eff_l[2][i][0][bin_no]*ee_lam[i]**2 *b_frac    *luminosity_e_mu for bin_no in range(num_bin)]
            ni = [ni[bin_no] + ee_cs[3][i]*ee_eff_l[3][i][0][bin_no]*ee_lam[i]**2               *luminosity_e_mu for bin_no in range(num_bin)]
            nt = [nt[bin_no] + ee_cs[4][i]*ee_eff_l[4][i][0][bin_no]*ee_lam[i]**4               *luminosity_e_mu for bin_no in range(num_bin)]
        ntc_cntr = 0
        for i in range(ee_lambdas_len):
            for j in range(i+1, ee_lambdas_len):
                "use cross-terms"
                ntc = [ntc[bin_no] + cs_ee_t_ct[ntc_cntr]*ee_eff_t_ct[ntc_cntr][0][bin_no]* ee_lam[i]**2 * ee_lam[j]**2            *luminosity_e_mu for bin_no in range(num_bin)]
                ntc_cntr += 1
    elif tag==5:
        nq = [nq[bin_no] + mumu_cs[0][0]*mumu_eff_l[0][0][0][bin_no] * b_frac**2 *luminosity_e_mu for bin_no in range(num_bin)]
        for i in range(mumu_lambdas_len):
            np = [np[bin_no] + mumu_cs[1][i]*mumu_eff_l[1][i][0][bin_no]*mumu_lam[i]**4 *b_frac**2 *luminosity_e_mu for bin_no in range(num_bin)]
            ns = [ns[bin_no] + mumu_cs[2][i]*mumu_eff_l[2][i][0][bin_no]*mumu_lam[i]**2 *b_frac    *luminosity_e_mu for bin_no in range(num_bin)]
            ni = [ni[bin_no] + mumu_cs[3][i]*mumu_eff_l[3][i][0][bin_no]*mumu_lam[i]**2            *luminosity_e_mu for bin_no in range(num_bin)]
            nt = [nt[bin_no] + mumu_cs[4][i]*mumu_eff_l[4][i][0][bin_no]*mumu_lam[i]**4            *luminosity_e_mu for bin_no in range(num_bin)]
        ntc_cntr = 0
        for i in range(mumu_lambdas_len):
            for j in range(i+1, mumu_lambdas_len):
                "use cross-terms"
                ntc = [ntc[bin_no] + cs_mumu_t_ct[ntc_cntr]*mumu_eff_t_ct[ntc_cntr][0][bin_no]* mumu_lam[i]**2 * mumu_lam[j]**2            *luminosity_e_mu for bin_no in range(num_bin)]
                ntc_cntr += 1
    chi_ind = 0.0
    for bin_no in range(num_bin):
        if ignorePairSingle:
            chi_ind += ((nq[bin_no]+ni[bin_no]+nt[bin_no]+ntc[bin_no]+nsm[bin_no]-nd[bin_no])**2)/denominator[bin_no] 
        else:
            chi_ind += ((nq[bin_no]+np[bin_no]+ns[bin_no]+ni[bin_no]+nt[bin_no]+ntc[bin_no]+nsm[bin_no]-nd[bin_no])**2)/denominator[bin_no] 
    chi_sq_tag = sym.Add(chi_ind)
    return chi_sq_tag

# compute the polynomial by getting partial polynomials
def get_chi_square_symb(mass, all_lam, cs_list, eff_list, br_frac, ignorePairSingle, margin):    
    ee_chi = 0
    mumu_chi = 0
    hhbt_chi = 0
    hhbv_chi = 0
    lhbt_chi = 0
    lhbv_chi = 0
    if (len(all_lam[0]) > 0):
        ee_chi = sym.simplify(get_chisq_ind(4, mass, all_lam, cs_list, eff_list, br_frac[0], ignorePairSingle, margin))
        print("Dielectron contributions computed.")
    if (len(all_lam[1]) > 0):
        mumu_chi = sym.simplify(get_chisq_ind(5, mass, all_lam, cs_list, eff_list, br_frac[1], ignorePairSingle, margin))
        print("Dimuon contributions computed.")
    if (len(all_lam[2]) > 0):
        hhbt_chi = sym.simplify(get_chisq_ind(0, mass, all_lam, cs_list, eff_list, br_frac[2], ignorePairSingle, margin))
        print("Ditau: HHbT contributions computed.")
        hhbv_chi = sym.simplify(get_chisq_ind(1, mass, all_lam, cs_list, eff_list, br_frac[2], ignorePairSingle, margin))
        print("Ditau: HHbV contributions computed.")
        lhbt_chi = sym.simplify(get_chisq_ind(2, mass, all_lam, cs_list, eff_list, br_frac[2], ignorePairSingle, margin))    
        print("Ditau: LHbT contributions computed.")
        lhbv_chi = sym.simplify(get_chisq_ind(3, mass, all_lam, cs_list, eff_list, br_frac[2], ignorePairSingle, margin))
        print("Ditau: LHbV contributions computed.")
    return sym.Add(ee_chi,mumu_chi,hhbt_chi, hhbv_chi, lhbt_chi, lhbv_chi)


# Use the lambdified function (numpy_chisq) to calculate chi square for the given query input
def get_delta_chisq(lam_vals, lam_vals_original, chisq_min, numpy_chisq, num_lam):
    validity_list = []
    delta_chisq = []
    for lam_val,lam_val_copy in zip(lam_vals,lam_vals_original):
        temp = [float(x) for x in lam_val]
        all_zeroes = True
        for x in temp:
            if x:
                all_zeroes = False
                break
        if all_zeroes:
            temp = [0.00000001]*len(lam_val)
        chisq_given_vals = numpy_chisq(*flatten(temp))
        delta_chisq.append(chisq_given_vals - chisq_min)
        if chisq_given_vals - chisq_min <= chi_sq_limits[num_lam-1]:
            validity_list.append("Yes")
        else:
            validity_list.append("No")
    return delta_chisq, validity_list

# print red/cyan
def prCyan(skk): print("\033[96m{}\033[00m" .format(skk), end="")
def prRed(skk): print("\033[91m{}\033[00m" .format(skk), end="")

# Initiate procedure if non-interactive
def initiate_with_files(card, vals, output_yes, output_no):
    try:
        with open(card) as c:
            c_lines = c.readlines()
    except OSError:
        sys.exit(f"Card file {card} does not exist. Exiting.")
    try:
        with open(vals) as v:
            v_lines = v.readlines()
    except OSError:
        sys.exit(f"Values file {vals} does not exist. Exiting.")
    if (len(c_lines) < 5):
        sys.exit(f"Number of lines in file: {len(c_lines)}, expected 5. Exiting.")
    mass_f = c_lines[0].split('#')[0].strip()
    lambdas_f = c_lines[1].split('#')[0].strip()
    ignore_f = c_lines[2].split('#')[0].strip()
    sigma_f = c_lines[3].split('#')[0].strip()
    margin_f = c_lines[4].split('#')[0].strip()
    global chi_sq_limits
    if (sigma_f=="1"):
        chi_sq_limits = chi_sq_limits_1
    elif (sigma_f=="2"):
        chi_sq_limits = chi_sq_limits_2
    else:
        sys.exit("[Sigma Error]: Line 4 of input card must contain either 1 or 2 as the sigma value. Exiting.")
    if not (ready_to_initiate(mass_f, lambdas_f, ignore_f, margin_f)):
        sys.exit("[Input Error]: Syntax Error encountered in input card. Exiting.")
    home(*parse(mass_f, lambdas_f, ignore_f, margin_f, v_lines), False, output_yes, output_no)

# Checking if the inputs are in the correct format
def ready_to_initiate(mass_f, lambdas_f, ignore_f, margin_f):
    try:
        mass = int(mass_f)
        if (mass < 1000 or mass > 3000):
            prRed("[Mass Error]: Acceptable mass values in GeV: integers between 1000 and 3000.\n")
            return False
    except ValueError:
        prRed("[Mass Error]: Mass value should be an integer (between 1000 and 3000).\n")
        return False
    try:
        margin = float(margin_f)
        if (margin < 0 or margin > 1):
            prRed("[Systematic-Error Error]: Acceptable systematic error values: float values between 0 and 1.\n")
            return False
    except ValueError:
        prRed("[Systematic-Error Error]: Systematic error value should be a float (between 0 and 1).\n")
        return False
    l = lambdas_f.split()
    if len(l)==0:
        return False
    for i in range(len(l)):
        if (l[i][:2]!="LM"):
            return False
        if not (l[i][2]=='1' or l[i][2]=='2' or l[i][2]=='3'):
            return False
        if not (l[i][3]=='1' or l[i][3]=='2' or l[i][3]=='3'):
            return False
        if (l[i][4]!='L' and l[i][4]!='R'):
            return False
    if not (ignore_f.lower()=="yes" or ignore_f.lower()=="no" or ignore_f.lower()=="n" or ignore_f.lower()=="y"):
        prRed("ignore_single_pair takes input either 'yes'/'y' or 'no'/'n'\n")
        return False
    return True

# Check if queries are in correct form
def lam_val_ok(lam_val_f, num_lam):
    lam_vals = lam_val_f.replace(',', ' ').split()
    if (len(lam_vals) != num_lam):
        return False
    try:
        for i in range(num_lam):
            a = float(lam_vals[i])
    except ValueError:
        return False
    return True

# Initiate procedure if interactive
def initiate_interactive():
    mass_f = ""
    lambdas_f = ""
    ignore_f = "yes"
    sigma_limit = 2
    margin_f = "0.1"
    lam_values_f = []
    print("Commands available: ", end="")
    prCyan("mass=, couplings=, systematic_error=, ignore_single_pair=(yes/no), significance=(1/2), import_model=, status, initiate, help\n")
    print("Default Model loaded: U1")
    print("Couplings available: LM11L, LM12L, LM13L, LM21L, LM22L, LM23L, LM31L, LM32L, LM33L, LM11R, LM12R, LM13R, LM21R, LM22R, LM23R, LM31R, LM32R, LM33R")
    global chi_sq_limits
    while(True):
        prCyan("icalq > ")
        s = input().split("=")
        slen = len(s)
        if (s[0].strip()=="mass" and slen==2):
            mass_f = s[1].strip()
        elif (s[0].strip()==""):
            continue
        elif (s[0].strip()=="couplings" and slen>1):
            lambdas_f = s[1].strip()
        elif (s[0].strip()=="ignore_single_pair" and slen==2):
            ignore_f = s[1].strip()
        elif (s[0].strip()=="systematic-error" and slen==2):
            margin_f = s[1].strip()
        elif (s[0].strip()=="import_model" and slen==2):
            print("Currently only U1 model is available.")
        elif (s[0].strip()=="significance" and slen==2):
            if (s[1].strip()=="1"):
                sigma_limit = 1
                chi_sq_limits = chi_sq_limits_1
            elif (s[1].strip()=="2"):
                sigma_limit = 2
                chi_sq_limits = chi_sq_limits_2
            else:
                prRed("Allowed values of 'significance': 1 or 2\n")
        elif (s[0].strip()=="status"):
            print(f"\nMass: {mass_f}\nCouplings: {lambdas_f}\nIgnore Single & Pair = {ignore_f}\nSignificance = {sigma_limit}\nSystematic-Error = {margin_f}")
        elif (s[0].strip()=="help"):
            print("Commands available: ",end="")
            prCyan("mass=, couplings=, systematic-error=, ignore_single_pair=(yes/no), significance=(1/2), import_model=, status, initiate, help\n")
            print("Couplings available: LM11L, LM12L, LM13L, LM21L, LM22L, LM23L, LM31L, LM32L, LM33L, LM11R, LM12R, LM13R, LM21R, LM22R, LM23R, LM31R, LM32R, LM33R")
            print("commands with '=' expect appropriate value. Read README.md for more info on individual commands.\n")
        elif (s[0].strip()=="initiate"):
            if not ready_to_initiate(mass_f, lambdas_f, ignore_f, margin_f):
                prRed("[Lambda Error]: Example of a valid input - 'LM23L LM33R LM12R'\n")
                continue
            num_lam = len(lambdas_f.split())
            lam_values_f = [" ".join(["0"]*num_lam)]
            home(*parse(mass_f, lambdas_f, ignore_f, margin_f, lam_values_f), True)
        elif (s[0].strip()=="exit" or s[0].strip()=="q" or s[0].strip()=="exit()" or s[0].strip()==".exit"):
            return
        else:
            prRed(f"Command {s[0]} not recognised. Please retry or enter 'q' to exit.\n")

def home(mass, lambdastring, original_lambdastring, ignorePairSingle, lam_vals, original_lam_vals, margin, interactive, output_yes="icalq_yes.csv", output_no="icalq_no.csv"):
    num_lam = len(lambdastring)
    lam = sym.symbols(lambdastring)
    cs_list = get_cs(mass, lambdastring, num_lam)
    closest_mass = get_closest_mass(mass)
    eff_list = get_efficiencies(closest_mass, lambdastring, num_lam, cs_list)
    mass_dict = make_mass_dict(lambdastring, num_lam)
    all_lam, all_ls = get_lam_separate(lam)
    br_frac = branching_fraction(all_ls, all_lam, mass_dict, mass, 0)
    chisq_symb = get_chi_square_symb(mass, all_lam, cs_list, eff_list, br_frac, ignorePairSingle, margin)
    # print("Lambdifying...")
    numpy_chisq=lambdify(flatten(lam),chisq_symb, modules='numpy')
    startLambda = 0.5
    startLambdas = np.array([startLambda for x in range(num_lam)])
    print("Minimizing...")
    minima = optimize.minimize(lambda x: numpy_chisq(*flatten(x)),startLambdas, method='Nelder-Mead',options={'fatol':0.0001})
    minima_list_1 = [optimize.minimize(lambda x: numpy_chisq(*flatten(x)),randarr, method='Nelder-Mead',options={'fatol':0.0001}) for randarr in np.random.rand(3, num_lam)]
    minima_list_2 = [optimize.minimize(lambda x: numpy_chisq(*flatten(x)),randarr) for randarr in np.random.rand(3, num_lam)]
    for m in minima_list_1 + minima_list_2:
        if m.fun < minima.fun:
            print("New Minimum Found!")
            minima = m
    chisq_min = minima.fun
    opt_lambda_list = minima.x
    print("Minimum Chi Square at values:", end="")
    print(*[f"\n{lambdastring[i]} : {opt_lambda_list[i]}" for i in range(num_lam)])
    if (interactive):
        print("Input query values in the following order: " ,end='\t')
        for x in original_lambdastring:
            print(x, end = '\t')
        while(True):
            print("\n > ", end="")
            lam_val_f = input()
            if (lam_val_f == "done"):
                return
            if not lam_val_ok(lam_val_f, num_lam):
                prRed(f"[Query Error]: Please input {num_lam} float input/s.\n")
                # prRed(f"[Query Error]: Query values for lambdas are either not {num_lam} (number of lambdas) in count or not convertible to float.\n")
                print("Type 'done' to continue to icalq prompt.")
                continue
            lam_vals, original_lam_vals = parse_lam(original_lambdastring, lam_val_f)
            delta_chisq, validity_list = get_delta_chisq(lam_vals, original_lam_vals, chisq_min, numpy_chisq, num_lam)
            print(f"Delta Chi Square: {delta_chisq[0]}\nAllowed: {validity_list[0]}")
            if (delta_chisq[0] < 0):
                print("A negative value should imply precision less than 1e-4 while calculating minima and can be considered equal to 0. Try initiating again to randomize minimization.")
    delta_chisq, validity_list = get_delta_chisq(lam_vals, original_lam_vals, chisq_min, numpy_chisq, num_lam)
    yes_list = [i for i in range(len(validity_list)) if validity_list[i]=='Yes']
    no_list = [i for i in range(len(validity_list)) if validity_list[i]=='No']
    print("\nYes List:\n")
    with open(output_yes, "w") as f:
        for x in original_lambdastring:
            print(x, end='\t\t')
            print(x, end=',', file=f)
        print("Delta_chisquare")
        print("Delta_chisquare", file=f)
        for i in yes_list:
            for x in original_lam_vals[i]:
                print(x, end="\t\t")
                print(x, end=",", file=f)
            print(delta_chisq[i])
            print(delta_chisq[i], file=f)
    print("\nNo List:\n")
    with open(output_no, "w") as f:
        for x in original_lambdastring:
            print(x, end='\t\t')
            print(x, end=',', file=f)
        print("Delta_chisquare")
        print("Delta_chisquare", file=f)
        for i in no_list:
            for x in original_lam_vals[i]:
                print(x, end='\t\t')
                print(x, end=',', file=f)
            print(delta_chisq[i])
            print(delta_chisq[i], file=f)
    print(f"Output files {output_yes} and {output_no} written")


def welcome_message():
    try:
        with open("banner.txt") as f:
            contents = f.read()
            print("\n")
            print(contents)
            # print("LHC Dilepton Limits Calculator")
    except OSError:
        print("\niCaLQ: Indirect LHC-Limits Calculator for Leptoquark models\nAlpha version\n")

def help_message():
    print("iCaLQ Usage:")
    prCyan("\t\t\t > python3 icalq.py [options]\n")
    print("Options:")
    print("--help\t\t\t : Display this help message.")
    print("--input-card=[filename]\t : Input card file. Format is explained in README.txt")
    print("--input-values=[filename]: Input values to check from the given file. Format is explained in README.txt")
    print("--non-interactive, -ni\t : Run in non-interactive mode. This requires input-card and input-values to be specified")
    print("--no-banner, -nb\t : iCaLQ banner is not printed.")
    print("--output-yes=[filename]\t : Specify the name of output file (allowed values) (overwrites the existing file). Default: icalq_yes.csv")
    print("--output-no=[filename]\t : Specify the name of output file (disallowed values) (overwrites the existing file). Default: icalq_no.csv")
    print(("\n"))


# Starting function, accepts command line argument and passes control to further functions accordingly.
def run():
    n = len(sys.argv)
    input_card = ""
    input_vals = ""
    help_user = False
    non_interactive = False
    no_banner = False
    output_yes = "icalq_yes.csv"
    output_no = "icalq_no.csv"
    for i in range(1,n):
        if (sys.argv[i][:13]=="--input-card="):
            input_card = sys.argv[i][13:]
        elif (sys.argv[i][:15]=="--input-values="):
            input_vals = sys.argv[i][15:]
        elif (sys.argv[i]=="--help"):
            help_user = True
        elif (sys.argv[i]=="--non-interactive" or sys.argv[i] == "-ni"):
            non_interactive = True
        elif (sys.argv[i]=="--no-banner" or sys.argv[i]=="-nb"):
            no_banner = True
        elif (sys.argv[i][:13]=="--output-yes="):
            output_yes = sys.argv[i][13:]
        elif (sys.argv[i][:12]=="--output-no="):
            output_no = sys.argv[i][12:]
        else:
            prRed(f"Argument {sys.argv[i]} not recognised. Continuing without it.\n")
    
    if (not no_banner):
        welcome_message()
    
    if help_user:
        help_message()
        sys.exit(0)

    if (non_interactive):
        if (input_card==""):
            sys.exit("[Card Error]: Input Card file not specified in the expected format (mandatory for non-interactive mode). Exiting.\n")
        if (input_vals==""):
            sys.exit("[Values Error]: Input Values file not specified in the expected format (mandatory for non-interactive mode). Exiting.\n")
        print(f"Input Card file: {input_card}")            
        print(f"Input Values file: {input_vals}")
        print(f"Output Yes file: {output_yes}")
        print(f"Output No file: {output_no}")
        initiate_with_files(input_card, input_vals, output_yes, output_no)
    else:
        initiate_interactive()            


try:
    run()
except KeyboardInterrupt:
    sys.exit("\n\nKeyboardInterrupt recieved. Exiting.")



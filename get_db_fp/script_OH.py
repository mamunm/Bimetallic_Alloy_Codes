#!/usr/bin/env python
# -*-coding: utf-8 -*-

#generate_reconstruction_fp.py
#Osman Mamun
#LAST UPDATED: 09-09-2018

import numpy as np
from ase.db import connect
from catkit.mydb import FingerprintDB
import os
from code_gen_struct import get_struct
from catkit.pawprint import Fingerprinter
from catlearn.featurize.adsorbate_prep import autogen_info
from catlearn.featurize.setup import FeatureGenerator
import re
import ast
from glob import glob 
import sys

def get_sym(A):
    ordered_metals = ['Y', 'La', 'Sc', 'Zr', 'Hf', 'Ti', 'Ta', 'Nb', 'V',
                      'Cr', 'Mo', 'W', 'Re', 'Tc', 'Os', 'Ru', 'Ir','Rh',
                      'Ni', 'Co', 'Fe', 'Mn', 'Pt', 'Pd', 'Au', 'Ag', 'Cu',
                      'Zn',  'Cd', 'Hg', 'Al',  'Ga',  'In', 'Tl', 'Pb',
                      'Sn', 'Bi']
    if '3' in A:
        return A, 'L12'
    elif A in ordered_metals:
        return A, 'A1'
    else:
        sym = []
        for i in ordered_metals:
            if i in A:
                sym += [i]

        if not len(sym) == 2:
            A = 'ERROR!!'
        else:
            return '2'.join(sym) + '2', 'L10'

def get_catlearn_fp(A, Name=False):

    if not Name:
        A = autogen_info([A])
    gen = FeatureGenerator()
    if isinstance(A, list):
        gen._get_atom_types(A)
    else:
        gen._get_atom_types([A])

    train_fpv = [gen.mean_chemisorbed_atoms,
                 gen.bulk,
                 gen.term,
                 gen.strain,
                 gen.mean_surf_ligands,
                 gen.mean_site,
                 gen.median_site,
                 gen.max_site,
                 gen.min_site,
                 gen.sum_site,
                 #gen.generalized_cn,
                 #gen.bag_cn,
                 #gen.en_difference_ads,
                 #gen.en_difference_chemi,
                 #gen.en_difference_active,
                 #gen.bag_atoms_ads,
                 #gen.bag_connections_ads,
                 gen.count_chemisorbed_fragment]
    
    if not Name:
        gen.normalize_features(A)
        return gen.return_vec(A, train_fpv)
    else:
        return gen.return_names(train_fpv)

#DONE: N O H C S CH CH2  CH3 NH OH
data = ['OH']
for d in data:
    txt_f = glob('{}_*.txt'.format(d))
    db_f = glob('{}_*.db'.format(d))
    for fff in txt_f + db_f:
        if os.path.exists(fff):
            os.remove(fff)

parameters_list_nonlocal = []
operation_list_nonlocal = []
parameters_list_local = []
operation_list_local = []

#peroidic fp for slabs only
parameters_list_nonlocal += [[
                        'atomic_number',
                        'atomic_radius',
                        'dband_center_slab',
                        'dband_width_slab',
                        'dband_skewness_slab',
                        'dband_kurtosis_slab',
                        'dipole_polarizability',
                        'electron_affinity',
                        'heat_of_formation',
                        'specific_heat'
                                ]]

parameters_list_nonlocal += [parameters_list_nonlocal[-1]]

parameters_list_nonlocal += [[
                        'atomic_number',
                        'atomic_radius',
                        'dband_center_slab',
                        'dband_width_slab',
                        'dipole_polarizability',
                        'electron_affinity',
                        'heat_of_formation',
                        'specific_heat',
                        'en_allen'
                                ]]

operation_list_nonlocal += [
                      'periodic_convolution',
                      ['periodic_convolution', {'d': 1}],
                      'bimetal_fp'
                      ]

#fp for ads_slab_system
parameters_list_local += [[
                       'atomic_number',
                       'atomic_radius',
                       'dband_center_slab',
                       'dband_width_slab',
                       'boiling_point',
                       'dipole_polarizability',
                       'electron_affinity',
                       'en_allen',
                       'melting_point',
                       'heat_of_formation',
                       'vdw_radius'
                             ]]

parameters_list_local += [parameters_list_local[-1]]
parameters_list_local += [parameters_list_local[-1]]

operation_list_local += [
                  'bonding_convolution',
                  'layered_sum',
                  'local_ads_metal_fp'
                 ]


for datum in data:
    #generate atoms object
    #get data obtained from cathub
    w_data = np.load('../get_data_cathub/{}_data.npy'.format(datum))[()]
    #connect the database to store the objects
    db = connect('{}_atoms.db'.format(datum))
    fp_file = '{}_fp.db'.format(datum)

    with FingerprintDB(fp_file) as fpd:
        for i in range(2):
            fpd.parameter_entry('AN_PC{}'.format(i),
                                'Atomic Number Slab: pc {}'.format(i))
            fpd.parameter_entry('AR_PC{}'.format(i),
                                'Atomic Radius Slab: pc {}'.format(i))
            fpd.parameter_entry('DBC_PC{}'.format(i),
                                'Dband center Slab: pc {}'.format(i))
            fpd.parameter_entry('DBW_PC{}'.format(i),
                                'Dband Width Slab: pc {}'.format(i))
            fpd.parameter_entry('DBS_PC{}'.format(i),
                                'Dband Skewness Slab: pc {}'.format(i))
            fpd.parameter_entry('DBK_PC{}'.format(i),
                                'Dband Kurtosis Slab: pc {}'.format(i))
            fpd.parameter_entry('DP_PC{}'.format(i),
                                'Dipole Polarizability Slab: pc {}'.format(i))
            fpd.parameter_entry('EA_PC{}'.format(i),
                                'Electron Affinity Slab: pc {}'.format(i))
            fpd.parameter_entry('HOF_PC{}'.format(i),
                                'Heat Of Formation Slab: pc {}'.format(i))
            fpd.parameter_entry('SH_PC{}'.format(i),
                                'Specific Heat Slab: pc {}'.format(i))

        fpd.parameter_entry('AN_BM{}'.format(i),
                            'Atomic Number Slab: bm {}'.format(i))
        fpd.parameter_entry('AR_BM{}'.format(i),
                            'Atomic Radius Slab: bm {}'.format(i))
        fpd.parameter_entry('DBC_BM{}'.format(i),
                            'Dband center Slab: bm {}'.format(i))
        fpd.parameter_entry('DBW_BM{}'.format(i),
                            'Dband Width Slab: bm {}'.format(i))
        fpd.parameter_entry('DP_BM{}'.format(i),
                            'Dipole Polarizability Slab: bm {}'.format(i))
        fpd.parameter_entry('EA_BM{}'.format(i),
                            'Electron Affinity Slab: bm {}'.format(i))
        fpd.parameter_entry('HOF_BM{}'.format(i),
                            'Heat Of Formation Slab: bm {}'.format(i))
        fpd.parameter_entry('SH_BM{}'.format(i),
                            'Specific Heat Slab: bm {}'.format(i))
        fpd.parameter_entry('EN_BM{}'.format(i),
                            'Electronegativity: bm {}'.format(i))

        if datum in ['C', 'H', 'N', 'S', 'O']:
            tempfl = ['BC', 'LS_1', 'LS_2', 'LS_3', 'LS_4', 'LAMFP']
        else:
            tempfl = ['BC', 'LS_1', 'LS_2', 'LS_3', 'LS_4', 'LS_5', 'LAMFP']

        for i in tempfl:
            fpd.parameter_entry('AN_{}'.format(i),
                                'Atomic Number Ads Slab: {}'.format(i))
            fpd.parameter_entry('AR_{}'.format(i),
                                'Atomic Radius Ads Slab: {}'.format(i))
            fpd.parameter_entry('DBC_{}'.format(i),
                                'Dband center Ads Slab: {}'.format(i))
            fpd.parameter_entry('DBW_{}'.format(i),
                                'Dband width Ads Slab: {}'.format(i))
            fpd.parameter_entry('BP_{}'.format(i),
                                'Boiling Point: {}'.format(i))
            fpd.parameter_entry('DP_{}'.format(i),
                                'Dipole Polarizability Ads Slab: {}'.format(i))
            fpd.parameter_entry('EA_{}'.format(i),
                                'Electron Affinity Ads Slab: {}'.format(i))
            fpd.parameter_entry('EN_{}'.format(i),
                                'Electronegativity Ads Slab: {}'.format(i))
            fpd.parameter_entry('MP_{}'.format(i),
                                'Melting Point Ads Slab: {}'.format(i))
            fpd.parameter_entry('HOF_{}'.format(i),
                                'Heat Of Foramtion Ads Slab: {}'.format(i))
            fpd.parameter_entry('VDW_{}'.format(i),
                                'VDW Radius Ads Slab: {}'.format(i))
        '''
        #Catlearn generated fingerprints
        d_slab, d_structure = get_struct(symbol='Pt',
                                     adsorbate='H',
                                     site='top',
                                     site_type='Pt',
                                     SB_symbol='L12')
        CLnames = get_catlearn_fp(d_structure, Name=True)
        for n, cln in enumerate(CLnames):
            a = ''.join([i[0] for i in cln.split('_')])
            a += a + str(n)
            fpd.parameter_entry(a, cln)
        '''

        fpd.parameter_entry('E*', 'Adsorption energy [eV]')

        fpd.metadata_params_entry('Adsorbate', 'Adsorbate_name')
        fpd.metadata_params_entry('Slab', 'slab_name')
        fpd.metadata_params_entry('Site', 'Site_name')
        fpd.metadata_params_entry('Site_type', 'Site_type_name')
        fpd.metadata_params_entry('sb_symbol', 'Strukterbericht symbol')

        par = fpd.get_parameters()
        met_par = fpd.get_metadata_params()

        for i, d in enumerate(w_data['reactions']['edges']):

            if datum in ['N', 'H']:
                if not '0.5' in d['node']['reactants']:
                    continue
            #sdict = ast.literal_eval(d['node']['sites'])
            #adsorbate = list(dict.keys())[0]
            #site, site_type = list(dict.values())[0].split('|', 1)
            adsorbate = re.findall(r'"([A-Z0-4]+)":', d['node']['sites'])[0]
            reactionenergy = d['node']['reactionEnergy']
            symbol, sb_symbol = get_sym(d['node']['surfaceComposition'])
            if datum == 'S':
                if '"H":' in d['node']['sites']:
                    print('Problematic stuffs!\n\n')
                    print(d['node']['sites'])
                    continue 
            site, site_type = re.findall(r'(?<=: )"(.+)"',
                                         d['node']['sites'])[0].split('|', 1)

            if site not in ['top', 'bridge', 'hollow']:
                print('Site {} not a valid site_type.'.format(site))
                continue

            if sb_symbol == 'A1':
                site_type = site_type.replace('A', symbol)
            elif sb_symbol == 'L12':
                a, b = symbol.split('3')
                site_type = site_type.replace('A', 'aaa').replace('B', 'bbb')
                site_type = site_type.replace('aaa', a).replace('bbb', b)
            else:
                a, b, _ = symbol.split('2')
                if a > b:
                    a, b = b, a
                site_type = site_type.replace('A', 'aaa').replace('B', 'bbb')
                site_type = site_type.replace('aaa', a).replace('bbb', b)

            with open('{}_debug.txt'.format(datum), 'a') as ff:
                ff.write('Params: {} {} {} {}\n'.format(symbol,
                                                      adsorbate,
                                                      site,
                                                      site_type))
            if ':' in site_type:
                print('This is a problematic site type.\n\n\n')
                continue
            slab, structure = get_struct(symbol=symbol,
                                         adsorbate=adsorbate,
                                         site=site,
                                         site_type=site_type,
                                         SB_symbol=sb_symbol)

            if slab is None or structure is None:
                with open('debug.txt', 'a') as f:
                    f.write('Problem found in the file. Specification:\n')
                    f.write('symbol: {}\n'.format(symbol))
                    f.write('adsorbate: {}\n'.format(adsorbate))
                    f.write('site: {}\n'.format(site))
                    f.write('site_type: {}\n'.format(site_type))
                continue
            #write slab to db
            db.write(slab,
                     symbol=symbol,
                     sb_symbol=sb_symbol,
                     site=site,
                     site_type=site_type,
                     adsorbate=adsorbate,
                     reaction_energy=reactionenergy,
                     grepid=i)

            #write structure to db
            db.write(structure,
                     symbol=symbol,
                     sb_symbol=sb_symbol,
                     site=site,
                     site_type=site_type,
                     adsorbate=adsorbate,
                     reaction_energy=reactionenergy,
                     grepid=i)

            fpd.image_entry(i)
            fpd.fingerprint_entry(i, par[-1], reactionenergy)

            fpd.metadata_entry(i, met_par[0], adsorbate)
            fpd.metadata_entry(i, met_par[1], symbol)
            fpd.metadata_entry(i, met_par[2], site)
            fpd.metadata_entry(i, met_par[3], site_type)
            fpd.metadata_entry(i, met_par[4], sb_symbol)

            fp_nonlocal = Fingerprinter(slab)
            fingerprints_nonlocal = fp_nonlocal.get_fp(parameters_list_nonlocal,
                                                       operation_list_nonlocal)
            fp_local = Fingerprinter(structure)
            fingerprints_local = fp_local.get_fp(parameters_list_local,
                                                 operation_list_local)
            '''
            try:
                catlearn_fp = get_catlearn_fp(structure)
            except AttributeError as e:
                print('Got problem in catlearn.')
                catlearn_fp = 'I\'m outta here.'
            '''
            if fingerprints_local == 'I\'m outta here.':
                print('Got problem in local.')
                continue
            if fingerprints_nonlocal == 'I\'m outta here.':
                print('Got problem in nonlocal.')
                continue
            '''
            if catlearn_fp == 'I\'m outta here.':
                continue
            '''
            for j, n in enumerate(fingerprints_nonlocal.tolist()[0]):
                fpd.fingerprint_entry(i, par[j], n)

            ind_cnt = len(fingerprints_nonlocal.tolist()[0])

            for j, n in enumerate(fingerprints_local.tolist()[0]):
                fpd.fingerprint_entry(i, par[ind_cnt+j], n)

            with open('{}_record.txt'.format(datum), 'a+') as fw:
                fw.write('Data entry id:{}\n'.format(i))
                fw.write('Metadata: ads: {} sym: {} site: {} stype: {}\n'.format(
                    adsorbate, symbol, site, site_type))
                fw.write('Number of parameter: {}\n'.format(len(par)))
                fw.write('Number of nonlocal fp: {}\n'.format(
                     len(fingerprints_nonlocal.tolist()[0])))
                fw.write('Number of local fp: {}\n'.format(
                     len(fingerprints_local.tolist()[0])))
                '''
                fw.write('Number of catlearn fp: {}'.format(
                     len(catlearn_fp.tolist()[0])))
                '''
            '''
            ind_cnt += len(fingerprints_local.tolist()[0])
            for j, n in enumerate(catlearn_fp.tolist()[0]):
                fpd.fingerprint_entry(i, par[ind_cnt+j], n)
            '''



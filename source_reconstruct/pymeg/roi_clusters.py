'''
This module contains definitions of clusters of areas that
are of interest for the analysis.

The top part of this module defines standard clusters, the
bottom contains functions to average over clusters.

This module also contains scripts to generate an image for
clusters of labels.
'''

import os
if 'DISPLAY' in os.environ.keys():
    try:
        from surfer import Brain
    except:
        Brain = None
        print('No pysurfer support')

import numpy as np
import pandas as pd

from conf_analysis.behavior import metadata
from joblib import Memory

memory = Memory(cachedir=os.environ['PYMEG_CACHE_DIR'], verbose=0)

# Define visual field clusters:
#   keys are cluster names, values uniquely identify ROIs in Wang
#   and Destrieux atlas.

visual_field_clusters = {
    'vfcvisual': (u'wang2015atlas.V1d', u'wang2015atlas.V1v',
                  u'wang2015atlas.V2d', u'wang2015atlas.V2v',
                  u'wang2015atlas.V3d', u'wang2015atlas.V3v',
                  u'wang2015atlas.hV4'),
    'vfcVO': (u'wang2015atlas.VO1', u'wang2015atlas.VO2',),
    'vfcPHC': (u'wang2015atlas.PHC1', u'wang2015atlas.PHC2'),
    'vfcV3ab': (u'wang2015atlas.V3A', u'wang2015atlas.V3B'),
    'vfcTO': (u'wang2015atlas.TO1', u'wang2015atlas.TO2'),
    'vfcLO': (u'wang2015atlas.LO1', u'wang2015atlas.LO2'),
    'vfcIPS_occ': (u'wang2015atlas.IPS0', u'wang2015atlas.IPS1'),
    'vfcIPS_dorsal': (u'wang2015atlas.IPS2', u'wang2015atlas.IPS3',
                      u'wang2015atlas.IPS4', u'wang2015atlas.IPS5'),
    'vfcSPL': (u'wang2015atlas.SPL1',),
    'vfcFEF': (u'wang2015atlas.FEF')
}


# ROIs from JW de Gee's eLife paper
jwrois = {'IPS_Pces': (u'JWDG.lr_IPS_PCes',),
          'M1': (u'JWDG.lr_M1',),
          'aIPS1': (u'JWDG.lr_aIPS1',)}

# Various frontal areas.
frontal = {'ACC': ('G&S_cingul-Ant-lh', 'G&S_cingul-Mid-Ant'),
           'frontomargin': ('G&S_frontomargin-',),
           'frontopol': ('G&S_transv_frontopol'),
           'f_inf_opercular': ('G_front_inf-Opercular',),
           'f_inf_orbital': ('G_front_inf-Orbital',),
           'f_inf_Triangul': ('G_front_inf-Triangul',),
           'Gf_middle': ('G_front_middle',),
           'Gf_sup': ('G_front_sup',),
           'Sf_inf': ('S_front_inf',),
           'Sf_middle': ('S_front_middle',),
           'Sf_sup': ('S_front_sup',)}

glasser = {'LIP': ('LIPd', 'LIPv'),
           'Area6_dorsal_medial': ('6mp', '6d'),
           'Area6_anterior': ('6ma', '6a'),
           'A10': ('10pp', '10r', '10v'),
           'A6si': ('s6-8', 'i6-8',)}
for A in ['PEF', '55b', '8Av', '8Ad', '9p', '8Bl', '8C', 'p9-46v', '46', '9-46d', '9a',
          'p10p', 'a10p', 'a47r', 'p32', 's32', 'a24', '9m', 'd32', 'a32pr', '8BM',
          'p24', 'a24', 'p32pr', '24dv', 'p24pr']:
    glasser[A] = A


rows = {'visual': visual_field_clusters,
        'choice': jwrois,
        'frontal': frontal}

all_clusters = {}
all_clusters.update(visual_field_clusters)
all_clusters.update(jwrois)
all_clusters.update(glasser)


def labels_to_clusters(labels, clusters, hemi='lh'):
    '''
    Sort MNE label files into clusters.
    '''
    rois2clusters = {}
    for name, areas in clusters.items():
        rois2clusters[name] = []
        for label in labels:
            for area in ensure_iter(areas):
                if hemi in label.name:
                    if area not in label.name:
                        continue
                    if ('wang' in label.name) or ('JWDG' in label.name):
                        rois2clusters[name].append(label)
                    elif ('L_'+area+'_ROI') in label.name:
                        rois2clusters[name].append(label)
                    elif ('R_'+area+'_ROI') in label.name:
                        rois2clusters[name].append(label)
    return rois2clusters


def rh(columns):
    '''
    Returns only right hemisphere column nanmes

    Arguments
    ---------
      columns : list of column names

    Returns
    -------
      List of right hemisphere column names
    '''
    return [x for x in columns if (x.startswith('rh') | x.endswith('rh'))]


def lh(columns):
    '''
    Returns only left hemisphere column nanmes

    Arguments
    ---------
      columns : list of column names

    Returns
    -------
      List of left hemisphere column names
    '''
    return [x for x in columns if (x.startswith('lh') | x.endswith('lh'))]


def filter_cols(columns, select):
    '''
    Filter columns such that they contain a keyword in select.


    Arguments
    ---------
      columns : list of column names
      select : list of str selectors

    Returns
    -------
      List of names (from columns) that contain at least one of the 
      strings in select

    '''
    return [x for x in columns if any([y.lower() in x.lower() for y in ensure_iter(select)])]


def reduce(df, all_clusters=all_clusters):
    '''
    Reduce ROIs to visual field clusters

    Arguments
    ---------
      df : DataFrame
        DataFrame that has areas as columns.
      all_clusters : dict
        This dict defines visual field clusters. It contains
        cluster names as keys and a list of column names of df
        that belong to a cluster as values. Entries in values
        need not be full column names but should be unique
        substrings.

    Returns
    -------
      The resulting data frame averaged over labels within a
      ROI cluster.
    '''
    columns = df.columns.values
    clusters = []
    for hemi, hcolumns in zip(['-lh', '-rh'], [lh(columns), rh(columns)]):
        for name, cols in all_clusters.items():
            cols = filter_cols(hcolumns, cols)
            if any([x in df.columns for x in ensure_iter(cols)]):
                cluster = df.loc[:, cols].mean(1)
                cluster.name = name + hemi
                clusters.append(cluster)
    clusters = pd.concat(clusters, 1)
    return clusters


def lateralize(data, ipsi, contra, suffix='_Lateralized'):
    '''
    Lateralize set of rois.

    Arguments
    ---------
      data : DataFrame
        Needs to contain ROIs as columns
      ipsi : list of str
        ROIs that are ipsilateral to event of interest
      contra : list of str
        ROIS that contralateral to event of interest

    Returns
    -------
      Input data frame with added columns that are contralateal -
      ipsilateral. These columns will have _Lateralized appended,
      but still have 'lh' in them.
    '''
    if not len(ipsi) == len(contra):
        raise RuntimeError('Ipsi and Contra lists must have same length')
    ipsi, contra = sorted(ipsi), sorted(contra)
    out = []
    for i, c in zip(ipsi, contra):
        out.append(data.loc[:, c] - data.loc[:, i])
        out[-1].name = i.replace('rh', 'lh').replace('PCeS',
                                                     'PCes') + suffix
    return pd.concat(out, 1)


'''
Following: Utility functions and plotting of ROIS on brain.
'''


def to_edges(vals):
    delta = np.diff(vals)[0]
    edges = vals - delta / 2.
    edges = np.array(list(edges) + [edges[-1] + delta])
    return edges


@memory.cache
def plot_roi(hemi, labels, colors, view='parietal',
             fs_dir=os.environ['SUBJECTS_DIR']):
    import os
    subject_id = "fsaverage"
    surf = "inflated"
    brain = Brain(subject_id, hemi, surf, offscreen=True)
    for label, color in zip(labels, colors):
        label_file = os.path.join(fs_dir, subject_id, 'label',
                                  (label.replace('-rh', '.label')
                                   .replace('-lh', '.label')
                                   .replace('&', '_and_')
                                   .replace('_Havg', '')
                                   .replace('_Lateralized', '')))
        brain.add_label(label_file, color=color)
    brain.show_view(view)
    return brain.screenshot()


def ensure_iter(input):
    if isinstance(input, str):
        yield input
    else:
        try:
            for item in input:
                yield item
        except TypeError:
            yield input

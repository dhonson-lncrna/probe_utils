from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import CompoundLocation

import pandas as pd
import numpy as np

def exon_intron(genbank, gene, _type_="ncRNA"):
    """Takes a genbank file and returns intron and exon sequences"""
    # Read in and extract the sequence
    record = SeqIO.read(genbank, "genbank")

    sequence = record.seq

    # Extract the desired feature
    feature = []
    
    for i in record.features:
        if i.type == _type_ and i.qualifiers["gene"] == gene:
            feature.append(i)

    if len(feature) < 1:
        raise ValueError("Feature not found")

    elif len(feature) > 1:
        raise ValueError("More than one feature of given description")

    else:
        location = feature[0].location

    # Use the exon indices to extract the mature transcript
    mat_rna = location.extract(sequence)
    
    # Make a list of exon indices
    exon_ls = list(location)
    
    # Make a list of intron indices using np.s_
    intron_slices = []
    
    for i, v in enumerate(exon_ls):

        if v != exon_ls[i - 1] + 1 and v > exon_ls[i - 1]:
            intron_slices.append(np.s_[exon_ls[i - 1] : v])

        else:
            pass

    # Extract the intron sequences
    introns = [sequence[i] for i in intron_slices]

    return mat_rna, introns

def rap_probes(seq, gene, probe_length = 90):
    '''Takes a sequence and makes probes of a given length'''
    # Extract indices of the desired probe length
    inds = np.arange(0, len(seq), probe_length)
    
    s_list = []
    
    for i in range(len(inds)-1):
        s_list.append(np.s_[inds[i]:inds[i+1]])
    
    # Use those indices to make probes
    s_seq = [seq[i] for i in s_list]
    
    # If there is more than a quarter probe of gene left uncovered, add one last probe 
    if len(seq) - inds[-1] > probe_length / 4 : 
        s_seq.append(seq[-90:])
    
    else:
        pass
    
    s_seq = [str(i.reverse_complement()) for i in s_seq]
    
    # Name the probes and return a dataframe
    prb_nms = [gene + str(i+1) for i in range(len(s_seq))]
    
    return pd.DataFrame({'Name':prb_nms,
                        'Sequence':s_seq})
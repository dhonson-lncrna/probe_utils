import pandas as pd
from Bio.Seq import Seq

neb_dict = {'AatII': 'GACGTC', 'Acc65I': 'GGTACC', 'AclI': 'AACGTT', 'AfeI': 'AGCGCT', 'AflII': 'CTTAAG', 'AgeI': 'ACCGGT', 'AluI': 'AGCT', 'ApaI': 'GGGCCC', 'ApaLI': 'GTGCAC', 'AscI': 'GGCGCGCC', 'AseI': 'ATTAAT', 'AsiSI': 'GCGATCGC', 'AvrII': 'CCTAGG', 'BamHI': 'GGATCC', 'BclI': 'TGATCA', 'BfaI': 'CTAG', 'BglII': 'AGATCT', 'BmtI': 'GCTAGC', 'BsiWI': 'CGTACG', 'BsmBI-v2': 'CGTCTC', 'BspEI': 'TCCGGA', 'BspHI': 'TCATGA', 'BsrGI': 'TGTACA', 'BssHII': 'GCGCGC', 'BstBI': 'TTCGAA', 'BstUI': 'CGCG', 'BstZ17I': 'GTATAC', 'ClaI\xa0BspDI': 'ATCGAT', 'CviAII': 'CATG', 'CviQI': 'GTAC', 'DpnI': 'GATC', 'DraI': 'TTTAAA', 'EagI': 'CGGCCG', 'Eco53kI': 'GAGCTC', 'EcoRI': 'GAATTC', 'EcoRV': 'GATATC', 'FatI': 'CATG', 'FseI': 'GGCCGGCC', 'FspI': 'TGCGCA', 'HaeIII': 'GGCC', 'HhaI': 'GCGC', 'HindIII': 'AAGCTT', 'HinP1I': 'GCGC', 'HpaI': 'GTTAAC', 'HpyCH4IV': 'ACGT', 'HpyCH4V': 'TGCA', 'KasI': 'GGCGCC', 'KpnI': 'GGTACC', 'MboI': 'GATC', 'MfeI': 'CAATTG', 'MluCI': 'AATT', 'MluI': 'ACGCGT', 'MscI': 'TGGCCA', 'MseI': 'TTAA', 'MspI HpaII': 'CCGG', 'NaeI': 'GCCGGC', 'NarI': 'GGCGCC', 'Nb.BbvCI': 'CCTCAGC', 'Nb.BsmI': 'GAATGC', 'Nb.BsrDI': 'GCAATG', 'Nb.BssSI': 'CACGAG', 'Nb.BtsI': 'GCAGTG', 'NcoI': 'CCATGG', 'NdeI': 'CATATG', 'NgoMIV': 'GCCGGC', 'NheI': 'GCTAGC', 'NlaIII': 'CATG', 'NotI': 'GCGGCCGC', 'NruI': 'TCGCGA', 'NsiI': 'ATGCAT', 'PacI': 'TTAATTAA', 'XhoI': 'CTCGAG', 'PciI': 'ACATGT', 'PluTI': 'GGCGCC', 'PmeI': 'GTTTAAAC', 'PmlI': 'CACGTG', 'PsiI-v2': 'TTATAA', 'PspOMI': 'GGGCCC', 'PstI': 'CTGCAG', 'PvuI': 'CGATCG', 'PvuII': 'CAGCTG', 'RsaI': 'GTAC', 'SacI': 'GAGCTC', 'SacII': 'CCGCGG', 'SalI': 'GTCGAC', 'SbfI': 'CCTGCAGG', 'ScaI': 'AGTACT', 'SfoI': 'GGCGCC', 'SmaI': 'CCCGGG', 'SnaBI': 'TACGTA', 'SpeI': 'ACTAGT', 'SphI': 'GCATGC', 'SrfI': 'GCCCGGGC', 'SspI': 'AATATT', 'StuI': 'AGGCCT', 'SwaI': 'ATTTAAAT', 'TaqI-v2': 'TCGA', 'XbaI': 'TCTAGA', 'XmaI': 'CCCGGG', 'ZraI': 'GACGTC'}

def rerap_adapt(re_name, 
                re_dict, 
                top_seq = 'GTGCCAGCCTTGGTCGAGAGCCGAGCCCGTGATGC',
                bot_seq = 'GCATCACGGGCTCGGCTCTCGACCAAGGCTGGCACaaaaaaaaaaaaaaaaaaaaaaaa',
                adapt = 'TGACTTG'):
    
    dig_site = Seq(re_dict[re_name])
    
    return (str(top_seq + dig_site), str(adapt + dig_site.reverse_complement() + bot_seq))

def rerap_bulk(re_ls,
              re_dict = neb_dict):
    site_ls = []
    name_ls = []
    
    for i in re_ls:
        adap_tup = rerap_adapt(i, re_dict=re_dict)
        site_ls.append(adap_tup[0])
        site_ls.append(adap_tup[1])
        name_ls.append(i + '_top')
        name_ls.append(i + '_bot')
        
    return pd.DataFrame({'Name':name_ls, 'Sequence':site_ls})
import pandas as pd
from Bio.Seq import Seq

def rerap_adapt(re_name, 
                re_dict, 
                top_seq = 'GTGCCAGCCTTGGTCGAGAGCCGAGCCCGTGATGC',
                bot_seq = 'GCATCACGGGCTCGGCTCTCGACCAAGGCTGGCACaaaaaaaaaaaaaaaaaaaaaaaa',
                adapt = 'TGACTTG'):
    
    dig_site = Seq(re_dict[re_name])
    
    return (str(top_seq + dig_site), str(adapt + dig_site.reverse_complement() + bot_seq))

def rerap_bulk(re_ls,
              re_dict):
    site_ls = []
    name_ls = []
    
    for i in re_ls:
        adap_tup = rerap_adapt(i, re_dict=re_dict)
        site_ls.append(adap_tup[0])
        site_ls.append(adap_tup[1])
        name_ls.append(i + '_top')
        name_ls.append(i + '_bot')
        
    return pd.DataFrame({'Name':name_ls, 'Sequence':site_ls})
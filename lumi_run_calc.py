#!/usr/bin/env python

"""
Calculate total lumi per Run for given triggers

Needed for scaling, and for run-dependant SF
"""


import pandas as pd


# run numbers for each run period, each run period includes the latter entry
RUN_DICT = {
    "B": [272007, 275376],
    "C": [275657, 276283],
    "D": [276315, 276811],
    "E": [276831, 277420],
    "F": [277772, 278808],
    "G": [278820, 280385],
    "H": [280919, 284044],
}


def calc_lumi_per_run(df):
    """Get dict of total record lumi per run period"""
    run_results = {}
    for k, (run_start, run_end) in RUN_DICT.items():
        run_results[k] = df[(df['run'] >= run_start) & (df['run'] <= run_end)]['recorded(/pb)'].sum()
    return run_results


def print_BtoF_GtoH(run_results):
    """Print lumi for B-F, G-H, and the fraction of total that's in G-H"""
    sum_BtoF = sum([v for k,v in run_results.items() if k in ['B', 'C', 'D', 'E', 'F']])
    sum_GtoH = sum([v for k,v in run_results.items() if k in ['G', 'H']])
    print("BtoF:", sum_BtoF)
    print("GtoH:", sum_GtoH)
    print("GtoH fraction:", sum_GtoH / (sum_BtoF + sum_GtoH))
    print("Total lumi:", sum_BtoF + sum_GtoH)


def read_brilcalc_csv(filename):
    """Read bricalc lumi CSV into pandas dataframe"""
    # awful way to figure out how many rubbish lines at bottom...can't we just
    # parse this thing once it's been read in?
    # FIXME use StringIO with read_csv, or use read_table
    with open(filename) as f:
        contents = f.readlines()
    len_footer = 0
    count_footer = False
    for line in contents:
        if "#Summary:" in line:
            count_footer = True
        if count_footer:
            len_footer += 1

    # #run:fill,time,ncms,hltpath,delivered(/pb),recorded(/pb)
    df = pd.read_csv(filename, 
                     usecols=[0, 3, 5], 
                     skiprows=1, 
                     skipfooter=len_footer, 
                     engine='python')
    # extract run number from #run:fill column, drop the original
    df['run'] = df['#run:fill'].str.split(":", expand=True).loc[:,0].astype('int32')
    df['hltpath'] = df['hltpath'].str.replace('_v.*', '', case=False)
    df.drop('#run:fill', axis=1, inplace=True)
    # print(df.head())
    # print(df.dtypes)
    return df


def do_results_per_trig(df):
    """Calculate results & printout per trigger path

    Have to do it this way, since a given run might have several triggers firing,
    so don't want to double count its lumi
    """
    for name, group in df.groupby('hltpath'):
        print("----", name, "----")
        results = calc_lumi_per_run(group)
        print(results)
        print_BtoF_GtoH(results)


if __name__ == "__main__":
    # brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json \
    #  -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt \
    #  -u /pb --hltpath 'HLT_ZeroBias_v*' -o zero_bias.csv
    df_zb = read_brilcalc_csv("zero_bias.csv")
    do_results_per_trig(df_zb)

    # brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json \
    #  -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt \
    #  -u /pb --hltpath 'HLT_PFJet*_v*' -o pfjet.csv
    df_jetht = read_brilcalc_csv("pfjet.csv")
    do_results_per_trig(df_jetht)

    # brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/Normtags/normtag_PHYSICS.json \
    #  -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt \
    #  -u /pb --hltpath 'HLT_Iso*Mu24_v*' -o single_mu.csv
    df_singlemu = read_brilcalc_csv("single_mu.csv")
    do_results_per_trig(df_singlemu)

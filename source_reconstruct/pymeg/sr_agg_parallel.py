
subjects = {'DCB': 3,
            'DHB': 3,
            'ECB': 3,
            'EMB': 3,
            'EXF': 3,
            'EXG': 3,
            'GSB': 3,
            'HBC': 3,
            'JTB': 3,
            'KSV': 3,
            'NIF': 3,
            'OMF': 3,
            'PDP': 3,
            'QNV': 4,
            'TFD': 3,
            'TNB': 3}
 #           'TSJ': 3}
            
#subjects = {'TSJ': 3}
            

def submit_aggregates(cluster='uke'):
    from pymeg import parallel
    for subject, session in subjects.items():
        for sessnum in range(1,session+1):
            for datatype in ['G']: #for datatype in ['F','BB','G']:
                parallel.pmap(aggregate, [(subject, sessnum, datatype)],
                      name='agg' + str(sessnum) + str(subject) + datatype,
                          tasks=3, memory=30, walltime='30:00:00', env="mne" # 3,30 for G, 5,50 for F,BB
                      )


def aggregate(subject, session, datatype):
    from pymeg import aggregate_sr as asr
    from os.path import join
    data = (
        '/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne_induced/%s-SESS%i-*%s*-lcmv.hdf' % (
            subject, session, datatype))

    if datatype == 'F':    # time-frequency
        agg = asr.aggregate_files(data, data, (-0.4, -0.2), to_decibels=True)
    elif datatype == 'BB':    # broadband
        agg = asr.aggregate_files(data, data, (-0.2, 0), to_decibels=False)
    elif datatype == 'G':    # gamma (optimized)
        agg = asr.aggregate_files(data, data, (-0.1, -0.05), to_decibels=True)

    filename = join(
        '/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne_induced/agg/',
        'S%s_SESS%i_%s_agg.hdf' % (subject, session, datatype))
    asr.agg2hdf(agg, filename)


if __name__=="__main__":
    submit_aggregates()

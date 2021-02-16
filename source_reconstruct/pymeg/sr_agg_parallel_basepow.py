
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
            'TNB': 3,
            'TSJ': 3}
            
#subjects = {'TSJ': 3}
            

def submit_aggregates(cluster='uke'):
    from pymeg import parallel
    for subject, session in subjects.items():
        for sessnum in range(1,session+1):
            for datatype in ['BB']: # only apply to BB
                parallel.pmap(aggregate, [(subject, sessnum, datatype)],
                      name='agg' + str(sessnum) + str(subject) + datatype,
                          tasks=3, memory=30, walltime='30:00:00', env="mne"
                      )


def aggregate(subject, session, datatype):
    from pymeg import aggregate_sr_basepow as asr
    from os.path import join
    data = (
        '/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/%s-SESS%i-*%s*-lcmv.hdf' % (
            subject, session, datatype))
    
    agg = asr.aggregate_files(data, data, (0, 5.5), hemis=['Averaged', 'Pair'])

    filename = join(
        '/home/pmurphy/Surprise_accumulation/Analysis/MEG/Conv2mne/agg/',
        'S%s_SESS%i_stimpow_agg.hdf' % (subject, session))
    asr.agg2hdf(agg, filename)


if __name__=="__main__":
    submit_aggregates()

#!/usr/bin/env python
'''
Parallelize functions with simple arguments via torque
'''

import errno
import os
import subprocess
import tempfile


def submit(walltime, memory, cwd, tmpdir,
           script, name, nodes='1:ppn=1',
           shellfname=None, env=None, ssh_to=None):
    '''
    Submit a script to torque
    '''

    cmd_top = '''
    #!/bin/bash
    # walltime: defines maximum lifetime of a job
    # nodes/ppn: how many nodes (usually 1)? how many cores?

    #PBS -q batch
    #PBS -l walltime={walltime}
    #PBS -l nodes={nodes}
    #PBS -l mem={memory}gb
    #PBS -N {name}

    cd {cwd}
    mkdir -p cluster
    chmod a+rwx cluster

    #### set journal & error options
    #PBS -o {cwd}/$PBS_JOBID.o
    #PBS -e {cwd}/$PBS_JOBID.e

    '''.format(**{'walltime': walltime,
                  'nodes': nodes,
                  'memory': memory,
                  'cwd': cwd,
                  'script': script,
                  'name': name})

    if env is not None:
        cmd_top += 'source activate %s\n' % env

    cmd_bottom = '''
    # FILE TO EXECUTE
    {script} 1> {cwd}/$PBS_JOBID.out 2> {cwd}/$PBS_JOBID.err
    '''.format(**{'walltime': walltime,
                  'nodes': nodes,
                  'memory': memory,
                  'cwd': cwd,
                  'script': script,
                  'name': name})
    command = cmd_top + cmd_bottom
    with tempfile.NamedTemporaryFile(delete=False, dir=tmpdir,
                                     prefix='delete_me_tmp') as shellfname:
        shellfname.write(command.encode('utf-8'))
        shellfname = shellfname.name
    if ssh_to is None:
        command = "qsub %s" % (shellfname)
    else:
        command = "ssh %s 'qsub %s'" % (ssh_to, shellfname)
    output = subprocess.check_output(
        command,
        stderr=subprocess.STDOUT,
        shell=True)
    return output


def slurm_submit(walltime, memory, tmpdir, logdir, script, name,
                 nodes=1, tasks=16, email=None, env=None,
                 shellfname=None):
    '''
    Submit a script to torque
    '''
    print('script in submit {}'.format(script))
    sbatch_directives = '''#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --nodes={nodes}
#SBATCH --tasks-per-node={tasks}
#SBATCH --time={walltime}
#SBATCH --export=NONE
#SBATCH --mem={memory}GB
#SBATCH --partition=std
    '''.format(walltime=walltime,
               nodes=nodes,
               memory=memory,
               tasks=tasks,
               name=name)
    if email is not None:
        sbatch_directives += '''
#SBATCH --mail-user={email}
#SBATCH --mail-type=ALL
        '''.format(email=email)
    sbatch_directives += '''
#SBATCH --error={logdir}/slurm_%j.out
#SBATCH --output={logdir}/slurm_%j.err
source /sw/modules/rrz-modules.sh
    '''.format(logdir=logdir)

    environment_variables = '''
module purge
module load env
module load site/hummel
source ~/.bashrc

{script}
    '''.format(script=script)
    command = sbatch_directives + environment_variables
    with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=tmpdir,
                                     prefix='sbatch_script') as shellfname:
        shellfname.write(command)
        shellfname = shellfname.name
    command = "sbatch %s" % (shellfname)
    output = subprocess.check_output(
        command,
        stderr=subprocess.STDOUT,
        shell=True)
    return output


def to_script(func, tmpdir, *args):
    '''
    Write a simple stub python function that calls this function.
    '''

    with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=tmpdir,
                                     prefix='py_submit_script') as script:
        code = """
print('Parameters:', '{function}', {qargs})
from {module} import {function}
{function}{args}
        """.format(module=func.__module__,
                   function=func.__name__,
                   args=args, qargs=str(args))
        script.write(code)
        return str(script.name)


def pmap(func, args, cluster='PBS', walltime=12, memory=10, logdir=None, tmpdir=None,
         name=None, nodes=1, tasks=1, verbose=True, env=None, email=None,
         ssh_to='node028', home=None):
    from os.path import expanduser, join
    if type(walltime) == type(int):
        walltime = '%i:00:00'
    if name is None:
        name = func.__name__
    if logdir is None:
        if home is None:
            home = expanduser("~")
        logdir = join(home, 'cluster_logs', func.__name__)
        mkdir_p(logdir)
    if tmpdir is None:
        if home is None:
            home = expanduser("~")
        tmpdir = join(home, 'cluster_logs', 'tmp')
        mkdir_p(tmpdir)
    out = []
    for arg in args:
        script = 'ipython ' + to_script(func, tmpdir, *arg)
        if verbose:
            print(arg, '->', script)
        if cluster.upper() == 'PBS':
            node_statement = '%i:ppn=%i' % (nodes, tasks)
            pid = submit(walltime, memory, logdir, tmpdir,
                         script, name, env=env, nodes=node_statement)
        elif cluster.upper() == 'SLURM':
            pid = slurm_submit(walltime, memory, logdir, tmpdir, script, name, env=env,
                               email=email, nodes=nodes, tasks=tasks)
        out.append(pid)
    return out


def status(pid):
    output = subprocess.check_output(
        "ssh node028 'qstat %s'" % pid.replace('\n', ''),
        stderr=subprocess.STDOUT,
        shell=True)
    if " C " in output.split('\n')[-2]:
        return True
    elif " E " in output.split('\n')[-2]:
        raise RuntimeError('Job %s failed')
    else:
        return False


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

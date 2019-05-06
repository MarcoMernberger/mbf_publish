#!/usr/bin/python3
import sys
import os
import pickle
import subprocess

###This is the register with website code"""
def load_meta_data():
    with open("web/scb/metaddata.dat") as op:
        return pickle.load(op)



def _write_scb_rsync_file():
    meta_data = load_meta_data()
    paths_to_copy = set()
    for group in meta_data:
        for entry in meta_data[group]:
            for x in ['table_path', 'path_bigbed', 'path_table', 'path_bam', 'path']:
                if x in entry:  # genes
                    paths_to_copy.add(entry[x])
            if 'path_bam' in entry:
                    paths_to_copy.add(entry['path_bam'] + '.bai')
    for fn in os.listdir('web/scb'):
        paths_to_copy.add(os.path.join('web','scb', fn))
    op = open("web/scb/rsync_list.txt", 'wb')
    for x in sorted(paths_to_copy):
        if os.path.abspath(x).startswith(os.path.abspath(os.path.join(os.path.dirname(__file__),'..', '..'))):
            op.write(x)
            op.write("\n")
    op.close()

def _rsync_to_server(project_name):
    cmd = [
    'sudo', 
    '-u',
     'ffs',
     'rsync',
     os.path.abspath('.') + '/',
     'ffs@mf:/mf/scb/%s/@chmod_after@chmod=o+rwX,g+rwX@chown=1000:2000' % project_name,
     '--files-from=web/scb/rsync_list.txt',
     '-e', "ssh -p 223",
     '--rsync-path=rprsync', '-r', 
     '-P',
     '-t', # otherwise we retrigger sync because of the chmod after transfer...
    ]
    print(cmd)
    subprocess.check_call(cmd)



def _register_with_server(accepted_server, path, revision):
    import requests
    auth = requests.auth.HTTPBasicAuth('feed', 'feed')
    if 'scb_server' in os.environ:
        top_level_url = os.environ['scb_server']
    else:
        top_level_url = accepted_server
    url = accepted_server + "/register/%s?revision=%s"  % (path, revision)
    req = requests.get(url, auth=auth)
    print('registered')



def get_vid_info(vid):
    import requests
    url = "http://imt.flofloflo.de/scb/vid_info/%s" % vid
    auth = requests.HTTPBasicAuth(
        os.environ('MBF_AUTH_USER'),
        os.environ('MBF_AUTH_PASSWORD'),
    )
    return requests.get(url).text

accepted_servers = {
        'scb': 'http://mbf.imt.uni-marburg.de/scb',
        'scb_dev': 'http://mbf.imt.uni-marburg.de/scb_dev',
        'localhost': 'http://localhost:8080/scb',
        }
def print_usage():
    print("Usage:")
    print('scb_submit.py')
    sys.exit()

def get_current_repo_revision():
    """Does not require auto commit"""
    x = subprocess.check_output(['hg','log','-r','tip','-q'])
    return x[:x.find(":")]

def main():
    path = os.path.abspath(os.getcwd())
    components = path.split("/")
    ok = components[-2] == 'e' and components[-3] == 'imt'
    ok |= components[-2] == 'e' and components[-3] == 'ffs'
    if not ok:
        raise ValueError("This must be run in a server/imt/e  or ffs/e folder")
    project_name = components[-1]
    print("submitting %s to scb..." % project_name)
    print('now rsyncing')
    _write_scb_rsync_file()
    _rsync_to_server(project_name)

    print('calling webserver')
    _register_with_server('http://mbf.imt.uni-marburg.de/scb', project_name, get_current_repo_revision())

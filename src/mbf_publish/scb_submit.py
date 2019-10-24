#!/usr/bin/python3
import sys
import os
import pickle
import subprocess
import json

###This is the register with website code"""
def load_meta_data():
    try:
        with open("web/scb/metadata.json") as op:
            return json.load(op)
    except OSError:
        raise ValueError(
            "web/scb/metadata.json not found - make sure your run script calls mbf_publish.scb.prep_scb()"
        )


def _rsync_to_server(project_name):
    cmd = [
        "sudo",
        #'-u',
        #'ffs',
        "rsync",
        os.path.abspath(".") + "/",
        "ffs@mbf.imt.uni-marburg.de:/mf/scb/%s/@chmod_after@chmod=o+rwX,g+rwX@chown=1000:2000"
        % (project_name),
        "--files-from=web/scb/rsync_list.txt",
        "-e",
        "ssh -p 223 -i /.ffs_ssh/id_rsa -o StrictHostKeyChecking=no",
        "--rsync-path=rprsync",
        "-r",
        "-P",
        "-t",  # otherwise we retrigger sync because of the chmod after transfer...
        "-v",
    ]
    print(cmd)
    subprocess.check_call(cmd)


def _register_with_server(accepted_server, path, revision):
    import requests

    auth = requests.auth.HTTPBasicAuth("feed", "feed")
    if "scb_server" in os.environ:
        top_level_url = os.environ["scb_server"]
    else:
        top_level_url = accepted_server
    url = accepted_server + "/register/%s?revision=%s" % (path, revision)
    print(url)
    req = requests.get(url, auth=auth)
    if req.status_code == 200:
        print("registered")
    else:
        print("error registring")
        print(req.text)


accepted_servers = {
    "scb": "http://mbf.imt.uni-marburg.de/scb",
    "scb_dev": "http://mbf.imt.uni-marburg.de/scb_dev",
    "localhost": "http://localhost:8080/scb",
}


def print_usage(msg=""):
    print("Usage:")
    print("scb_submit.py")
    if msg:
        print(msg)
    sys.exit()


def get_current_repo_revision():
    """Does not require auto commit"""
    x = subprocess.check_output(["hg", "log", "-r", "tip", "-q"]).decode("utf-8")
    return x[: x.find(":")]


def main():
    try:
        path = os.environ["ANYSNAKE_PROJECT_PATH"]
        project_name = path.split("/")[-1]
    except KeyError:
        print_usage("Must be run from inside anysnake")
    print("submitting %s to scb..." % project_name)
    print("now rsyncing")
    print("sudo password is test123 in anysnake container!")
    _rsync_to_server(project_name)

    print("calling webserver")
    _register_with_server(
        "http://mbf.imt.uni-marburg.de/scb", project_name, get_current_repo_revision()
    )

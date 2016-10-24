import os
import re
import sys
import time
import subprocess
from os.path import join
from os.path import dirname
from os.path import basename
from os.path import exists
from urlparse import urlparse
from itertools import chain
from collections import defaultdict
import xml.etree.ElementTree as ET

from cutlass.aspera import aspera as asp
from anadama.util import addtag

from . import ssh
from .serialize import indent
from .serialize import to_xml
from .util import reportnum
from .update import update_osdf_from_report


class NoEqual(object):
    def __eq__(self, other):
        return False

    
def fsize(fname):
    return os.stat(fname).st_size


def parse_fasp_url(u):
    parsed = urlparse(u)
    return parsed.netloc, parsed.path


identity = lambda x: x
def groupby(keyfunc=identity, seq=[]):
    grouped = defaultdict(list)
    for item in seq:
        grouped[keyfunc(item)].append(item)
    return grouped


def _sequences(sample_records):
    def _s():
        for sample in sample_records:
            for prep, seq in sample.prepseqs:
                yield seq

    grouped = groupby(lambda s: s.urls[0], _s())
    return [ grp[0] for grp in grouped.itervalues() ]

def _completeparse(fname):
    with open(fname) as f:
        ret = []
        for line in f:
            ret.append( line.split('\t') )
    return ret


def untar(fname):
    proc = subprocess.Popen(["tar", "-xvf", fname,], 
                            stdout=subprocess.PIPE)
    out = proc.communicate()[0]
    names = filter(bool, map(str.strip, out.split('\n')))
    return names, filter(os.path.isfile, names)


class DownUpUpToDate(object):
    def __init__(self, seq, ssh_session):
        self.seq = seq
        self.ssh_session = ssh_session

    def __call__(self, task, values):
        cf = task.targets[0]
        t = re.sub(r'\....\.complete$', '', cf)
        if not exists(t) or not exists(cf):
            return False
        if not os.stat(t).st_size == self.seq.size:
            return False
        with open(cf) as f:
            for name, size in [ s.split("\t") for s in map(str.strip, f) ]:
                key = os.path.join(self.ssh_session.remote_path, name)
                remote_size = self.ssh_session.file_cache.get(key, NoEqual())
                if not remote_size == int(size):
                    return False
        return True
    

def download_upload(recs_16s, cached_16s_files, recs_wgs, 
                    cached_wgs_files, dcc_user, dcc_pw, ncbi_srv, 
                    ncbi_path, ncbi_user, ncbi_keyfile, products_dir):

    ssh_session = ssh.SSHConnection(ncbi_user, ncbi_srv, ncbi_keyfile, ncbi_path)

    cached_dir_16s = dirname(cached_16s_files[0]) if cached_16s_files else products_dir
    cached_dir_wgs = dirname(cached_wgs_files[0]) if cached_wgs_files else products_dir
    six_files = set(map(basename, cached_16s_files))
    wgs_files = set(map(basename, cached_wgs_files))
    
    def _du(url, local_dir, local_cached, remote_size, namespace):
        def _actually_du():
            srv, remote_path = parse_fasp_url(url)
            bn = basename(remote_path)
            local_file = join(local_dir, bn)
            skip = (bn in local_cached 
                    and os.stat(local_file).st_size == remote_size)
            if skip == False:
                ret = asp.download_file(srv, dcc_user, dcc_pw, 
                                        remote_path, local_dir)
                if not ret:
                    raise Exception("Download failed: "+url)
            to_rm, files_to_upload = untar(local_file)
            for i, f in enumerate(files_to_upload):
                new_f = addtag(f, namespace)
                os.rename(f, new_f)
                files_to_upload[i] = new_f
            names_sizes = [(basename(f), os.stat(f).st_size) 
                           for f in files_to_upload]
            for f in files_to_upload:
                ret = asp.upload_file(ncbi_srv, ncbi_user, None, f,
                                      ncbi_path, keyfile=ncbi_keyfile)
            with open(local_file+"."+namespace+".complete", 'w') as f:
                for name_size in names_sizes:
                    print >> f, "\t".join(map(str, name_size))
            for f in reversed(to_rm):
                try:
                    os.rmdir(f) if os.path.isdir(f) else os.remove(f)
                except:
                    print >> sys.stderr, "Unable to remove "+f
        return _actually_du

    complete_16s, complete_wgs, tasks = [],[], []
                
    args = ([six_files, cached_dir_16s, recs_16s, complete_16s, "16s"],
            [wgs_files, cached_dir_wgs, recs_wgs, complete_wgs, "wgs"])
    for local_files, local_dir, recs, result_container, namespace in args:
        for seq in _sequences(recs):
            remote_fname = basename(seq.urls[0])
            target = join(local_dir, remote_fname)
            
            if not seq.urls:
                raise Exception("Sequence ID %s has no urls"%(seq.id))
            tasks.append(
                { "name": "serialize:download_upload: "+remote_fname+"."+namespace,
                  "actions": [_du(seq.urls[0], local_dir, 
                                  local_files, seq.size, namespace)],
                  "file_dep": [],
                  "uptodate": [DownUpUpToDate(seq, ssh_session)],
                  "targets": [target+"."+namespace+".complete"] }
                )
            result_container.append(target+"."+namespace+".complete")
    return complete_16s, complete_wgs, tasks
        

def serialize(session, study, records_16s, files_16s, records_wgs, files_wgs,
              unsequenced_records, submission_fname, ready_fname, products_dir, 
              dcc_user, dcc_pw, study_id=None, release_date=None, 
              bioproject_id=None):
    """
    Download raw sequence files and serialize metadata into xml for a
    cutlass.Study

    :param dcc_user: the user used for the cutlass.iHMPSession

    :param dcc_pw: String; the password used for the cutlass.iHMPSession

    :param study_id: String; OSDF-given ID for the study you want to serialize
    """


    def _write_xml():
        tardict = {}
        for complete_fname in chain(files_16s, files_wgs):
            seqtype = re.sub(r'.*\.(...)\.complete$', r'\1', complete_fname)
            key = (basename(re.sub(r'\....\.complete$', '', complete_fname)), seqtype)
            tardict[key] = _completeparse(complete_fname)
        samples = list(records_16s)+list(records_wgs)+list(unsequenced_records)
        xml = to_xml(study, samples, tardict, release_date, bioproject_id)
        indent(xml)
        et = ET.ElementTree(xml)
        et.write(submission_fname)

    yield {
        "name": "serialize:xml: "+submission_fname,
        "actions": [_write_xml],
        "file_dep": list(files_16s)+list(files_wgs),
        "targets": [submission_fname]
    }

    yield {
        "name": "serialize:ready_file: "+ready_fname,
        "actions": [lambda *a, **kw: open(ready_fname, 'w').close()],
        "file_dep": [],
        "targets": [ready_fname]
    }


def kickoff(sub_fname, ready_fname, complete_fnames, keyfile,
            remote_path, remote_srv, user, products_dir):
    """Upload raw sequence files and xml.

    :param keyfile: String; absolute filepath to private SSH keyfile for
    access to NCBI's submission server

    :param remote_path: String; the directory on the NCBI submission
    server where to upload data. If unset, the remote_path is
    automatically determined.

    :param remote_srv: String; TLD of NCBI's submission server

    :param user: String; username used to access NCBI's submission server

    """

    def _upload(local_fname, complete_fname, blithely=False):
        def _u():
            ret = asp.upload_file(remote_srv, user, None, local_fname,
                                  remote_path, keyfile=keyfile)
            if blithely or ret:
                open(complete_fname, 'w').close()
            return blithely or ret # return True if blithely is True
        return _u

    yield {
        "name": "upload: "+basename(sub_fname),
        "actions": [_upload(sub_fname, sub_fname+".complete")],
        "file_dep": complete_fnames,
        "targets": [sub_fname+".complete"]
    }

    yield {
        "name": "upload: "+basename(ready_fname),
        "actions": [_upload(ready_fname, ready_fname+".complete", True)],
        "file_dep": complete_fnames+[sub_fname+".complete"],
        "targets": [ready_fname+".complete"]
    }


def report(session, ready_complete_fname, user, remote_srv,
           remote_path, keyfile):
    reports_dir = dirname(ready_complete_fname)
    def _download():
        c = ssh.SSHConnection(user, remote_srv, keyfile, remote_path)
        for _ in range(60*20*2): # 20 minutes in half-seconds
            report_fnames = [basename(n) for n in c.files()
                             if re.search(r'report\.[\d.]*xml', n)
                             and not exists(join(reports_dir, basename(n)))]
            if report_fnames:
                break
            else:
                time.sleep(.5)
        for n in report_fnames:
            if exists(join(reports_dir, n)):
                continue
            asp.download_file(remote_srv, user, None, join(remote_path, n),
                              reports_dir, keyfile=keyfile)
        if not report_fnames:
            print >> sys.stderr, "Timed out waiting for report xml files."
            return False
        most_recent_report = max(report_fnames, key=reportnum)
        update_osdf_from_report(session, join(reports_dir, most_recent_report))

    yield {
        "name": "report:get_reports",
        "actions": [_download],
        "file_dep": [ready_complete_fname],
        "uptodate": [False],
        "targets": [],
    }
    



import os
import re
import sys
import time
from os.path import join
from os.path import dirname
from os.path import basename
from os.path import exists
from urlparse import urlparse
from collections import defaultdict
import xml.etree.ElementTree as ET

from anadama.util import new_file
import cutlass.aspera as asp

from . import ssh
from .serialize import indent
from .serialize import to_xml
from .util import reportnum
from .update import update_osdf_from_report

def fsize(fname):
    return os.stat(fname).st_size

def parse_fasp_url(u):
    parsed = urlparse(u)
    return parsed.netloc, parsed. path

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


def serialize(session, study, records_16s, files_16s, records_wgs, files_wgs,
              unsequenced_records, submission_fname, ready_fname, products_dir, 
              dcc_user, dcc_pw, study_id=None, release_date=None):
    """
    Download raw sequence files and serialize metadata into xml for a
    cutlass.Study

    :param dcc_user: the user used for the cutlass.iHMPSession

    :param dcc_pw: String; the password used for the cutlass.iHMPSession

    :param study_id: String; OSDF-given ID for the study you want to serialize
    """


    cached_dir_16s = dirname(files_16s[0]) if files_16s else products_dir
    cached_dir_wgs = dirname(files_wgs[0]) if files_wgs else products_dir
    local_fnames_16s = set()
    local_fnames_wgs = set()
    for fname in files_16s:
        local_fnames_16s.add(basename(fname))
    for fname in files_wgs:
        local_fnames_wgs.add(basename(fname))

    def _download(url, local_dir):
        def _d():
            srv, remote_path = parse_fasp_url(url)
            ret = asp.download_file(srv, dcc_user, dcc_pw,
                                    remote_path, local_dir)
            return ret
        return _d

    args = ([local_fnames_16s, cached_dir_16s, records_16s, files_16s],
            [local_fnames_wgs, cached_dir_wgs, records_wgs, files_wgs])
    for skip_these, local_dir, rec_container, file_container in args:
        for seq in _sequences(rec_container):
            remote_fname = basename(seq.urls[0])
            target = join(local_dir, remote_fname)
            u2d = lambda *a, **kw: exists(target) and fsize(target) == seq.size
            if not seq.urls:
                raise Exception("Sequence ID %s has no urls"%(seq.id))
            if remote_fname in skip_these and fsize(target) == seq.size:
                continue
            yield {
                "name": "serialize:download: "+remote_fname,
                "actions": [_download(seq.urls[0], local_dir)],
                "file_dep": [],
                "uptodate": [u2d],
                "targets": [target]
            }
            file_container.append(target)

    def _write_xml():
        samples = list(records_16s)+list(records_wgs)+list(unsequenced_records)
        xml = to_xml(study, samples, release_date)
        indent(xml)
        et = ET.ElementTree(xml)
        et.write(submission_fname)

    yield {
        "name": "serialize:xml: "+submission_fname,
        "actions": [_write_xml],
        "file_dep": [],
        "targets": [submission_fname]
    }

    yield {
        "name": "serialize:ready_file: "+ready_fname,
        "actions": [lambda *a, **kw: open(ready_fname, 'w').close()],
        "file_dep": [],
        "targets": [ready_fname]
    }


def upload(files_16s, files_wgs, sub_fname, ready_fname, keyfile,
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

    to_upload = [ f for f in list(files_16s)+list(files_wgs)
                  if not f.endswith(".complete") ]
    ssh_session = ssh.SSHConnection(user, remote_srv, keyfile, remote_path)
    uptodate = [ssh_session.uptodate]

    def _upload(local_fname, complete_fname, blithely=False):
        def _u():
            ret = asp.upload_file(remote_srv, user, None, local_fname,
                                  remote_path, keyfile=keyfile)
            if blithely or ret:
                open(complete_fname, 'w').close()
            return blithely or ret # return True if blithely is True
        return _u

    complete_fnames = [new_file(f+".complete", basedir=products_dir) 
                       for f in to_upload]
    for f, complete_fname in zip(to_upload, complete_fnames):
        yield {
            "name": "upload: "+basename(f),
            "actions": [_upload(f, complete_fname)],
            "file_dep": [f],
            "uptodate": uptodate,
            "targets": [complete_fname]
        }

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
    



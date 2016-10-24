import os
import getpass

import cutlass
import anadama.pipelines

from . import workflows
from . import SubmitRecord
from . import PrepSeq

def get_prepseqs(preps):
    def _ps():
        for prep in preps:
            meth = getattr(prep, "child_seq_sets", None)
            if not meth:
                meth = getattr(prep, "raw_seq_sets", None)
            for seq in meth():
                yield prep, seq
    return map(PrepSeq._make, _ps())

def filter_unsequenced(records_wgs, records_16s):
    recs_16s, recs_wgs, unsequenced = list(), list(), set()
    for rec in records_16s:
        if not rec.prepseqs:
            unsequenced.add(rec.sample)
        else:
            recs_16s.append(rec)
    for rec in records_wgs:
        if not rec.prepseqs:
            unsequenced.add(rec.sample)
        else:
            recs_wgs.append(rec)
    return [SubmitRecord(s, []) for s in unsequenced], recs_16s, recs_wgs
        

def _remote_path(options):
    study_id = options['serialize']['study_id']
    return "/submit/Production/{}/".format(study_id)


class DCCSRAPipeline(anadama.pipelines.Pipeline):
    """Pipeline for submitting metadata from the iHMP DCC's OSDF instance
    to NCBI's SRA.

    Steps:

    1. Query OSDF for all samples, preps, and raw sequence sets for a
       given study

    2. For each raw sequence, download the raw sequence file if it's
    not available locally

    3. Serialize all metadata useful for SRA from OSDF into a
    submission.xml file

    4. Create an empty submit.empty file.

    5. Upload raw sequence files to SRA as necessary

    6. Upload submission.xml and submit.ready file

    Workflows used:

    * :py:func:`dcc_sra.workflows.download_upload`
    * :py:func:`dcc_sra.workflows.serialize`
    * :py:func:`dcc_sra.workflows.kickoff`
    * :py:func:`dcc_sra.workflows.report`

    """


    name = "DCCSRA"
    products = {
        "cached_16s_files": list(),
        "cached_wgs_files": list()
    }

    default_options = {
        "serialize": {
            "dcc_user": None,
            "dcc_pw": None,
            "study_id": None,
            "release_date": None,
            "bioproject_id": None,
        },
        "upload": {
            "keyfile": "/home/rschwager/test_data/broad_metadata/dcc_sra/iHMP_SRA_key",
            "remote_path": None,
            "remote_srv" : "upload.ncbi.nlm.nih.gov",
            "user": "asp-hmp2",
        },
        "report": {
            "products_dir": "reports"
        }
    }

    workflows = {
        "serialize": workflows.serialize,
        "kickoff": workflows.kickoff,
        "download_upload": workflows.download_upload,
    }

    def __init__(self, cached_16s_files=list(),
                 cached_wgs_files=list(),
                 products_dir=str(),
                 workflow_options=dict(),
                 *args, **kwargs):

        """Initialize the pipeline.

        :keyword cached_16s_files: List of strings; raw 16S sequence
        files (fasta or fastq format) already downloaded.

        :keyword cached_16s_files: List of strings;raw WGS sequence
        files (fasta or fastq format) already downloaded.

        :keyword products_dir: String; Directory path for where outputs will 
                               be saved.

        :keyword workflow_options: Dictionary; **opts to be fed into the 
                                   respective workflow functions.

        """

        super(DCCSRAPipeline, self).__init__(*args, **kwargs)

        self.options = self.default_options.copy()
        for k in self.options.iterkeys():
            self.options[k].update(workflow_options.get(k,{}))

        if not products_dir:
            products_dir = self.options['report']['products_dir']
        self.products_dir = os.path.abspath(products_dir)
        if not os.path.isdir(self.products_dir):
            os.mkdir(self.products_dir)

        if not self.options['serialize'].get('dcc_user', None):
            default = getpass.getuser()
            prompt = "Enter your DCC username: (%s)"%(default)
            entered = raw_input(prompt)
            if not entered:
                entered = default
            self.options['serialize']['dcc_user'] = entered

        if not self.options['serialize'].get('study_id', None):
            prompt = "Enter the study ID to submit: "
            self.options['serialize']['study_id'] = raw_input(prompt)

        if not self.options['serialize'].get('dcc_pw', None):
            prompt = "Enter your DCC password: "
            self.options['serialize']['dcc_pw'] = getpass.getpass(prompt)

        if not self.options['upload'].get('remote_path', None):
            self.options['upload']['remote_path'] = _remote_path(self.options)

        if not self.options['upload']['remote_path'].endswith('/'):
            self.options['upload']['remote_path'] += '/'

        self.add_products(
            cached_wgs_files = cached_wgs_files,
            cached_16s_files = cached_16s_files
        )


    def _configure(self):
        session = cutlass.iHMPSession(self.options['serialize']['dcc_user'],
                                      self.options['serialize']['dcc_pw'])
        study = cutlass.Study.load(self.options['serialize']['study_id'])

        records_wgs = list()
        records_16s = list()
        for subject in study.subjects():
            for visit in subject.visits():
                for sample in visit.samples():
                    prepseqs_16s = get_prepseqs(sample.sixteenSDnaPreps())
                    records_16s.append(SubmitRecord(sample, prepseqs_16s))
                    prepseqs_wgs = get_prepseqs(sample.wgsDnaPreps())
                    records_wgs.append(SubmitRecord(sample, prepseqs_wgs))
                    
        unsequenced, recs_16s, recs_wgs = filter_unsequenced(records_wgs,
                                                             records_16s)

        submission_file = os.path.join(self.products_dir, "submission.xml")
        ready_file = os.path.join(self.products_dir, "submit.ready")
        six_fnames, wgs_fnames, tasks = workflows.download_upload(
            recs_16s, self.cached_16s_files, 
            recs_wgs, self.cached_wgs_files, 
            dcc_user = self.options['serialize']['dcc_user'],
            dcc_pw = self.options['serialize']['dcc_pw'],
            ncbi_srv = self.options['upload']['remote_srv'],
            ncbi_path = self.options['upload']['remote_path'],
            ncbi_user = self.options['upload']['user'],
            ncbi_keyfile = self.options['upload']['keyfile'],
            products_dir = self.products_dir
            )
        for t in tasks:
            yield t

        yield workflows.serialize(session, study, recs_16s,
                                  six_fnames,
                                  recs_wgs,
                                  wgs_fnames,
                                  unsequenced,
                                  submission_file,
                                  ready_file,
                                  self.products_dir,
                                  **self.options['serialize'])

        yield workflows.kickoff(submission_file, ready_file,
                                six_fnames+wgs_fnames,
                                products_dir=self.products_dir,
                                **self.options['upload'])

        yield workflows.report(session, ready_file+".complete",
                               **self.options['upload'])

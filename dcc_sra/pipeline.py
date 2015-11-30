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
            for seq in prep.raw_seq_sets():
                yield prep, seq
    return map(PrepSeq._make, _ps())

def filter_unsequenced(records_wgs, records_16s):
    recs_16s, recs_wgs, unsequenced = list(), list(), set()
    for rec in records_wgs:
        if not rec.prepseqs:
            unsequenced.add(rec.sample)
        else:
            recs_16s.append(rec)
    for rec in records_16s:
        if not rec.prepseqs:
            unsequenced.add(rec.sample)
        else:
            recs_wgs.append(rec)
    return [SubmitRecord(s, []) for s in unsequenced], recs_16s, recs_wgs
        

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

    * :py:func:`dcc_sra.workflows.serialize`
    * :py:func:`dcc_sra.workflows.upload`

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
        },
        "upload": {
            "keyfile": "/home/rschwager/test_data/broad_metadata/dcc_sra/iHMP_SRA_key",
            "remote_path": "/submit/Test/Fulltest",
            "remote_srv" : "upload.ncbi.nlm.nih.gov",
            "user": "asp-hmp2",
            "ftp_pass": None
        },
        "report": {
            "products_dir": "reports"
        }
    }

    workflows = {
        "serialize": workflows.serialize,
        "upload": workflows.upload,
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
        self.options.update(workflow_options)

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
        yield workflows.serialize(session, study, records_16s,
                                  self.cached_16s_files,
                                  records_wgs,
                                  self.cached_wgs_files,
                                  unsequenced,
                                  submission_file,
                                  ready_file,
                                  self.products_dir,
                                  **self.options['serialize'])

        yield workflows.upload(self.cached_16s_files,
                               self.cached_wgs_files,
                               submission_file, ready_file,
                               **self.options['upload'])

        # yield workflows.report(self.products_dir, **self.options['upload'])

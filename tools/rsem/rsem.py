"""
RSEM datatypes
"""
import logging
import os
import os.path

from galaxy.datatypes.images import Html
from galaxy.datatypes.metadata import MetadataElement
from galaxy.datatypes.sniff import get_headers
from galaxy.datatypes.tabular import Tabular


log = logging.getLogger(__name__)


class RsemIsoformsResults(Tabular):
    file_ext = "rsem.isoforms.results"
    """
    required columns:
    transcript_id gene_id length effective_length expected_count TPM
    FPKM IsoPct
    optional columns:
    pme_expected_count pme_TPM pme_FPKM IsoPct_from_pme_TPM TPM_ci_lower_bound
    TPM_ci_upper_bound FPKM_ci_lower_bound FPKM_ci_upper_bound
    """

    def __init__(self, **kwd):
        Tabular.__init__(self, **kwd)
        """Initialize RsemResults datatype"""
        self.comment_lines = 1

    def sniff(self, filename):
        headers = get_headers(filename, '\n', count=1)
        return (len(headers) > 0 and len(headers[0]) >= 8 and
                headers[0][0] == "transcript_id" and
                headers[0][1] == "gene_id" and headers[0][6] == "FPKM")

    def set_meta(self, dataset, **kwd):
        Tabular.set_meta(self, dataset, skip=None, **kwd)


class RsemGenesResults(Tabular):
    file_ext = "rsem.genes.results"
    """
    required columns:
    gene_id transcript_id(s) length effective_length expected_count TPM FPKM
    optional columns:
    pme_expected_count pme_TPM pme_FPKM TPM_ci_lower_bound TPM_ci_upper_bound
    FPKM_ci_lower_bound FPKM_ci_upper_bound
    """

    def __init__(self, **kwd):
        Tabular.__init__(self, **kwd)
        """Initialize RsemResults datatype"""
        self.comment_lines = 1

    def sniff(self, filename):
        headers = get_headers(filename, '\n', count=1)
        return (len(headers) > 0 and len(headers[0]) >= 7 and
                headers[0][0] == "gene_id" and
                headers[0][1].startswith("transcript_id") and
                headers[0][6] == "FPKM")

    def set_meta(self, dataset, **kwd):
        Tabular.set_meta(self, dataset, skip=None, **kwd)


class RsemReference(Html):
    """Class describing an RSEM reference"""
    MetadataElement(name='reference_name', default='rsem_ref',
                    desc='RSEM Reference Name', readonly=True, visible=True,
                    set_in_upload=True, no_value='rsem_ref')
    file_ext = 'rsem_ref'
    allow_datatype_change = False
    composite_type = 'auto_primary_file'

    def __init__(self, **kwd):
        Html.__init__(self, **kwd)
        """
        Expecting files:
        extra_files_path/<reference_name>.grp
        extra_files_path/<reference_name>.ti
        extra_files_path/<reference_name>.seq
        extra_files_path/<reference_name>.transcripts.fa
        Optionally includes files:
        extra_files_path/<reference_name>.chrlist
        extra_files_path/<reference_name>.idx.fa
        extra_files_path/<reference_name>.1.ebwt
        extra_files_path/<reference_name>.2.ebwt
        extra_files_path/<reference_name>.3.ebwt
        extra_files_path/<reference_name>.4.ebwt
        extra_files_path/<reference_name>.rev.1.ebwt
        extra_files_path/<reference_name>.rev.2.ebwt
        """
        self.add_composite_file('%s.grp', description='Group File',
                                substitute_name_with_metadata='reference_name',
                                is_binary=False)
        self.add_composite_file('%s.ti', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=False)
        self.add_composite_file('%s.seq', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=False)
        self.add_composite_file('%s.transcripts.fa', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=False)
        self.add_composite_file('%s.chrlist', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=False, optional=True)
        self.add_composite_file('%s.idx.fa', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=False, optional=True)
        self.add_composite_file('%s.1.ebwt', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=True, optional=True)
        self.add_composite_file('%s.2.ebwt', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=True, optional=True)
        self.add_composite_file('%s.3.ebwt', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=True, optional=True)
        self.add_composite_file('%s.4.ebwt', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=True, optional=True)
        self.add_composite_file('%s.rev.1.ebwt', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=True, optional=True)
        self.add_composite_file('%s.rev.2.ebwt', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=True, optional=True)

    def generate_primary_file(self, dataset=None):
        """
        This is called only at upload to write the file
        cannot rename the datasets here - they come with
        the default unfortunately
        """

    def regenerate_primary_file(self, dataset):
        """
        cannot do this until we are setting metadata
        """
        link_to_exts = ['.grp', '.ti', '.seq', '.fa', '.chrlist', '.log']
        ref_name = dataset.metadata.reference_name
        efp = dataset.extra_files_path
        flist = os.listdir(efp)
        rval = ['<html><head><title>%s</title></head><body><p/>RSEM \
                 Reference   %s   files:<p/><ul>' % (dataset.name, ref_name)]
        rvalb = []
        for i, fname in enumerate(flist):
            sfname = os.path.split(fname)[-1]
            f, e = os.path.splitext(fname)
            if e in link_to_exts:
                rval.append('<li><a href="%s">%s</a></li>' % (sfname, sfname))
            else:
                rvalb.append('<li>%s</li>' % (sfname))
        if len(rvalb) > 0:
            rval += rvalb
        rval.append('</ul></body></html>')
        fh = open(dataset.file_name, 'w')
        fh.write("\n".join(rval))
        fh.write('\n')
        fh.close()

    def set_meta(self, dataset, **kwd):
        Html.set_meta(self, dataset, **kwd)
        efp = dataset.extra_files_path
        flist = os.listdir(efp)
        for i, fname in enumerate(flist):
            if fname.endswith('.grp'):
                dataset.metadata.reference_name = fname[:-4]
                break
        self.regenerate_primary_file(dataset)


class RsemBt2Reference(Html):
    """Class describing an RSEM reference"""
    MetadataElement(name='reference_name', default='rsem_bt2_ref',
                    desc='RSEM Bowtie2 Reference Name', readonly=True,
                    visible=True, set_in_upload=True, no_value='rsem_bt2_ref')
    file_ext = 'rsem_bt2_ref'
    allow_datatype_change = False
    composite_type = 'auto_primary_file'

    def __init__(self, **kwd):
        Html.__init__(self, **kwd)
        """
        Expecting files:
        extra_files_path/<reference_name>.grp
        extra_files_path/<reference_name>.ti
        extra_files_path/<reference_name>.seq
        extra_files_path/<reference_name>.transcripts.fa
        Optionally includes files:
        extra_files_path/<reference_name>.chrlist
        extra_files_path/<reference_name>.idx.fa
        extra_files_path/<reference_name>.n2g.idx.fa
        extra_files_path/<reference_name>.1.bt2
        extra_files_path/<reference_name>.2.bt2
        extra_files_path/<reference_name>.3.bt2
        extra_files_path/<reference_name>.4.bt2
        extra_files_path/<reference_name>.rev.1.bt2
        extra_files_path/<reference_name>.rev.2.bt2
        """
        self.add_composite_file('%s.grp', description='Group File',
                                substitute_name_with_metadata='reference_name',
                                is_binary=False)
        self.add_composite_file('%s.ti', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=False)
        self.add_composite_file('%s.seq', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=False)
        self.add_composite_file('%s.transcripts.fa', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=False)
        self.add_composite_file('%s.chrlist', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=False, optional=True)
        self.add_composite_file('%s.idx.fa', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=False, optional=True)
        self.add_composite_file('%s.n2g.idx.fa', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=False, optional=True)
        self.add_composite_file('%s.1.bt2', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=True, optional=True)
        self.add_composite_file('%s.2.bt2', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=True, optional=True)
        self.add_composite_file('%s.3.bt2', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=True, optional=True)
        self.add_composite_file('%s.4.bt2', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=True, optional=True)
        self.add_composite_file('%s.rev.1.bt2', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=True, optional=True)
        self.add_composite_file('%s.rev.2.bt2', description='',
                                substitute_name_with_metadata='reference_name',
                                is_binary=True, optional=True)

    def generate_primary_file(self, dataset=None):
        """
        This is called only at upload to write the file
        cannot rename the datasets here - they come with
        the default unfortunately
        """

    def regenerate_primary_file(self, dataset):
        """
        cannot do this until we are setting metadata
        """
        link_to_exts = ['.grp', '.ti', '.seq', '.fa', '.chrlist', '.log']
        ref_name = dataset.metadata.reference_name
        efp = dataset.extra_files_path
        flist = os.listdir(efp)
        rval = ['<html><head><title>%s</title></head><body><p/>RSEM \
                 Reference   %s   files:<p/><ul>' % (dataset.name, ref_name)]
        rvalb = []
        for i, fname in enumerate(flist):
            sfname = os.path.split(fname)[-1]
            f, e = os.path.splitext(fname)
            if e in link_to_exts:
                rval.append('<li><a href="%s">%s</a></li>' % (sfname, sfname))
            else:
                rvalb.append('<li>%s</li>' % (sfname))
        if len(rvalb) > 0:
            rval += rvalb
        rval.append('</ul></body></html>')
        fh = open(dataset.file_name, 'w')
        fh.write("\n".join(rval))
        fh.write('\n')
        fh.close()

    def set_meta(self, dataset, **kwd):
        Html.set_meta(self, dataset, **kwd)
        efp = dataset.extra_files_path
        flist = os.listdir(efp)
        for i, fname in enumerate(flist):
            if fname.endswith('.grp'):
                dataset.metadata.reference_name = fname[:-4]
                break
        self.regenerate_primary_file(dataset)

"""
Hmmer classes
"""

import data
import logging
import re
import string
from cgi import escape
from galaxy.datatypes.metadata import MetadataElement
from galaxy.datatypes import metadata
import galaxy.model
from galaxy import util
from sniff import *

log = logging.getLogger(__name__)

class Hmm( data.Text ):
    """Class for hmmer database files"""

    file_ext = 'hmm'

    def init_meta( self, dataset, copy_from=None ):
        data.Text.init_meta( self, dataset, copy_from=copy_from )

class HmmPressed( Hmm ):
    """Class describing a hmmer database produced by hmmpress"""

    file_ext = 'hmmPressed'
    composite_type='basic'

    MetadataElement( readonly=True, optional=True, visible=False, no_value=0 )

    def __init__(self,**kwd):
        data.Data.__init__(self, **kwd)
        self.add_composite_file('hmm.h3m')
        self.add_composite_file('hmm.h3i')
        self.add_composite_file('hmm.h3f')
        self.add_composite_file('hmm.h3p')
    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek  = "Folder of multiple files"
            dataset.blurb = "Folder of multiple files"
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'
    def display_peek( self, dataset ):
        try:
            return dataset.peek
        except:
            return "Folder of multiple files"
    def get_mime(self):
        """Returns the mime type of the datatype"""
        return 'text/plain'

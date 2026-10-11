""" FeatureXML Filter """
from typing import Union
from pyopenms import FeatureMap, ConsensusMap, PeptideIdentification

class FeatureXMLFilter():
    """ FeatureXML Filter """
    def __init__(self, ):
        """ Constructor """
        self.fmap:FeatureMap = FeatureMap()
        self.fmap_filtered:FeatureMap = None
        self.blacklist:Union[ConsensusMap, PeptideIdentification] = None

    def init_fmap_filtered(self):
        """ Initialize the filtered fmap. """
        self.fmap_filtered = FeatureMap(self.fmap)
        self.fmap_filtered.clear(False)
        self.fmap_filtered.setProteinIdentifications([])
        self.fmap_filtered.setUnassignedPeptideIdentifications([])

    def filter_feature_with_internal_id(self):
        """ Filter features with peptide identifications. """
        self.init_fmap_filtered()

        for feature in self.fmap:
            if not feature.getPeptideIdentifications():
                self.fmap_filtered.push_back(feature)

        self.fmap_filtered.ensureUniqueId()

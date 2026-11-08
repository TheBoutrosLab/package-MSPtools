""" Unit test cases for ProteinQuantCSV """
from unittest import TestCase
from io import StringIO
import numpy as np
from pyopenms import (
    AASequence, BaseFeature, ConsensusFeature, ConsensusMap, PeptideHit,
    PeptideIdentification
)
from msptools.ProteinQuantCSV.ProteinQuantCSV import ProteinQuantCSV, PeptideQuantRecord
from msptools.ProteinQuantCSV.ProteinQuantCSVMetadata import ProteinQuantCSVMetadata

class TestProteinQuantCSV(TestCase):
    """ Test cases for ProteinQuantCSV """
    @staticmethod
    def peptide_charge_map():
        """ A feature and its identification can have different charges. """
        cmap = ConsensusMap()
        cf = ConsensusFeature()
        cf.setCharge(3)
        hit = PeptideHit()
        hit.setSequence(AASequence.fromString('PEPTIDEK'))
        hit.setCharge(2)
        identification = PeptideIdentification()
        identification.setHits([hit])
        cf.setPeptideIdentifications([identification])
        feature = BaseFeature()
        feature.setIntensity(10)
        cf.insert(0, feature)
        cmap.push_back(cf)
        return cmap

    @staticmethod
    def peptide_quant_record(charge, abundance):
        """ Construct one peptide row with two channels in the same plex. """
        return PeptideQuantRecord(
            abundance=abundance, header=['ch0', 'ch1'], peptide='PEPTIDEK',
            protein='P1', n_proteins=1, charge=charge, fraction='plex'
        )

    def test_remove_unquantified_peptides_uses_identified_hit_charge(self):
        """ Keep quantified abundance when only the peptide-hit charge matches. """
        record = self.peptide_quant_record(2, [10, 0])
        quant = ProteinQuantCSV(
            records=[record],
            metadata=ProteinQuantCSVMetadata(
                header_map={1: ('plex', 'ch0'), 2: ('plex', 'ch1')}
            )
        )
        summary = []
        quant.remove_unquantified_peptides(self.peptide_charge_map(), summary)
        self.assertEqual(record.abundance[0], 10)
        self.assertTrue(np.isnan(record.abundance[1]))
        self.assertEqual(summary, [('ch0', 0), ('ch1', 1)])

    def test_remove_unquantified_peptides_reports_missing_mapping(self):
        for peptide, charge in [('PEPTIDEK', 4), ('UNKNOWN', 2)]:
            record = self.peptide_quant_record(charge, [10, 0])
            record.peptide = peptide
            quant = ProteinQuantCSV(
                records=[record],
                metadata=ProteinQuantCSVMetadata(
                    header_map={1: ('plex', 'ch0'), 2: ('plex', 'ch1')}
                )
            )
            with self.subTest(peptide=peptide, charge=charge):
                with self.assertRaisesRegex(ValueError,
                        f'No peptide/charge mapping found.*{peptide}, charge: {charge}'):
                    quant.remove_unquantified_peptides(self.peptide_charge_map())

    def test_remove_unquantified_peptides_checks_nonzero_unquantified_channel(self):
        record = self.peptide_quant_record(2, [10, 5])
        quant = ProteinQuantCSV(
            records=[record],
            metadata=ProteinQuantCSVMetadata(
                header_map={1: ('plex', 'ch0'), 2: ('plex', 'ch1')}
            )
        )
        with self.assertRaisesRegex(ValueError,
                'non-zero intensity value.*PEPTIDEK, charge: 2'):
            quant.remove_unquantified_peptides(self.peptide_charge_map())

    def test_protein_metadata_origin_prefix(self):
        """ Preserve the protein abundance prefix when rewriting metadata. """
        with open('test/files/test_protein_quant_csv_reheader_csv_1.csv', 'rt') as handle:
            metadata = ProteinQuantCSVMetadata.parse(handle)
        self.assertEqual(metadata.to_lines()[0],
            '# Protein abundances computed from file CPTAC_LUAD_01_merged.consensusXML'
        )

    def test_parse_header_map_protein(self):
        """ Ensure the header map is read correctly from protein csv. """
        with open('test/files/test_protein_quant_csv_reheader_csv_2.csv', 'rt') as handle:
            peptide_quant = ProteinQuantCSV.read(handle, 'protein')
        header = [
            'abundance_1',
            'abundance_2'
        ]
        for record in peptide_quant.records:
            self.assertEqual(record.header, header)

    def test_parse_header_map_peptides(self):
        """ Ensure the header map is read correctly from peptide csv. """
        with open('test/files/test_protein_quant_csv_reheader_csv_3.csv', 'rt') as handle:
            peptide_quant = ProteinQuantCSV.read(handle, 'peptide')
        header = [
            'abundance_1',
            'abundance_2'
        ]
        for record in peptide_quant.records:
            self.assertEqual(record.header, header)

    def test_openms_abundance_headers_define_plexes(self):
        """ New OpenMS CSVs encode the plex and channel in each column name. """
        source = StringIO(
            '# Peptide abundances computed from file \'merged.consensusXML\'\n'
            '# Parameters (relevant only): best_charge_and_fraction=true\n'
            '"peptide"\t"protein"\t"n_proteins"\t"charge"\t'
            '"abundance|fraction_1|ch1"\t"abundance|fraction_1|ch2"\t'
            '"abundance|fraction_2|ch1"\t"abundance|fraction_2|ch2"\t"fraction"\n'
            '"PEPTIDE"\t"PROTEIN"\t1\t2\t0\t0\t5\t0\t1\n'
        )
        quant = ProteinQuantCSV.read(source, 'peptide')
        self.assertEqual(quant.metadata.get_plex_2_column(-1),
            {'fraction_1': [0, 1], 'fraction_2': [2, 3]})
        quant.remove_zero_plex()
        self.assertTrue(np.isnan(quant.records[0].abundance[0]))
        self.assertTrue(np.isnan(quant.records[0].abundance[1]))
        self.assertEqual(quant.records[0].abundance[2:], [5.0, 0.0])

    def test_openms_protein_headers_without_metadata(self):
        source = StringIO(
            '# Protein abundances computed from file \'merged.consensusXML\'\n'
            '"protein"\t"n_proteins"\t"protein_score"\t"n_peptides"\t'
            '"abundance|sample|ch1"\t"abundance|sample|ch2"\n'
            '"P1"\t1\t0.9\t2\t0\t0\n'
        )
        quant = ProteinQuantCSV.read(source, 'protein')
        self.assertEqual(quant.metadata.get_plex_2_column(-1), {'sample': [0, 1]})
        quant.remove_zero_plex()
        self.assertTrue(all(np.isnan(x) for x in quant.records[0].abundance))

""" Test IdXML FilterProteinGroups """
import argparse
from test.integration import IntegrationTestCase
from pyopenms import IdXMLFile
from msptools import cli


class TestIdXMLFilterProteinGroups(IntegrationTestCase):
    """ Test the msptools::IdXML::filterProteinGroups cli module """
    def create_base_args(self) -> argparse.Namespace:
        """ Create base arguments """
        args = argparse.Namespace()
        args.input_path = None
        args.output_path = str(self.work_dir/'output.idXML')
        args.accept_list = None
        args.probability_cutoff = 0
        args.grouping_algorithm = 'ProteinProphet'
        args.remove_decoy = False
        args.decoy_string = 'DECOY_'
        args.decoy_string_position = 'prefix'
        return args

    def test_protein_quant_csv_reheader_case_1(self):
        """ Test IdXML filterProteinGroups case 1 """
        args = self.create_base_args()
        args.input_path = 'test/files/test_filter_protein_groups_1.idXML'
        cli.IdXML.filter_protein_groups.main(args)

        peptides = []
        proteins = []
        IdXMLFile().load(args.output_path, proteins, peptides)
        self.assertTrue(all(x.probability > 0 for x in proteins[0].getProteinGroups()))
        self.assertTrue(any(any(acc.decode('utf-8').startswith('DECOY_')
            for acc in pg.accessions) for pg in proteins[0].getProteinGroups()))

    def test_protein_quant_csv_reheader_case_2(self):
        """ Test IdXML filterProteinGroups case 2 """
        args = self.create_base_args()
        args.remove_decoy = True
        args.input_path = 'test/files/test_filter_protein_groups_1.idXML'
        args.accept_list = 'test/files/test_filter_protein_groups_accept_list.consensusXML'
        cli.IdXML.filter_protein_groups.main(args)

        peptides = []
        proteins = []
        IdXMLFile().load(args.output_path, proteins, peptides)
        self.assertTrue(all(x.probability > 0 for x in proteins[0].getProteinGroups()))
        protein_groups = [tuple(x.decode('utf-8') for x in pg.accessions)
            for pg in proteins[0].getProteinGroups()]
        self.assertTrue(all(not any(x.startswith('DECOY_') for x in pg)
            for pg in protein_groups))
        self.assertEqual(
            set(protein_groups),
            {('sp|CATA_HUMAN|',), ('ENSP00000311399.6',), ('ENSP00000503982.1',)}
        )

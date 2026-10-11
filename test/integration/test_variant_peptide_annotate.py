""" Integration ntest for ProteinQuantCSV annotate """
import argparse
from pathlib import Path
from test.integration import IntegrationTestCase
from msptools import cli


class TestProteinQuantCSVAnnotate(IntegrationTestCase):
    """ Integration ntest for ProteinQuantCSV annotate """
    def test_csv_input_1(self):
        """ test TMT """
        args = argparse.Namespace()
        args.quant_csv = Path('test/files/annotate/tmt/quant_csv.csv')
        args.idxml_csv = None
        args.database_csv = Path('test/files/annotate/tmt/split_fasta.csv')
        args.individual_fasta_csv = Path('test/files/annotate/tmt/individual_fasta.csv')
        args.gvf_csv = Path('test/files/annotate/tmt/gvf.csv')
        args.index_dir = Path('test/files/annotate/tmt/index')
        args.novel_orf_fasta = Path('test/files/annotate/tmt/noncoding.fasta')
        args.order_source = 'Coding,NovelORF,NovelORF-ALT'
        args.group_source = ['Coding:Mutation,Fusion', 'ALT:CodonReassign,SECT']
        args.allowed_sources = None
        args.label_free = False
        args.output_dir = self.work_dir
        cli.VariantPeptide.annotate.main(args)

    def test_idxml_input_1(self):
        """ test tmt with idXML """
        args = argparse.Namespace()
        args.quant_csv = None
        args.idxml_csv = Path('test/files/annotate/tmt/idxml.csv')
        args.database_csv = Path('test/files/annotate/tmt/split_fasta.csv')
        args.individual_fasta_csv = Path('test/files/annotate/tmt/individual_fasta.csv')
        args.gvf_csv = Path('test/files/annotate/tmt/gvf.csv')
        args.index_dir = Path('test/files/annotate/tmt/index')
        args.novel_orf_fasta = Path('test/files/annotate/tmt/noncoding.fasta')
        args.order_source = 'Coding,NovelORF,NovelORF-ALT'
        args.group_source = ['Coding:Mutation,Fusion', 'ALT:CodonReassign,SECT']
        args.allowed_sources = None
        args.label_free = False
        args.output_dir = self.work_dir
        cli.VariantPeptide.annotate.main(args)

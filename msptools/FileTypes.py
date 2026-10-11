""" File Types """
from enum import Enum, unique


@unique
class FileTypes(Enum):
    """ Proteomic data file tyes.
    https://github.com/OpenMS/OpenMS/blob/develop/src/openms/source/FORMAT/FileTypes.cpp
    """
    UNKNOWN           = 'unknown'      # unknown file extension
    DTA               = 'dta'          # dta raw data file
    DTA2D             = 'dta2d'        # dta2d raw data file
    MZDATA            = 'mzData'       # mzData raw data file
    MZXML             = 'mzXML'        # mzXML raw data file
    FEATUREXML        = 'featureXML'   # OpenMS feature map
    IDXML             = 'idXML'        # OpenMS peptide identification  file
    CONSENSUSXML      = 'consensusXML' # OpenMS consensus feature map
    MGF               = 'mgf'          # mascot generic format file
    INI               = 'ini'          # OpenMS parameter file
    TOPPAS            = 'toppas'       # OpenMS TOPPAS pipeline
    TRANSFORMATIONXML = 'trafoXML'     # RT transformation file
    MZML              = 'mzML'         # mzML raw data file
    CACHEDMZML        = 'cachedMzML'   # cachedMzML raw data file
    MS2               = 'ms2'          # ms2 file
    PEPXML            = 'pepXML'       # pepXML file
    PROTXML           = 'protXML'      # protXML file
    MZIDENTML         = 'mzid'         # mzIdentML file
    MZQUANTML         = 'mzq'          # mzQuantML file
    QCML              = 'qcml'         # quality control file
    MZQC              = 'mzqc'         # quality control file in json format
    GELML             = 'gelML'        # gelML file
    TRAML             = 'traML'        # transition file
    MSP               = 'msp'          # NIST spectra library file format
    OMSSAXML          = 'omssaXML'     # omssaXML file
    MASCOTXML         = 'mascotXML'    # mascotXML file
    PNG               = 'png'          # portable network graphics file
    XMASS             = 'fid'          # XMass analysis file
    TSV               = 'tsv'          # tab-separated file
    MZTAB             = 'mzTab'        # mzTab file
    PEPLIST           = 'peplist'      # SpecArray file
    HARDKLOER         = 'hardkloer'    # hardkloer file
    KROENIK           = 'kroenik'      # kroenik file
    FASTA             = 'fasta'        # FASTA file
    EDTA              = 'edta'         # enhanced dta file
    CSV               = 'csv'          # comma-separated values file
    TXT               = 'txt'          # generic text file
    OBO               = 'obo'          # controlled vocabulary file
    HTML              = 'html'         # any HTML file
    ANALYSISXML       = 'analysisXML'  # analysisXML file
    XSD               = 'xsd'          # XSD schema format
    PSQ               = 'psq'          # NCBI binary blast db
    MRM               = 'mrm'          # SpectraST MRM list
    SQMASS            = 'sqMass'       # SQLite format for mass and chromatograms
    PQP               = 'pqp'          # pqp file
    MS                = 'ms'           # SIRIUS file
    OSW               = 'osw'          # OpenSwath output files
    PSMS              = 'psms'         # Percolator tab-delimited output (PSM level)
    PIN               = 'pin'          # Percolator tab-delimited input (PSM level)
    PARAMXML          = 'paramXML'     # OpenMS internal XML file
    SPLIB             = 'splib'        # SpectraST binary spectral library file
    NOVOR             = 'novor'        # Novor custom parameter file
    XQUESTXML         = 'xquest.xml'   # xquest.xml file
    SPECXML           = 'spec.xml'     # spec.xml file
    JSON              = 'json'         # JavaScript Object Notation file
    RAW               = 'raw'          # (Thermo) Raw data file
    OMS               = 'oms'          # OpenMS SQLite file
    EXE               = 'exe'          # Windows executable
    BZ2               = 'bz2'          # bzip2 compressed file
    GZ                = 'gz'           # gzip compressed file
    XML               = 'xml'          # any XML file

    # pylint: disable=E1101
    @classmethod
    def _missing_(cls, value:str):
        """ Tries to find the file type if it's missing. This turns it into
        case insensitive. """
        name = value.replace('.', '').upper()
        if name in cls.__members__:
            return cls.__members__[name]
        return cls.UNKNOWN

    def __str__(self):
        """ str """
        return self.value

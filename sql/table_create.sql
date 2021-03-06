# standard denovo_db table
CREATE TABLE denovo_db ( 
    SampleID VARCHAR(125) NOT NULL,
    StudyName VARCHAR(125) NOT NULL,
    PubmedID VARCHAR(125) NOT NULL,
    NumProbands VARCHAR(125) NOT NULL,
    NumControls VARCHAR(125) NOT NULL,
    SequenceType VARCHAR(125) NOT NULL,
    PrimaryPhenotype VARCHAR(125) NOT NULL,
    Validation VARCHAR(125) NOT NULL,
    Chromosome VARCHAR(125) NOT NULL,
    Position BIGINT(15) NOT NULL,
    Variant VARCHAR(500) NOT NULL,
    rsID VARCHAR(125) NOT NULL,
    DbsnpBuild VARCHAR(125) NOT NULL,
    AncestralAllele VARCHAR(125) NOT NULL,
    1000GenomeCount VARCHAR(125) NOT NULL,
    ExacFreq VARCHAR(125) NOT NULL,
    EspAaFreq VARCHAR(125) NOT NULL,
    EspEaFreq VARCHAR(125) NOT NULL,
    Transcript VARCHAR(125) NOT NULL,
    codingDnaSize VARCHAR(125) NOT NULL,
    Gene VARCHAR(125) NOT NULL,
    FunctionClass VARCHAR(125) NOT NULL,
    cDnaVariant VARCHAR(125) NOT NULL,
    ProteinVariant VARCHAR(125) NOT NULL,
    Exon_Intron VARCHAR(125) NOT NULL,
    PolyPhen_HDiv VARCHAR(125) NOT NULL,
    PolyPhen_HVar VARCHAR(125) NOT NULL,
    SiftScore VARCHAR(125) NOT NULL,
    CaddScore VARCHAR(125) NOT NULL,
    LofScore VARCHAR(125) NOT NULL,
    LrtScore VARCHAR(125) NOT NULL,
    KEY (Gene),
    KEY (Transcript)
) ENGINE=InnoDB;

# filtered denovo_db table with scores
CREATE TABLE scored_denovo_db (
    PrimaryPhenotype VARCHAR(125) NOT NULL,
    Gene VARCHAR(125) NOT NULL,
    Transcript VARCHAR(125) NOT NULL,
    Chromosome VARCHAR(125) NOT NULL,
    Position BIGINT(15) NOT NULL,
    Variant VARCHAR(500) NOT NULL,
    SpliceScore DOUBLE PRECISION NOT NULL,
    PathogenScore DOUBLE PRECISION NOT NULL,
    P1Risk DOUBLE PRECISION NOT NULL,
    P2Risk DOUBLE PRECISION NOT NULL,
    P3Risk DOUBLE PRECISION NOT NULL,
    P4Risk DOUBLE PRECISION NOT NULL,
    P5Risk DOUBLE PRECISION NOT NULL,
    P6Risk DOUBLE PRECISION NOT NULL,
    P7Risk DOUBLE PRECISION NOT NULL,
    P8Risk DOUBLE PRECISION NOT NULL,
    KEY (Gene),
    KEY (Transcript),
    KEY (PrimaryPhenotype)
) ENGINE=InnoDB;


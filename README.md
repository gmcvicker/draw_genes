# draw_genes.py

This document describes how to use draw_genes.py, a flexible python program for generating figures for 
regions of the genome. For example, we used this program to make 
Figure 3 in [this paper](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1003036#pgen-1003036-g003) 
(note: some aspects of the Figure, such as labels were adjusted using Adobe Illustrator).


Once the program is configured it is easy to generate plots like this for hundreds of regions at a time.


## Setup

draw_genes.py depends on the the following:
* [PyTables](http://www.pytables.org/moin), 
* [rpy2](http://rpy.sourceforge.net/rpy2.html) and the 
* [genome python library](https://github.com/gmcvicker/genome). 

Once these are installed, obtain the source code for draw_genes.py. If you are using git, you can use the 
following command (replacing src with whatever directory you would like to use):

    mkdir ~/src
    cd ~/src
    git clone https://github.com/gmcvicker/draw_genes/

## Usage

Run `python draw_genes.py <config_file>` to generate figures.


## Configuration

draw_genes.py takes a path to a configuration file as a command line argument. 
This configuration file uses a standard syntax that can be parsed by python's 
[ConfigParser](http://docs.python.org/2/library/configparser.html). Several example configuration files 
are included in the repository under the draw_genes/conf directory. 
Let's take a look at the file draw_genes.mnase_example_regions.conf:

    [DEFAULT]
    # base directory that can be referenced by other options in this file:
    HOME=/mnt/lustre/home/gmcvicker
    # magnification used for labels in the plots
    CEX=1.0
    
    [MAIN]
    # defines the spacing between tracks
    WINDOW_MARGIN=0.4
    # if true a grid is drawn behind the tracks
    DRAW_GRID=false
    # if true a line is draw at the center of the plot
    DRAW_MIDLINE=false
    # type of regions drawn, can one of: RANDOM, GENE, COORD, BEDFILE
    REGION_TYPE=COORD
    # use this genome assembly
    ASSEMBLY=hg18
    
    
    # settings for output files
    OUTPUT_DIR=%(HOME)s/data/MNase/mnase_example_regions/
    OUTPUT_PREFIX=region_mnase_example
    OUTPUT_FORMAT=pdf
    # plot width and height in inches 
    # height of 0 means it will be set dynamically to accommodate the number of tracks
    WINDOW_WIDTH=10
    WINDOW_HEIGHT=6
    # use single or multiple output files?
    SINGLE_FILE=false
    
    # include genes at the top of the plots?
    DRAW_GENES=true

    # can output in PNG format instead of PDF. In this case dimensions are pixels.
    # OUTPUT_FORMAT=png
    # WINDOW_WIDTH=1024
    # WINDOW_HEIGHT=600
    # SINGLE_FILE=false
    
    
    # Name of genes track to use
    # GENES=CCDS
    GENES=ENSEMBL

    # Names of tracks that are plotted 
    TRACKS=CENTIPEDE,DNASE_SMOOTH10,MNASE_SMOOTH30,MNASE_SIM_PE18_SMOOTH30,KAPLAN_OCCUPANCY,KAPLAN_SCORE_SMOOTH30,YRI_READ_DEPTH
    
    
    [REGION_RANDOM]
    # Config section used if REGION_TYPE=RANDOM
    # Plots regions around randomly sampled genes
    N_REGION=10
    SEED=1234
    FLANKING=30000
    
    [REGION_GENE]
    # Config section used if REGION_TYPE=GENE
    # names of genes to plot
    GENE_NAMES=CCDS41238.1
    FLANKING=30000
    
    [REGION_COORD]
    # Config section used if REGION_TYPE=COORD
    # coordinates of regions to plot:
    COORDS=chr10:103527989-103534540,chr14:61934039-61936315,chr12:34420000-34430000
    
    
    [REGION_BEDFILE]
    # config section used if REGION_TYPE=BEDFILE
    # draw regions from the following BED file
    PATH=%(HOME)s/data/MNase/mnase_longarray/hi.hiMapRemoveHiCov.dnaseQ99.9.dist.nearest100.txt
    # expand regions to this many bp if they are smaller
    MIN_REGION_SIZE=2000
    # if this is non-zero, sample a random subset of this many regions from the bed file
    RANDOM_SUBSET=0
    SEED=1234

Hopefully the comments make it clear what most of the options are for. 
The TRACKS option in the [MAIN] section names the tracks that are plotted in each figure. 
The tracks themselves are specified in another configuration file named conf/tracks.conf. 
This file is automatically read by the draw_genes.py and contains a master list of all of 
the tracks that can be drawn.

Let's take a look inside the conf/tracks.conf configuration file. 
There are numerous subsections, most of which begin with the prefix TRACK. 
Each of these subsections describes a track that can be drawn. Here are some examples:

    [TRACK_ERNST_STATE_LCL]
    TYPE=ErnstStateTrack
    TRACK=ernst_states/Gm12878/ernst_state_label
    TRACK_LABEL=GM12878

    [TRACK_ERNST_STATE_HEPG2]
    TYPE=ErnstStateTrack
    TRACK=ernst_states/Hepg2/ernst_state_label
    TRACK_LABEL=HepG2

    [TRACK_DNASE_SMOOTH10]
    COLOR="#E41A1C"
    TYPE=ReadDepthTrack
    TRACK=dnase/dnase_all_combined
    TRACK_LABEL=DNaseI (smooth10)
    SMOOTH=10
    MIN_VAL=0
    # scale is used to put counts into units of reads / billion mapped reads
    SCALE_FACTOR=1000000000
    HEIGHT=1.0
    MAX_VAL=20
    DRAW_BORDER=false
    
    [TRACK_MNASE_SMOOTH30]
    COLOR="#377EB8"
    TYPE=ReadDepthTrack
    TRACK=mnase/q10/mnase_mids_combined
    TRACK_LABEL=MNase
    SMOOTH=30
    MIN_VAL=0
    SCALE_FACTOR=1000000000
    MAX_VAL=2.0
    HEIGHT=1.0
    DRAW_BORDER=false

    [TRACK_H3K4ME3_SMOOTH30]
    COLOR="#4DAF4A"
    TYPE=ReadDepthTrack
    TRACK=H3K4me3/H3K4me3_combined
    TRACK_LABEL=H3K4me3
    SMOOTH=30
    MIN_VAL=0
    SCALE_FACTOR=1000000000
    MAX_VAL=25
    
    [TRACK_CENTIPEDE]
    COLOR=grey
    TYPE=SegmentTrack
    TRACK=centipede_p99_max_sites
    TRACK_LABEL=Centipede
    HEIGHT=0.5
    
    [GENE_ENSEMBL]
    PATH=/data/share/genes/hg18/ens_gene.txt
    COLOR="#08306B"
    UTR_COLOR="#DEEBF7"
    TYPE=GenesTrack
    # a height of 0 indicates that it should be set dynamically
    HEIGHT=0.0

Most of these configuration sections have an option like TRACK=H3K4me3/H3K4me3_combined. 
This gives the name of a track in the genome.db that data are retrieved from (in this example 
the track is named H3K4me3/H3K4me3_combined) . This means that if you want to add a completely 
new data track, you usually have to import the data (i.e. convert it to HDF5 format) as described
in the genome.db document. It is also possible to read data from wiggle files (see ReadDepthTrack below), 
but this will cause the program to run much slower.

*Note:* tracks can also be added to the bottom of the configuration file that is specified on 
the command line. It is often more useful to add them to conf/tracks.conf, however, so 
that they can be referenced by multiple configuration files.

Each track section also has a TYPE option, which specifies which drawing class should be 
used to plot the track. Each class draws data in a different style, and has its own particular 
configuration options. The following describe the different drawing classes that can be used:

#### ReadDepthTrack
This class draws data (such as read depths) as a continuous function along the genome. By default 
data are read from the genome.db and the name of the track is specified by the TRACK option. The 
ReadDepthTrack class can also draw data from a wiggle file by providing two option lines: 
SOURCE=wig and PATH=/path/to/wiggle.txt.

#### LLRTrack
This is similar to the ReadDepthTrack, but is intended to plot a mixture of positive and 
negative values (mirrored around an axis at y=0). The positive and negative values can be 
plotted with different colors.

#### SegmentTrack
This class draws non-overlapping genomic features as segments. The start and end coordinates 
of the features are read from an HDF5 table with the specified track name. *Note:* this drawing 
class does not work with 1D arrays like some of the other classes do. Instead it uses tables 
containing start and end coordinates. These tables can be created with the load_bed.py script 
that is part of the [genome repo](https://github.com/gmcvicker/genome).

#### FeatureTrack
This class is similar to the SegmentTrack, but the features can be non-overlapping and are drawn 
less-compactly. Additionally, if DRAW_LABELS=true, then the names of the features will be plotted. 
The features are read from an HDF5 table with the specified track name. *Note:* feature tables can be 
created with the load_bed.py script that is included in the [genome repo](https://github.com/gmcvicker/genome).

### StateTrack
Plots segments and labels corresponding to a discrete set of states. The states are specified by a 
track that contains 1D arrays of of integers ranging from 0 to N_STATE. The colors and labels for each 
state must be specified as configuration options for the track.

#### ErnstStateTrack
This is an extension of the StateTrack class that is specifically for drawing Ernst chromatin states. 
The advantage of using this class is that the colors and labels for the states are set automatically 
and do not need to be specified in the configuration file.

#### GenesTrack
Genes are a special type of drawing class. The configuration section for genes does not begin with 
TRACK, but instead begins with GENE. Genes are currently read from a text file instead of from the 
database. Because UCSC provides each type of gene as a slightly different type of format, I first 
convert their files to a common format. These files are located under /data/share/genes/. I would 
like to change the way genes work so that they behave more like other tracks and so that they can 
be read from the database (this would be much faster).

#### GenotypeReadDepthTrack
This is a ReadDepthTrack that combines data across several individuals with the same genotype. 
It can be used to plot data separately for homozygous major, homozygous minor, and heterozygotes 
individuals. The genotypes for the individuals in the regions to be drawn first need to be specified. 
As a consequence, tracks with this drawing class are more complicated to setup--if you would like to 
use it, talk to me and I will help you get started.

## Future directions

I would like to change the plotting of genes so they behave more like other tracks. I would also like to
load genes into HDF5 tables so they can be read and plotted more quickly.

Currently all of the plotting is performed using rpy2. I would like to replace the drawing 
code with matplotlib.

Other drawing classes could be added as needed. I have created one that draws splice junctions,
but it needs to be updated.

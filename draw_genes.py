import sys

from ConfigParser import SafeConfigParser

import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

import genome.db
import genome.transcript
import genome.gene
import genome.coord

from draw.window import Window
from draw.genestrack import GenesTrack
from draw.llrtrack import LLRTrack
from draw.featuretrack import FeatureTrack
from draw.readdepthtrack import ReadDepthTrack
from draw.genotypereaddepthtrack import GenotypeReadDepthTrack
from draw.statetrack import StateTrack
from draw.ernststatetrack import ErnstStateTrack
from draw.segmenttrack import SegmentTrack
from draw.gccontenttrack import GCContentTrack
from draw.normreaddepthtrack import NormReadDepthTrack


grdevices = importr('grDevices')


def get_track_types():
    return {"GenesTrack" : GenesTrack,
            "LLRTrack" : LLRTrack,
            "FeatureTrack" : FeatureTrack,
            "ReadDepthTrack" : ReadDepthTrack,
            "GenotypeReadDepthTrack" : GenotypeReadDepthTrack,
            "StateTrack" : StateTrack,
            "ErnstStateTrack" : ErnstStateTrack,
            "SegmentTrack" : SegmentTrack,
            "GCContentTrack" : GCContentTrack,
            "NormReadDepthTrack" : NormReadDepthTrack}



def sample(elements, n, replace=False):    
    n_elem = len(elements)

    if n > n_elem and replace == False:
        raise ValueError("cannot sample %d elements "
                         "without replacement from list of "
                         "%d elements" % (n, n_elem))

    if replace:
        # get random indices into list of elements(repeated values allowed)
        idx = np.randint(0, n_elem, n)
    else:
        # get randomly permuted indices into list of elements
        idx = np.random.permutation(n_elem)
    
    # take elements corresponding to first n permutted indices
    return [elements[i] for i in idx[0:n]]
    



def gene_filter(g):
    """Used for filtering genes. Returns true if gene is on a non-random
    non-haplotype autosome or chrX"""
    return((g.chrom.is_x) or (g.chrom.is_auto) and not
           (g.chrom.is_rand or g.chrom.is_hap))



def get_coord_regions(config, chrom_dict):
    coord_strs = config.get("REGION_COORD", "COORDS").split(",")
    regions = []

    for coord_str in coord_strs:
        words = coord_str.split(":")
        
        if len(words) != 2:
            raise ValueError("expected region format to be chr:start-end, "
                             "got '%s'" % coord_str)
        
        chrom_name = words[0]
        if chrom_name in chrom_dict:
            chrom = chrom_dict[chrom_name]
        else:
            raise ValueError("unknown chromosome %s" % chrom_name)
        
        start_end = words[1].split("-")
        
        if len(start_end) != 2:
            raise ValueError("expected region format to be chr:start-end")
        
        start = int(start_end[0])
        end   = int(start_end[1])
        if end < start:
            raise ValueError("expected start to be >= end")
        
        region = genome.coord.Coord(chrom, start, end, strand=0)
        regions.append(region)

    return regions


    

def get_bedfile_regions(config, chrom_dict):
    path = config.get("REGION_BEDFILE", "PATH")

    min_region_size = config.getint("REGION_BEDFILE", "MIN_REGION_SIZE")

    if config.has_option("REGION_BEDFILE", "HAS_HEADER"):
        has_header = config.getboolean("REGION_BEDFILE", "HAS_HEADER")
    else:
        has_header = False

    if config.has_option("REGION_BEDFILE", "REGION_ATTRIBUTES"):
        attrib_str = config.get("REGION_BEDFILE", "REGION_ATTRIBUTES")
        region_attrib = attrib_str.split(",")
    else:
        region_attrib = []
        
    regions = genome.coord.read_bed(path, chrom_dict, 
                                    min_region_size=min_region_size,
                                    has_header=has_header,
                                    other_attrib=region_attrib)
    
    # choose random set of genes to plot
    n_rand = config.getint("REGION_BEDFILE", "RANDOM_SUBSET")
    if n_rand > 0:
        sys.stderr.write("sampling random subset of %d regions\n" %
                         n_rand)
        seed = config.getint("REGION_BEDFILE", "SEED")
        np.random.seed(seed)
        regions = sample(regions, n_rand)

    return regions
        



def get_gene_regions(config, gene_dict):
    regions = []
    flanking = config.getint("REGION_GENE", "FLANKING")
    name_str = config.get("REGION_GENE", "GENE_NAMES")
    gene_names = [name.upper() for name in name_str.split(",")]

    # build a dictionary of genes, keyed on name
    lookup_dict = {}
    for gene_list in gene_dict.values():
        for g in gene_list:
            for tr in g.transcripts:
                tr_name = tr.name.upper()
                if tr_name in gene_dict:
                    lookup_dict[tr_name.upper()].add(g)
                else:
                    lookup_dict[tr_name.upper()] = set([g])

    # search for genes with the requested names
    for gene_name in gene_names:
        if gene_name in lookup_dict:
            for g in lookup_dict[gene_name]:
                # draw a region around each gene
                region = g.copy()
                region.expand(flanking)
                regions.append(region)
        else:
            sys.stderr.write("no gene with name %s found\n" % gene_name)

    return regions



def get_rand_regions(config, gene_dict):
    # choose random set of genes to plot
    n_rand = config.getint("REGION_RANDOM", "N_REGION")
    seed = config.getint("REGION_RANDOM", "SEED")
    flanking = config.getint("REGION_RANDOM", "FLANKING")

    combined_gene_list = []
    for gene_list in gene_dict.values():
        combined_gene_list.extend(gene_list)
    
    np.random.seed(seed)
    gene_set = sample(combined_gene_list, n_rand)

    regions = []
    for g in gene_set:
        # get region around gene
        region = g.copy()
        region.expand(flanking)
        regions.append(region)

    return regions



def get_regions(config, gene_dict, chrom_dict):
    """Returns list of regions to plot, based on config options"""
    region_type = config.get("MAIN", "REGION_TYPE")

    if region_type == "COORD":
        return get_coord_regions(config, chrom_dict)
    elif region_type == "RANDOM":
        return get_rand_regions(config, gene_dict)
    elif region_type == "GENE":
        return get_gene_regions(config, gene_dict)
    elif region_type == "BEDFILE":
        return get_bedfile_regions(config, chrom_dict)
    
            
    raise ValueError("unknown region type %s" % region_type)

    return []

                

def get_genes(config, chrom_dict):
    genes_str = config.get("MAIN", "GENES")
    gene_types = genes_str.split(",")

    gene_dict = {}
    tr_dict = {}

    for gene_type in gene_types:
        genes_label = "GENE_" + gene_type
        path = config.get(genes_label, "PATH")

        sys.stderr.write("reading transcripts from %s\n" % path)
        trs = genome.transcript.read_transcripts(path, chrom_dict)

        # sort transcripts and group them into genes
        sys.stderr.write("grouping transcripts into genes\n")
        genes = genome.gene.group_transcripts(trs)
        genes = filter(gene_filter, genes)

        # resort transcripts, ignoring strand this time
        sys.stderr.write("sorting transcripts\n")
        genome.coord.sort_coords(trs, use_strand=False)

        gene_dict[genes_label] = genes
        tr_dict[genes_label] = trs

    return tr_dict, gene_dict



def main():
    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s <config_file>\n" % sys.argv[0])
        exit(2)
    
    config = SafeConfigParser()
    config.read(['conf/tracks.conf', sys.argv[1]])

    if config.has_option("MAIN", "ASSEMBLY"):
        assembly = config.get("MAIN", "ASSEMBLY")
        gdb = genome.db.GenomeDB(assembly=assembly)
    else:
        gdb = genome.db.GenomeDB()

    sys.stderr.write("using assembly %s\n" % gdb.assembly)

    r = robjects.r
    chrom_dict = gdb.get_chromosome_dict()

    if config.getboolean("MAIN", "DRAW_GENES"):
        tr_dict, gene_dict = get_genes(config, chrom_dict)
    else:
        tr_dict = {}
        gene_dict = {}
    
    regions = get_regions(config, gene_dict, chrom_dict)
    track_types = get_track_types()
    
    output_prefix = config.get("MAIN", "OUTPUT_PREFIX")
    output_dir = config.get("MAIN", "OUTPUT_DIR")
    single_file = config.getboolean("MAIN", "SINGLE_FILE")
    width = config.getfloat("MAIN", "WINDOW_WIDTH")
    output_format = config.get("MAIN", "OUTPUT_FORMAT").lower()

    if output_dir:
        if not output_dir.endswith("/"):
            output_dir = output_dir + "/"
        output_prefix = "%s%s" % (output_dir, output_prefix)

    if single_file:
        # get output file parameters
        height = config.getfloat("MAIN", "WINDOW_HEIGHT")
        if height < 5.0:
            # make minimum height 5
            height = 5.0

        if output_format == "pdf":
            filename = "%s.pdf" % output_prefix
            grdevices.pdf(file=filename, width=width, height=height)
        elif output_format == "png":
            filename = "%s.png" % output_prefix
            grdevices.png(file=filename, width=width, height=height)
        else:
            raise ValueError("unknown output format %s" % output_format)

        sys.stderr.write("writing output to single file '%s'\n" % filename)


    plot_num = 0
    for i in range(len(regions)):
        region = regions[i]
        plot_num += 1
        sys.stderr.write("DRAWING REGION %d (%s)\n" %
                         (plot_num, str(region)))

        # create window for this region
        draw_grid = config.getboolean("MAIN", "DRAW_GRID")
        if config.has_option("MAIN", "DRAW_MIDLINE"):
            draw_midline = config.getboolean("MAIN", "DRAW_MIDLINE")
        else:
            draw_midline = False
        
        margin = config.getfloat("MAIN", "WINDOW_MARGIN")
        cex = config.getfloat("MAIN", "CEX")
        window = Window(region, draw_grid=draw_grid,
                        draw_midline=draw_midline,
                        margin=margin, cex=cex)

        # add gene tracks to window
        for genes_type in gene_dict.keys():
            sys.stderr.write("  adding genes track %s\n" % genes_type)
            options = dict(config.items(genes_type))
            track_class = track_types[options['type']]
            genes_track = track_class(tr_dict[genes_type], region, options)
            window.add_track(genes_track)

        # add other tracks to window
        track_names = config.get("MAIN", "TRACKS").split(",")
        for track_name in track_names:
            if track_name.strip() == "":
                continue
            sys.stderr.write("  adding track %s\n" % track_name)

            section_name = "TRACK_" + track_name
            if not config.has_section(section_name):
                sys.stderr.write("WARNING: no config section '%s' for "
                                 "track '%s'\n" % (section_name, track_name))
                continue

            options = dict(config.items(section_name))
            # add genome db to the options, so that it can
            # be optionally used by tracks
            options['gdb'] = gdb

            if 'type' not in options:
                sys.stderr.write("WARNING: track %s does not define "
                                 "TYPE in configuration file\n" % track_name)
                continue

            track_type = options['type']

            if track_type not in track_types:
                sys.stderr.write("WARNING: don't how to create "
                                 "track %s with type %s.\n"
                                 "         Known types are %s\n"
                                 % (track_name, track_type,
                                    ", ".join(track_types.keys())))
                continue
            
            try:
                track_class = track_types[track_type]
                track = track_class(region, options)
                window.add_track(track)
            except TypeError as err:
                sys.stderr.write("WARNING: could not init track %s of "
                                 "type %s:\n%s\n" %
                                 (track_name, track_type, str(err)))
            except ValueError as err:
                sys.stderr.write("WARNING: could not open track %s of "
                                 "type %s:\n%s\n" %
                                 (track_name, options['type'], str(err)))

        if single_file:
            # each region is a separate page of a single PDF
            window.draw(r)
        else:
            # make a separate PDF for each region            
            # get output file parameters
            height = config.getfloat("MAIN", "WINDOW_HEIGHT")
            if height <= 0.0:
                height = window.get_height() * 0.5
            if height < 5.0:
                # make minimum height 5 inches
                height = 5.0

            if output_format == "pdf":
                filename = "%s%d.pdf" % (output_prefix, plot_num)
                grdevices.pdf(file=filename, width=width, height=height)
            elif output_format == "png":
                filename = "%s%d.png" % (output_prefix, plot_num)
                grdevices.png(file=filename, width=width, height=height)
            else:
                raise ValueError("unknown output format %s" % output_format)

            # render window
            window.draw(r)
            grdevices.dev_off()


    if single_file:
        grdevices.dev_off()


main()

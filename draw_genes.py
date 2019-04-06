#
# TODO:
#  - remove GenomeDB
#  - read chromosomes from chromInfo or just use
#    chromosomes specified by tracks??
#  - why are gene tracks treated differently, and not just added
#    as another track?
#  - should instantiate tracks just based on class name from config
#    should not require import or get_track_types() function
#

import sys

from configparser import ConfigParser

import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import argparse
import traceback

import genome.track
import genome.transcript
import genome.gene
import genome.coord
import genome.chrom
import genome.gtf

import region

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
from draw.pointstrack import PointsTrack

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
            "NormReadDepthTrack" : NormReadDepthTrack,
            "PointsTrack" : PointsTrack}


                

def get_genes(config, chrom_dict):
    genes_str = config.get("MAIN", "GENES")
    gene_types = genes_str.split(",")

    gene_dict = {}
    gene_list_by_gene_type = {}
    
    for gene_type in gene_types:
        genes_label = "GENE_" + gene_type
        path = config.get(genes_label, "PATH")

        gene_list, gene_dict, tr_list, tr_dict = genome.gtf.parse_gtf(path, chrom_dict)
        sys.stderr.write("sorting genes\n")
        genome.coord.sort_coords(gene_list, use_strand=False)
        gene_list_by_gene_type[genes_label] = gene_list

    return gene_list_by_gene_type



def parse_args():
    parser = argparse.ArgumentParser(description="makes plots of genomic regions")

    parser.add_argument("--tracks_file", help="path to file containing "
                        "configuration information for drawing tracks",
                        default="conf/tracks.conf")

    parser.add_argument("config_file", help="path to file containing "
                        "all other config information, including which tracks to draw")

    return parser.parse_args()
    


def main():
    args = parse_args()
        
    config = ConfigParser()
    config.read([args.tracks_file, args.config_file])

    r = robjects.r
    
    chrom_dict = genome.chrom.parse_chromosomes_dict(config.get("MAIN",
                                                                "CHROM_INFO"))

    gene_types = []
    if config.getboolean("MAIN", "DRAW_GENES"):
        genes_str = config.get("MAIN", "GENES")
        gene_types = genes_str.split(",")
        gene_dict = get_genes(config, chrom_dict)
    else:
        gene_dict = {}
    
    regions = region.get_regions(config, gene_dict, chrom_dict)
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
        reg = regions[i]
        plot_num += 1
        sys.stderr.write("DRAWING REGION %d (%s)\n" %
                         (plot_num, str(reg)))

        # create window for this region
        draw_grid = config.getboolean("MAIN", "DRAW_GRID")
        if config.has_option("MAIN", "DRAW_MIDLINE"):
            draw_midline = config.getboolean("MAIN", "DRAW_MIDLINE")
        else:
            draw_midline = False

        # draw some vertical lines on this plot?
        if config.has_option("MAIN", "DRAW_VERTLINES"):
            vert_lines = [float(x) for x in config.get("MAIN","DRAW_VERTLINES").split(",")]
            vert_lines_col = ["black"] * len(vert_lines)
        else:
            vert_lines = []
            vert_lines_col = []

        # region attributes can also be used to specify locations of
        # vertical lines
        if config.has_option("MAIN", "VERTLINES_ATTRIBUTES"):
            attr_names = config.get("MAIN", "VERTLINES_ATTRIBUTES").split(",")

            sys.stderr.write("drawing vertical lines corresponding to "
                             "region attrs %s\n" %",".join(attr_names))
            
            for a in attr_names:
                if hasattr(reg, a):
                    pos = int(getattr(reg, a))
                    vert_lines.append(pos)
                else:
                    sys.stderr.write("region is missing attribute %s\n" % a)
                
            # colors can be specified for these lines
            if config.has_option("MAIN", "VERTLINES_COLORS"):
                cols = config.get("MAIN", "VERTLINES_COLORS").split(",")
                vert_lines_col.extend(cols)

        if len(vert_lines_col) < len(vert_lines):
            # set extra lines to color black
            diff = len(vert_lines) - len(vert_lines_col)
            vert_lines_col.extend(["black"] * diff)

        # sys.stderr.write("vert_lines: %s\n" % repr(vert_lines))
        # sys.stderr.write("vert_lines_col: %s\n" % repr(vert_lines_col))
        margin = config.getfloat("MAIN", "WINDOW_MARGIN")
        cex = config.getfloat("MAIN", "CEX")
        window = Window(reg, draw_grid=draw_grid,
                        draw_midline=draw_midline,
                        vert_lines=vert_lines,
                        vert_lines_col=vert_lines_col,
                        margin=margin, cex=cex)

        # add gene tracks to window
        for genes_type in gene_types:
            gene_label = "GENE_" + genes_type 
            sys.stderr.write("  adding genes track %s\n" % gene_label)
            options = dict(config.items(gene_label))
            track_class = track_types[options['type']]
            genes_track = track_class(gene_dict[gene_label], reg, options)
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
                                    ", ".join(list(track_types.keys()))))
                continue
            
            try:
                track_class = track_types[track_type]
                track = track_class(reg, options)
                window.add_track(track)
            except TypeError as err:
                sys.stderr.write(("-" * 60) + "\n") 
                sys.stderr.write("WARNING: could not init track %s of "
                                 "type %s:\n%s\n" %
                                 (track_name, track_type, str(err)))
                traceback.print_exc()
                sys.stderr.write(("-" * 60) + "\n") 
            except ValueError as err:
                sys.stderr.write(("-" * 60) + "\n") 
                sys.stderr.write("WARNING: could not open track %s of "
                                 "type %s:\n%s\n" %
                                 (track_name, options['type'], str(err)))
                traceback.print_exc()
                sys.stderr.write(("-" * 60) + "\n") 
        
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
                
                # turn off clipping
                r.par(xpd=True)

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

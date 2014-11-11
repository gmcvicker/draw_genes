
import sys

from draw.transcripttrack import TranscriptTrack
from draw.window import Window

import genome.db
import genome.transcript

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr




if len(sys.argv) < 3:
    sys.stderr.write("usage: %s <transcript_filename> <tr_name1> "
                     "[<tr_name2> ...]\n" % 
                     sys.argv[0])
    exit(2)

tr_path = sys.argv[1]

tr_names = sys.argv[2:]


gdb = genome.db.GenomeDB()
chrom_dict = gdb.get_chromosome_dict()

sys.stderr.write("reading transcripts\n")
trs = genome.transcript.read_transcripts(tr_path, chrom_dict)
tr_dict = dict([(tr.name, tr) for tr in trs])

for tr_name in tr_names:
    if tr_name not in tr_dict:
        sys.stderr.write("WARNING: could not find transcript %s\n" % tr_name)
        continue

    tr = tr_dict[tr_name]
    
    r = robjects.r

    grdevices = importr('grDevices')

    output_format = "pdf"
    output_filename = "%s.%s" % (tr_name, output_format)
    width=8
    height=5


    sys.stderr.write("drawing transcript (filename=%s)\n" % output_filename)

    grdevices.pdf(file=output_filename, width=width, height=height)

    region = tr
    options = {'color' : "#08306B",
               'utr_color' : '#DEEBF7', 
               'border' : 'false',
               'height' : 0.1}

    window = Window(region, draw_grid=False)

    tr_track = TranscriptTrack(tr, region, options)
    window.add_track(tr_track)

    window.draw(r)

    # top = 0
    # bottom = -1
    # tr_track.set_position(tr.start, tr.end, top, bottom)

    # tr_track.draw_track(r)

    grdevices.dev_off()


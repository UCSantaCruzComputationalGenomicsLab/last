#! /usr/bin/env python

# Read pair-wise alignments in MAF or LAST tabular format: write an
# "Oxford grid", a.k.a. dotplot.

# TODO: Currently, pixels with zero aligned nt-pairs are white, and
# pixels with one or more aligned nt-pairs are black.  This can look
# too crowded for large genome alignments.  I tried shading each pixel
# according to the number of aligned nt-pairs within it, but the
# result is too faint.  How can this be done better?

import sys, os, re, itertools, optparse
import Image, ImageDraw, ImageFont, ImageColor

my_name = os.path.basename(sys.argv[0])
usage = """
  %prog --help
  %prog [options] last-tabular-output dotplot.png
  %prog [options] last-tabular-output dotplot.gif
  etc."""
parser = optparse.OptionParser(usage=usage)
# Replace "width" & "height" with a single "length" option?
parser.add_option("-x", "--width", type="int", dest="width", default=1000,
                  help="maximum width in pixels (default: %default)")
parser.add_option("-y", "--height", type="int", dest="height", default=1000,
                  help="maximum height in pixels (default: %default)")
parser.add_option("-f", "--fontfile", dest="fontfile",
                  help="TrueType or OpenType font file")
parser.add_option("-s", "--fontsize", type="int", dest="fontsize", default=11,
                  help="TrueType or OpenType font size (default: %default)")
parser.add_option("-c", "--forwardcolor", dest="forwardcolor", default="red",
                  help="Color for forward alignments (default: %default)")
parser.add_option("-r", "--reversecolor", dest="reversecolor", default="blue",
                  help="Color for reverse alignments (default: %default)")
(opts, args) = parser.parse_args()
if len(args) != 2: parser.error("2 arguments needed")

if opts.fontfile:  font = ImageFont.truetype(opts.fontfile, opts.fontsize)
else:              font = ImageFont.load_default()

# Make these options too?
text_color = "black"
background_color = "white"
pix_tween_seqs = 2  # number of border pixels between sequences
border_shade = 239, 239, 239  # the shade of grey to use for border pixels
label_space = 5     # minimum number of pixels between axis labels

image_mode = 'RGB'
forward_color = ImageColor.getcolor(opts.forwardcolor, image_mode)
reverse_color = ImageColor.getcolor(opts.reversecolor, image_mode)
overlap_color = tuple([(i+j)//2 for i, j in zip(forward_color, reverse_color)])

def isGapless(alignmentColumn):
    return "-" not in alignmentColumn

def matchAndInsertLengths(alignmentColumns):
    for k, v in itertools.groupby(alignmentColumns, isGapless):
        if k:
            matchLength = sum(1 for i in v)
            yield str(matchLength)
        else:
            blockRows = itertools.izip(*v)
            insertLengths = (len(i) - i.count("-") for i in blockRows)
            yield ":".join(map(str, insertLengths))

def alignmentInput(lines):  # read alignments in either tabular or MAF format
    for line in lines:
        w = line.split()
        if line[0].isdigit():  # tabular format
            yield w
        elif line[0] == "a":  # MAF format
            sLines = []
        elif line[0] == "s":  # MAF format
            sLines.append(w)
            if len(sLines) == 2:
                alignmentRows = (i[6] for i in sLines)
                alignmentColumns = itertools.izip(*alignmentRows)
                blocks = ",".join(matchAndInsertLengths(alignmentColumns))
                yield sLines[0][0:6] + sLines[1][1:6] + [blocks]

seq_size_dic1 = {}  # sizes of the first set of sequences
seq_size_dic2 = {}  # sizes of the second set of sequences
alignments = []

f = open(args[0])
sys.stderr.write(my_name + ": reading alignments...\n")
for w in alignmentInput(f):
    seq1, pos1, strand1, size1 = w[1], int(w[2]), w[4], int(w[5])
    seq2, pos2, strand2, size2 = w[6], int(w[7]), w[9], int(w[10])
    blocks = w[11]
    seq_size_dic1[seq1] = size1
    seq_size_dic2[seq2] = size2
    aln = seq1, seq2, pos1, pos2, strand1, strand2, blocks
    alignments.append(aln)
sys.stderr.write(my_name + ": done\n")
f.close()

if not alignments:
    sys.exit(my_name + ": there are no alignments")

def natural_sort_key(my_string):
    '''Return a sort key for "natural" ordering, e.g. chr9 < chr10.'''
    parts = re.split(r'(\d+)', my_string)
    parts[1::2] = map(int, parts[1::2])
    return parts

def get_text_sizes(my_strings):
    '''Get widths & heights, in pixels, of some strings.'''
    if opts.fontsize == 0: return [(0, 0) for i in my_strings]
    image_size = 1, 1
    im = Image.new(image_mode, image_size)
    draw = ImageDraw.Draw(im)
    return [draw.textsize(i, font=font) for i in my_strings]

def get_seq_info(seq_size_dic):
    '''Return miscellaneous information about the sequences.'''
    seq_names = seq_size_dic.keys()
    seq_names.sort(key=natural_sort_key)
    seq_sizes = [seq_size_dic[i] for i in seq_names]
    name_sizes = get_text_sizes(seq_names)
    margin = max(zip(*name_sizes)[1])  # maximum text height
    return seq_names, seq_sizes, name_sizes, margin

seq_names1, seq_sizes1, name_sizes1, margin1 = get_seq_info(seq_size_dic1)
seq_names2, seq_sizes2, name_sizes2, margin2 = get_seq_info(seq_size_dic2)

def div_ceil(x, y):
    '''Return x / y rounded up.'''
    q, r = divmod(x, y)
    return q + (r != 0)

def tot_seq_pix(seq_sizes, bp_per_pix):
    '''Return the total pixels needed for sequences of the given sizes.'''
    return sum([div_ceil(i, bp_per_pix) for i in seq_sizes])

def get_bp_per_pix(seq_sizes, pix_limit):
    '''Get the minimum bp-per-pixel that fits in the size limit.'''
    seq_num = len(seq_sizes)
    seq_pix_limit = pix_limit - pix_tween_seqs * (seq_num - 1)
    if seq_pix_limit < seq_num:
        sys.exit(my_name + ": can't fit the image: too many sequences?")
    lower_bound = div_ceil(sum(seq_sizes), seq_pix_limit)
    for bp_per_pix in itertools.count(lower_bound):  # slow linear search
        if tot_seq_pix(seq_sizes, bp_per_pix) <= seq_pix_limit: break
    return bp_per_pix

sys.stderr.write(my_name + ": choosing bp per pixel...\n")
bp_per_pix1 = get_bp_per_pix(seq_sizes1, opts.width  - margin1)
bp_per_pix2 = get_bp_per_pix(seq_sizes2, opts.height - margin2)
bp_per_pix = max(bp_per_pix1, bp_per_pix2)
sys.stderr.write(my_name + ": bp per pixel = " + str(bp_per_pix) + "\n")

def get_seq_starts(seq_pix, pix_tween_seqs, margin):
    '''Get the start pixel for each sequence.'''
    seq_starts = []
    pix_tot = margin - pix_tween_seqs
    for i in seq_pix:
        pix_tot += pix_tween_seqs
        seq_starts.append(pix_tot)
        pix_tot += i
    return seq_starts

def get_pix_info(seq_sizes, margin):
    '''Return pixel information about the sequences.'''
    seq_pix = [div_ceil(i, bp_per_pix) for i in seq_sizes]
    seq_starts = get_seq_starts(seq_pix, pix_tween_seqs, margin)
    tot_pix = seq_starts[-1] + seq_pix[-1]
    return seq_pix, seq_starts, tot_pix

seq_pix1, seq_starts1, width  = get_pix_info(seq_sizes1, margin1)
seq_pix2, seq_starts2, height = get_pix_info(seq_sizes2, margin2)
seq_start_dic1 = dict(zip(seq_names1, seq_starts1))
seq_start_dic2 = dict(zip(seq_names2, seq_starts2))
hits = [0] * (width * height)  # the image data

sys.stderr.write(my_name + ": processing alignments...\n")
for aln in alignments:
    seq1, seq2, pos1, pos2, strand1, strand2, blocks = aln
    last1 = seq_size_dic1[seq1] - 1
    last2 = seq_size_dic2[seq2] - 1
    seq_start1 = seq_start_dic1[seq1]
    seq_start2 = seq_start_dic2[seq2]
    my_start = seq_start2 * width + seq_start1
    if strand1 == strand2: store_value = 1
    else:                  store_value = 2
    for i in blocks.split(","):
        if ":" in i:  # it's a gap region: skip over it
            insertLength1, insertLength2 = i.split(":")
            pos1 += int(insertLength1)
            pos2 += int(insertLength2)
        else:  # it's a match region: draw pixels for it
            matchLength = int(i)
            end1 = pos1 + matchLength
            end2 = pos2 + matchLength
            if strand1 == '+': j = xrange(pos1, end1)
            else:              j = xrange(last1 - pos1, last1 - end1, -1)
            if strand2 == '+': k = xrange(pos2, end2)
            else:              k = xrange(last2 - pos2, last2 - end2, -1)
            for real_pos1, real_pos2 in itertools.izip(j, k):
                pix1 = real_pos1 // bp_per_pix
                pix2 = real_pos2 // bp_per_pix
                hits[my_start + pix2 * width + pix1] |= store_value
            pos1 = end1
            pos2 = end2
sys.stderr.write(my_name + ": done\n")

def make_label(text, text_size, range_start, range_size):
    '''Return an axis label with endpoint & sort-order information.'''
    text_width  = text_size[0]
    label_start = range_start + (range_size - text_width) // 2
    label_end   = label_start + text_width
    sort_key    = text_width - range_size
    return sort_key, label_start, label_end, text

def get_nonoverlapping_labels(labels):
    '''Get a subset of non-overlapping axis labels, greedily.'''
    nonoverlapping_labels = []
    for i in labels:
        if True not in [i[1] < j[2] + label_space and j[1] < i[2] + label_space
                        for j in nonoverlapping_labels]:
            nonoverlapping_labels.append(i)
    return nonoverlapping_labels

def get_axis_image(seq_names, name_sizes, seq_starts, seq_pix):
    '''Make an image of axis labels.'''
    min_pos = seq_starts[0]
    max_pos = seq_starts[-1] + seq_pix[-1]
    height = max(zip(*name_sizes)[1])
    labels = [make_label(i, j, k, l) for i, j, k, l in
              zip(seq_names, name_sizes, seq_starts, seq_pix)]
    labels = [i for i in labels if i[1] >= min_pos and i[2] <= max_pos]
    labels.sort()
    labels = get_nonoverlapping_labels(labels)
    image_size = max_pos, height
    im = Image.new(image_mode, image_size, border_shade)
    draw = ImageDraw.Draw(im)
    for i in labels:
        position = i[1], 0
        draw.text(position, i[3], font=font, fill=text_color)
    return im

image_size = width, height
im = Image.new(image_mode, image_size, background_color)

for i in range(height):
    for j in range(width):
        store_value = hits[i * width + j]
        xy = j, i
        if   store_value == 1: im.putpixel(xy, forward_color)
        elif store_value == 2: im.putpixel(xy, reverse_color)
        elif store_value == 3: im.putpixel(xy, overlap_color)

if opts.fontsize != 0:
    axis1 = get_axis_image(seq_names1, name_sizes1, seq_starts1, seq_pix1)
    axis2 = get_axis_image(seq_names2, name_sizes2, seq_starts2, seq_pix2)
    axis2 = axis2.rotate(270)
    im.paste(axis1, (0, 0))
    im.paste(axis2, (0, 0))

for i in seq_starts1[1:]:
    box = i - pix_tween_seqs, margin2, i, height
    im.paste(border_shade, box)

for i in seq_starts2[1:]:
    box = margin1, i - pix_tween_seqs, width, i
    im.paste(border_shade, box)

im.save(args[1])

#! /usr/bin/env python

# Copyright 2011, 2012, 2013 Martin C. Frith

# This script reads alignments of DNA reads to a genome, and estimates
# the probability that each alignment represents the genomic source of
# the read.  It assumes that the reads come in pairs, where each pair
# is from either end of a DNA fragment.

# Seems to work with Python 2.x, x>=4.

# The --rna option makes it assume that the genomic fragment lengths
# follow a log-normal distribution (instead of a normal distribution).
# In one test with human RNA, log-normal was a remarkably good fit,
# but not perfect.  The true distribution looked like a mixture of 2
# log-normals: a dominant one for shorter introns, and a minor one for
# huge introns.  Thus, our use of a single log-normal fails to model
# rare, huge introns.  To compensate for that, the default value of
# --disjoint is increased when --rna is used.

# (Should we try to estimate the prior probability of disjoint mapping
# from the data?  But maybe ignore low-scoring alignments for that?
# Estimate disjoint maps to opposite strands of same chromosome = maps
# to same strand of same chromosome?)

import itertools, math, operator, optparse, os, signal, sys

def logSumExp(numbers):
    """Adds numbers, in log space, to avoid overflow."""
    n = list(numbers)
    if not n: return -1e99  # should be -inf
    m = max(n)
    s = sum(math.exp(i - m) for i in n)  # fsum is only Python >= 2.6.
    return math.log(s) + m

def warn(*things):
    prog = os.path.basename(sys.argv[0])
    text = " ".join(map(str, things))
    sys.stderr.write(prog + ": " + text + "\n")

def joinby(iterable1, iterable2, keyfunc):
    """Yields pairs from iterable1 and iterable2 that share the same key."""
    groups1 = itertools.groupby(iterable1, keyfunc)
    groups2 = itertools.groupby(iterable2, keyfunc)
    k1, v1 = groups1.next()
    k2, v2 = groups2.next()
    while 1:
        if k1 < k2:
            k1, v1 = groups1.next()
        elif k1 > k2:
            k2, v2 = groups2.next()
        else:
            v2 = list(v2)
            for i1 in v1:
                for i2 in v2:
                    yield i1, i2
            k1, v1 = groups1.next()
            k2, v2 = groups2.next()

class AlignmentParameters:
    """Parses the score scale factor, minimum score, and genome size."""

    def __init__(self):  # dummy values:
        self.t = -1  # score scale factor
        self.e = -1  # minimum score
        self.g = -1  # genome size

    def update(self, line):
        for i in line.split():
            if self.t == -1 and i.startswith("t="):
                self.t = float(i[2:])
                if self.t <= 0: raise Exception("t must be positive")
            if self.e == -1 and i.startswith("e="):
                self.e = float(i[2:])
                if self.e <= 0: raise Exception("e must be positive")
            if self.g == -1 and i.startswith("letters="):
                self.g = float(i[8:])
                if self.g <= 0: raise Exception("letters must be positive")

    def isValid(self):
        return self.t != -1 and self.e != -1 and self.g != -1

    def validate(self):
        if self.t == -1: raise Exception("I need a header line with t=")
        if self.e == -1: raise Exception("I need a header line with e=")
        if self.g == -1: raise Exception("I need a header line with letters=")

def printAlignmentWithMismapProb(alignment, prob, suf):
    lines = alignment[4]
    qName = alignment[5]
    if qName.endswith("/1") or qName.endswith("/2"): suf = ""
    p = "%.3g" % prob
    if len(lines) == 1:  # we have tabular format
        w = lines[0].split()
        w[6] += suf
        w.append(p)
        print "\t".join(w)
    else:  # we have MAF format
        print lines[0].rstrip() + " mismap=" + p
        pad = " " * len(suf)  # spacer to keep the alignment of MAF lines
        rNameEnd = len(alignment[0]) + 1  # where to insert the spacer
        qNameEnd = len(qName) + 2  # where to insert the suffix
        s = 0
        for i in lines[1:]:
            if i[0] in "sq":
                if i[0] == "s": s += 1
                if s == 1:    print i[:rNameEnd] + pad + i[rNameEnd:],
                else:         print i[:qNameEnd] + suf + i[qNameEnd:],
            elif i[0] == "p": print i[:1] + pad + i[1:]
            else:             print i,
        print  # each MAF block should end with a blank line

def headToHeadDistance(alignment1, alignment2):
    """The 5'-to-5' distance between 2 alignments on opposite strands."""
    length = alignment1[1] + alignment2[1]
    if length > alignment1[2]: length -= alignment1[2]  # for circular chroms
    return length

def conjointScores(aln1, alns2, fraglen, inner, isRna):
    for i in alns2:
        length = headToHeadDistance(aln1, i)
        if isRna:  # use a log-normal distribution
            if length <= 0: continue
            loglen = math.log(length)
            yield i[3] + inner * (loglen - fraglen) ** 2 - loglen
        else:      # use a normal distribution
            if (length > 0) != (fraglen > 0): continue  # ?
            yield i[3] + inner * (length - fraglen) ** 2

def probForEachAlignment(alignments1, alignments2, opts):
    x = opts.disjointScore + logSumExp(i[3] for i in alignments2)

    fraglen = opts.fraglen
    outer = opts.outer
    inner = opts.inner
    isRna = opts.rna

    groups2 = itertools.groupby(alignments2, operator.itemgetter(0))
    genomeStrand2 = " "  # assume this is < any genomeStrand1
    for aln1 in alignments1:
        genomeStrand1 = aln1[0]
        # get the items in alignments2 that have the same genomeStrand:
        if genomeStrand2 < genomeStrand1:
            for genomeStrand2, alns2 in groups2:
                if genomeStrand2 >= genomeStrand1:
                    alns2 = list(alns2)
                    break
            else:
                genomeStrand2 = "~"  # assume this is > any genomeStrand1
        if genomeStrand1 == genomeStrand2:
            y = outer + logSumExp(conjointScores(aln1, alns2, fraglen, inner, isRna))
            yield aln1[3] + logSumExp((x, y))
        else:  # no items in alignments2 have the same genomeStrand
            yield aln1[3] + x

def printAlnsForOneRead(alignments1, alignments2, opts, maxMissingScore, suf):
    if alignments2:
        zs = list(probForEachAlignment(alignments1, alignments2, opts))
        w = maxMissingScore + max(i[3] for i in alignments2)
    else:
        zs = [i[3] + opts.disjointScore for i in alignments1]
        w = maxMissingScore

    z = logSumExp(zs)
    zw = logSumExp((z, w))

    for i, j in itertools.izip(alignments1, zs):
        prob = 1 - math.exp(j - zw)
        if prob <= opts.mismap: printAlignmentWithMismapProb(i, prob, suf)

def unambiguousFragmentLength(alignments1, alignments2):
    """Returns the fragment length implied by alignments of a pair of reads."""
    old = None
    for i, j in joinby(alignments1, alignments2, operator.itemgetter(0)):
        new = headToHeadDistance(i, j)
        if old is None: old = new
        elif new != old: return None  # the fragment length is ambiguous
    return old

def unambiguousFragmentLengths(queryPairs):
    for i, j in queryPairs:
        length = unambiguousFragmentLength(i, j)
        if length is not None: yield length

def readHeaderOrDie(lines):
    params = AlignmentParameters()
    for line in lines:
        if line[0] == "#":
            params.update(line)
            if params.isValid():
                return params
        elif not line.isspace():
            break
    params.validate()  # die

def parseAlignment(score, rName, rStart, rSpan, rSize, qName, qStrand, text,
                   strand, scale, circularChroms):
    if qStrand == strand: genomeStrand = rName + "+"
    else:                 genomeStrand = rName + "-"

    rStart = int(rStart)
    rSize = int(rSize)

    if qStrand == "+":
        c = -rStart
    else:
        c = rStart + int(rSpan)
        if rName in circularChroms or "." in circularChroms: c += rSize

    scaledScore = float(score) / scale  # needed in 2nd pass

    return genomeStrand, c, rSize, scaledScore, text, qName

def parseMafScore(aLine):
    for i in aLine.split():
        if i.startswith("score="): return i[6:]
    raise Exception("missing score")

def parseMaf(lines, strand, scale, circularChroms):
    score = parseMafScore(lines[0])
    r, q = [i.split() for i in lines if i[0] == "s"]
    return parseAlignment(score, r[1], r[2], r[3], r[5], q[1], q[4], lines,
                          strand, scale, circularChroms)

def parseTab(line, strand, scale, circularChroms):
    w = line.split()
    return parseAlignment(w[0], w[1], w[2], w[3], w[5], w[6], w[9], [line],
                          strand, scale, circularChroms)

def readBatches(lines, strand, scale, circularChroms):
    """Yields alignment data from MAF or tabular format."""
    alns = []
    maf = []
    for line in lines:
        if line[0].isdigit():
            alns.append(parseTab(line, strand, scale, circularChroms))
        elif line[0].isalpha():
            maf.append(line)
        elif line.isspace():
            if maf: alns.append(parseMaf(maf, strand, scale, circularChroms))
            maf = []
        elif line.startswith("# batch "):
            if maf: alns.append(parseMaf(maf, strand, scale, circularChroms))
            maf = []
            yield alns  # might be empty
            alns = []
    if maf: alns.append(parseMaf(maf, strand, scale, circularChroms))
    yield alns  # might be empty

def readQueryPairs(in1, in2, scale1, scale2, circularChroms):
    batches1 = readBatches(in1, "+", scale1, circularChroms)
    batches2 = readBatches(in2, "-", scale2, circularChroms)
    for i, j in itertools.izip(batches1, batches2):
        i.sort()
        j.sort()
        yield i, j

def myRound(myFloat):
    """Round a real number to a moderate amount of significant figures."""
    return float("%g" % myFloat)

def estimateFragmentLengthDistribution(lengths, opts):
    if not lengths:
        raise Exception("can't estimate the distribution of distances")

    # Define quartiles in the most naive way possible:
    lengths.sort()
    sampleSize = len(lengths)
    quartile1 = lengths[sampleSize // 4]
    quartile2 = lengths[sampleSize // 2]
    quartile3 = lengths[sampleSize * 3 // 4]

    warn("distance sample size:", sampleSize)
    warn("distance quartiles:", quartile1, quartile2, quartile3)

    if opts.rna and quartile1 <= 0:
        raise Exception("too many distances <= 0")

    if opts.rna: thing = "ln[distance]"
    else:        thing = "distance"

    if opts.fraglen is None:
        if opts.rna: opts.fraglen = myRound(math.log(quartile2))
        else:        opts.fraglen = float(quartile2)
        warn("estimated mean %s: %s" % (thing, opts.fraglen))

    if opts.sdev is None:
        if opts.rna: iqr = math.log(quartile3) - math.log(quartile1)
        else:        iqr = quartile3 - quartile1
        # Normal Distribution: sdev = iqr / (2 * qnorm(0.75))
        opts.sdev = myRound(iqr / 1.34898)
        warn("estimated standard deviation of %s: %s" % (thing, opts.sdev))

def safeLog(x):
    if x == 0: return -1e99
    else:      return math.log(x)

def calculateScorePieces(opts, params1, params2):
    if opts.sdev == 0:
        if opts.rna: opts.outer = opts.fraglen
        else:        opts.outer = 0.0
        opts.inner = -1e99
    else:  # parameters for a Normal Distribution (of fragment lengths):
        opts.outer = -math.log(opts.sdev * math.sqrt(2 * math.pi))
        opts.inner = -1.0 / (2 * opts.sdev ** 2)

    opts.outer += safeLog(1 - opts.disjoint)

    if params1.g != params2.g: raise Exception("unequal genome sizes")
    # Multiply genome size by 2, because it has 2 strands:
    opts.disjointScore = safeLog(opts.disjoint) - math.log(params1.g * 2)

    # Max possible influence of an alignment just below the score threshold:
    maxLogPrior = opts.outer
    if opts.rna: maxLogPrior += opts.sdev ** 2 / 2 - opts.fraglen
    opts.maxMissingScore1 = (params1.e - 1) / params1.t + maxLogPrior
    opts.maxMissingScore2 = (params2.e - 1) / params2.t + maxLogPrior

def lastPairProbs(opts, args):
    fileName1, fileName2 = args

    if opts.fraglen is None or opts.sdev is None:
        in1 = open(fileName1)
        in2 = open(fileName2)
        qp = readQueryPairs(in1, in2, 1, 1, opts.circular)
        lengths = list(unambiguousFragmentLengths(qp))
        estimateFragmentLengthDistribution(lengths, opts)
        in1.close()
        in2.close()

    if not opts.estdist:
        in1 = open(fileName1)
        in2 = open(fileName2)
        params1 = readHeaderOrDie(in1)
        params2 = readHeaderOrDie(in2)
        calculateScorePieces(opts, params1, params2)
        printme = opts.fraglen, opts.sdev, opts.disjoint, params1.g
        print "# fraglen=%s sdev=%s disjoint=%s genome=%.17g" % printme
        qp = readQueryPairs(in1, in2, params1.t, params2.t, opts.circular)
        for i, j in qp:
            printAlnsForOneRead(i, j, opts, opts.maxMissingScore1, "/1")
            printAlnsForOneRead(j, i, opts, opts.maxMissingScore2, "/2")
        in1.close()
        in2.close()

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # avoid silly error message

    usage = """
  %prog --help
  %prog [options] alignments1 alignments2"""

    description = "Read alignments of paired DNA reads to a genome, and: (1) estimate the distribution of distances between paired reads, (2) estimate the probability that each alignment represents the genomic source of the read."

    op = optparse.OptionParser(usage=usage, description=description)
    op.add_option("-r", "--rna", action="store_true", help=
                  "assume the reads are from potentially-spliced RNA")
    op.add_option("-e", "--estdist", action="store_true",
                  help="just estimate the distribution of distances")
    op.add_option("-m", "--mismap", type="float", default=0.01, metavar="M",
                  help="don't write alignments with mismap probability > M (default: %default)")
    op.add_option("-f", "--fraglen", type="float", metavar="BP",
                  help="mean distance in bp")
    op.add_option("-s", "--sdev", type="float", metavar="BP",
                  help="standard deviation of distance")
    op.add_option("-d", "--disjoint", type="float",
                  metavar="PROB", help=
                  "prior probability of disjoint mapping (default: 0.02 if -r, else 0.01)")
    op.add_option("-c", "--circular", action="append", metavar="CHROM",
                  help="specifies that chromosome CHROM is circular (default: chrM)")
    (opts, args) = op.parse_args()
    if opts.disjoint is None:
        if opts.rna: opts.disjoint = 0.02
        else:        opts.disjoint = 0.01
    if opts.disjoint < 0: op.error("option -d: should be >= 0")
    if opts.disjoint > 1: op.error("option -d: should be <= 1")
    if opts.sdev and opts.sdev < 0: op.error("option -s: should be >= 0")
    if len(args) != 2: op.error("please give me two file names")
    if opts.circular is None: opts.circular = ["chrM"]

    try: lastPairProbs(opts, args)
    except KeyboardInterrupt: pass  # avoid silly error message
    except Exception, e:
        warn("error:", e)
        sys.exit(1)

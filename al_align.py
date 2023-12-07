#!/usr/bin/env python
# -*- coding: utf-8 -*-
#al_finder.py
#03/2014

__author__ = 'Igor Rodrigues da Costa'
__contact__ = 'igor.bioinfo@gmail.com'

''' This program will get ALs from the genomes from a sum file,
align them and make a nexus file.'''

import os
import argparse
import shlex
from hashlib import sha1
from random import sample, choice
from subprocess import Popen
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO import _FormatToWriter
from Bio.AlignIO.NexusIO import NexusWriter
from Bio.Nexus import Nexus


class NexusWriterInterleaved(NexusWriter):
    #Set interleave to always true.
    def write_alignment(self, alignment): 
        #Creates an empty Nexus object, adds the sequences, 
        #and then gets Nexus to prepare the output. 
        if len(alignment) == 0: 
            raise ValueError("Must have at least one sequence") 
        columns = alignment.get_alignment_length() 
        if columns == 0: 
            raise ValueError("Non-empty sequences are required") 
        minimal_record = "#NEXUS\n    begin data; dimensions ntax=0 nchar=0; " \
                         + "format datatype=dna; end;"  
        n = Nexus.Nexus(minimal_record) 
        for record in alignment: 
            n.add_sequence(record.id, str(record.seq)) 
        n.write_nexus_data(self.handle, interleave=True) 


def monkey_write_nexus_data(
    self,
    filename=None,
    matrix=None,
    exclude=(),
    delete=(),
    blocksize=None,
    interleave=False,
    interleave_by_partition=False,
    comment=None,
    omit_NEXUS=False,
    append_sets=True,
    mrbayes=False,
    codons_block=True,
):
    """Write a nexus file with data and sets block to a file or handle.

    Character sets and partitions are appended by default, and are
    adjusted according to excluded characters (i.e. character sets
    still point to the same sites (not necessarily same positions),
    without including the deleted characters.

    - filename - Either a filename as a string (which will be opened,
      written to and closed), or a handle object (which will
      be written to but NOT closed).
    - interleave_by_partition - Optional name of partition (string)
    - omit_NEXUS - Boolean.  If true, the '#NEXUS' line normally at the
      start of the file is omitted.

    Returns the filename/handle used to write the data.
    """
    if not matrix:
        matrix = self.matrix
    if not matrix:
        return
    if not filename:
        filename = self.filename
    if [t for t in delete if not self._check_taxlabels(t)]:
        raise NexusError(
            "Unknown taxa: %s"
            % ", ".join(set(delete).difference(set(self.taxlabels)))
        )
    if interleave_by_partition:
        if interleave_by_partition not in self.charpartitions:
            raise NexusError("Unknown partition: %r" % interleave_by_partition)
        else:
            partition = self.charpartitions[interleave_by_partition]
            # we need to sort the partition names by starting position
            # before we exclude characters
            names = _sort_keys_by_values(partition)
            newpartition = {}
            for p in partition:
                newpartition[p] = [c for c in partition[p] if c not in exclude]
    # how many taxa and how many characters are left?
    undelete = [
        taxon for taxon in self.taxlabels if taxon in matrix and taxon not in delete
    ]
    cropped_matrix = _seqmatrix2strmatrix(
        self.crop_matrix(matrix, exclude=exclude, delete=delete)
    )
    ntax_adjusted = len(undelete)
    nchar_adjusted = len(cropped_matrix[undelete[0]])
    if not undelete or (undelete and undelete[0] == ""):
        return

    with File.as_handle(filename, mode="w") as fh:
        if not omit_NEXUS:
            fh.write("#NEXUS\n")
        if comment:
            fh.write("[" + comment + "]\n")
        fh.write("begin data;\n")
        fh.write("    dimensions ntax=%d nchar=%d;\n" % (ntax_adjusted, nchar_adjusted))
        fh.write("    format datatype=" + self.datatype)
        if self.respectcase:
            fh.write(" respectcase")
        if self.missing:
            fh.write(" missing=" + self.missing)
        if self.gap:
            fh.write(" gap=" + self.gap)
        if self.matchchar:
            fh.write(" matchchar=" + self.matchchar)
        if self.labels:
            fh.write(" labels=" + self.labels)
        if self.equate:
            fh.write(" equate=" + self.equate)
        if interleave or interleave_by_partition:
            fh.write(" interleave")
        fh.write(";\n")
        # if self.taxlabels:
        #    fh.write('taxlabels '+' '.join(self.taxlabels)+';\n')
        if self.charlabels:
            newcharlabels = self._adjust_charlabels(exclude=exclude)
            clkeys = sorted(newcharlabels)
            fh.write(
                "charlabels "
                + ", ".join(
                    "%s %s" % (k + 1, safename(newcharlabels[k])) for k in clkeys
                )
                + ";\n"
            )
        fh.write("matrix\n")
        if not blocksize:
            if interleave:
                blocksize = 70
            else:
                blocksize = self.nchar
        # delete deleted taxa and ecxclude excluded characters...
        namelength = max(len(safename(t, mrbayes=mrbayes)) for t in undelete)
        if interleave_by_partition:
            # interleave by partitions, but adjust partitions with regard
            # to excluded characters
            seek = 0
            for p in names:
                fh.write("[%s: %s]\n" % (interleave_by_partition, p))
                if len(newpartition[p]) > 0:
                    for taxon in undelete:
                        fh.write(
                            safename(taxon, mrbayes=mrbayes).ljust(namelength + 1)
                        )
                        fh.write(
                            cropped_matrix[taxon][
                                seek : seek + len(newpartition[p])
                            ]
                            + "\n"
                        )
                    fh.write("\n")
                else:
                    fh.write("[empty]\n\n")
                seek += len(newpartition[p])
        elif interleave:
            for seek in range(0, nchar_adjusted, blocksize):
                for taxon in undelete:
                    fh.write(safename(taxon, mrbayes=mrbayes).ljust(namelength + 1))
                    fh.write(cropped_matrix[taxon][seek : seek + blocksize] + "\n")
                fh.write("\n")
        else:
            for taxon in undelete:
                if blocksize < nchar_adjusted:
                    fh.write(safename(taxon, mrbayes=mrbayes) + "\n")
                else:
                    fh.write(safename(taxon, mrbayes=mrbayes).ljust(namelength + 1))
                taxon_seq = cropped_matrix[taxon]
                for seek in range(0, nchar_adjusted, blocksize):
                    fh.write(taxon_seq[seek : seek + blocksize] + "\n")
                del taxon_seq
        fh.write(";\nend;\n")
        if append_sets:
            if codons_block:
                fh.write(
                    self.append_sets(
                        exclude=exclude,
                        delete=delete,
                        mrbayes=mrbayes,
                        include_codons=False,
                    )
                )
                fh.write(
                    self.append_sets(
                        exclude=exclude,
                        delete=delete,
                        mrbayes=mrbayes,
                        codons_only=True,
                    )
                )
            else:
                fh.write(
                    self.append_sets(
                        exclude=exclude, delete=delete, mrbayes=mrbayes
                    )
                )
    return filename
        
Nexus.write_nexus_data = monkey_write_nexus_data
_FormatToWriter['nexus'] = NexusWriterInterleaved

def argument_parser(hlp = False):
    '''al_align.py'''

    default_sum = os.getcwd() + '/al_blast.sum'
    
    parser = argparse.ArgumentParser(description = 'al_align.py',\
                                     argument_default = None, fromfile_prefix_chars = '@')
    parser.add_argument('-s', '--sum', nargs = '?', type = str, required = True,\
                        dest = 'sum', help = 'Path to sum file.')
    parser.add_argument('-o', '--outpath', nargs = '?', type = str, default = os.getcwd(),\
                        dest = 'outpath', help = 'Path where the aligned results will be saved. (default: %(default)s)')
    parser.add_argument('-l', '--log', nargs = '?', type = str, default = 'al_align.log',\
                        dest = 'log', help = 'Log file. (default: %(default)s)')
    parser.add_argument('-f', '--filter', nargs = '*', type = str,\
                        dest = 'filter', help = 'Folder to look for duplicated fasta to remove. (default: %(default)s)')
    parser.add_argument('-g', '--genomes', nargs = '*', type = str,\
                        dest = 'genomes', help = 'Path to all genomes used.')
    parser.add_argument('-c', '--chromossomes', nargs = '*', type = str,\
                        dest = 'excluded', help = 'Chromossomes to be excluded.')
    parser.add_argument('-a', '--min_align', nargs = '?', type = int, default = 500,\
                        dest = 'align_size', help = 'Minimum final alignment lenght.(default: %(default)s)')
    parser.add_argument('-d', '--distance', nargs = '?', type = int, default = 200000,\
                        dest = 'idist', help = 'Minimum distance between ALs.(default: %(default)s)')
    parser.add_argument('--minsize', nargs = '?', type = int, default = 500,\
                        dest = 'min_size', help = 'Minimum sequence size.(default: %(default)s)')
    parser.add_argument('--distance_file', nargs = '?', type = str, default = 'al_align.dist',\
                        dest = 'dist_file', help = 'File to save all distances.')
    parser.add_argument('--parts', nargs = '?', type = int, default = 0,\
                        dest = 'parts', help = 'Number of parts the locus was spliced.(default: %(default)s)')
    parser.add_argument('-p', '--pick', nargs = '?', type = int,\
                        dest = 'pick', help = 'Pick only N ALs.')
    parser.add_argument('--remove_gaps', action = 'store_true', dest = 'nogaps', help = 'Remove gaps from the final alignment.')
    parser.add_argument('--chromo_sep', action = 'store_true',  dest = 'chromo_sep', help = 'Separate ALs by chromossome.')
    parser.add_argument('--skip_alignment', action = 'store_true', dest = 'skip_align', help = 'Skip alignment with clustal, only filter for distances.')
    parser.add_argument('-v', '--verbose', action = 'store_true', dest = 'verbose', help = 'Verbose switch.')
    if hlp:
        args = parser.parse_args(['-h'])
    else:
        args = parser.parse_args().__dict__
    return args

def main(args):
    #args processing:
    if args['verbose']:
        a = open(args['log'], 'w')
        a.close()
        def vprint(*a):
            # Print only when verbose
            with open(args['log'], 'a') as log:
                log.write(' '.join(a))
                log.write('\n')
    else:
        def vprint(*a):
            return None
    for k in args:
        vprint(k, ' ', str(args[k]))
    if not args['genomes']:
        print('No genomes supplied. Check usage: al_align.py --help')
        return
    vprint('#######################')
    if args['outpath'][-1] != '/':
        args['outpath'] += '/'
    #Finished analysing arguments.
    nexus_files, sizes = [], []
    args['min_seqs'] = len(args['genomes'])
    filtered_als = filter_al(args, vprint) # Remove ALs too close to each other and filters excluded chromossomes.
    seqs = read_sum(args['genomes'], args['sum'], filtered_als, vprint, args['parts'], chromo_sep=args['chromo_sep']) #seqs[AL1] = [Homo:Seq_H, Gorilla:Seq_G, Pongo:Seq_P, Pan:Seq_C]
    separated_files = write_seqs(args['outpath'], seqs, args['min_seqs'], args['min_size'], vprint, chromo_sep=args['chromo_sep'])
    if args['filter']:
        for f in args['filter']:
            separated_files = check_for_duplicates(separated_files, f, vprint)
    if args['skip_align']:
        concatenate_fasta(args['outpath'], separated_files)
        return
    aligned_files = run_clustal(args['outpath'], separated_files)
    vprint(str(len(separated_files)), ' fasta files.')
    folder = '/'.join(args['outpath'].split('/')[:-1]) + '/'
    aligned_files = [folder + filename for filename in os.listdir(folder) if '.aln' in filename]
    vprint(str(len(aligned_files)), ' aligned files.')
    if args['chromo_sep']:
        join_nexus_by_chromo(args['outpath'], aligned_files, args['align_size'], args['nogaps'], vprint)
    else: 
        sizes, nexus_files = make_nexus(args['outpath'], aligned_files, args['align_size'], args['pick'], args['nogaps'], vprint)
        if args['pick']:
            assert len(nexus_files) <= args['pick']
        join_nexus(args['outpath'], sizes, nexus_files, ntax=len(args['genomes']))
        join_fasta(args['outpath'], nexus_files)
        make_phylip(args['outpath'], nexus_files)
        with open(args['outpath'] + 'sample.log', 'w') as sample_file:
            for n in nexus_files:
                sample_file.write(n + '\n')
#    if not sizes:
#        sizes = get_sizes(nexus_files)
#    if not nexus_files:
#        nexus_files = [f for f in os.listdir(args['outpath']) if f.endswith('.nexus') and 'all' not in f]
    
def read_sum(genome_files, sumf, filtered_als, vprint, parts, chromo_sep=False):
    #will open all genomes on memory!!!
    #todo: change to SeqIO.index to improve memory usage OR change BLAST outfile to get the alignments.
    #check pyfaidx
    #al1 db1 sbj start end
    #al1 db2 sbj start end
    #al1 db3 sbj start end
    #al1 db4 sbj start end
    #al2 db1 sbj start end
    
    seqs = {} #seqs[AL1][Gorilla] = Seq_Gorilla
    genomes = {} #genomes = {genome:{chr:SEQ}}
    for g in genome_files:
        genome_name = g.split('/')[-1]
        genomes[genome_name] = {}
        vprint('opening genome: ' + genome_name)
        with open(g, 'r') as genome_file: 
            for seq_rec in SeqIO.parse(genome_file, 'fasta'):
                chromo = seq_rec.description.split()[0]
                genomes[genome_name][chromo] = seq_rec.seq
    with open(sumf, 'r') as sumf:
        for l in sumf:
            a = l.split('\t')
            al = a[0]
            if a[1] in genomes:
                genome = a[1] #UGLY! Only works if file names of db and genome are the same...
            else:
                try:
                    genome = genome_files[int(a[1].split('.')[0])]
                    genome = genome.split('/')[-1]
                except:
                    raise('FORMATDB Database must have the same name as Genome fasta') 
            chromo = a[2]
            start = int(a[3])
            end = int(a[4])
            if al in filtered_als:
                if start < end:
                    seq = genomes[genome][chromo][start:end]
                else:
                    seq = genomes[genome][chromo][end:start]
                    seq = seq.reverse_complement()
                if chromo_sep: 
                    if al in seqs:
                        seqs[al][genome] = (seq, chromo)
                    else:
                        seqs[al] = {genome:(seq, chromo)}
                else:
                    if al in seqs:
                        seqs[al][genome] = seq
                    else:
                        seqs[al] = {genome:seq}
    if parts:
        joined_seqs = {}
        for al in seqs.keys():
            al_id = al.split('.')[0]
            if al_id not in joined_seqs:
                try:
                    if chromo_sep:
                        seq = seqs[al_id+'.1']
                        for n in range(1, parts):
                            for genome in seq.keys():
                                seq[genome] = (seq[genome][0] + seqs[al_id+'.'+str(n+1)][genome][0], seq[genome][1])
                        joined_seqs[al_id] = seq
                    else:
                        seq = seqs[al_id+'.1']
                        for n in range(1, parts):
                            for genome in seq.keys():
                                seq[genome] += seqs[al_id+'.'+str(n+1)][genome]
                        joined_seqs[al_id] = seq
                except KeyError:
                    continue #some part was not found
        seqs = joined_seqs    
    return seqs

def filter_al(args, vprint):
    ''' Filter ALs based on the distance between them.'''
    al_dict = {}
    filtered_als = []
    dist = args['idist']
    last_id = ''
    all_als = 0
    removed_chromo = 0
    with open(args['sum'], 'r') as sumf:
        for l in sumf:
            line = l.split('\t')
            if line[0] != last_id: #first line in sum file must be from reference genome!
                al_id = line[0]
                al_start = int(line[3])
                al_end = int(line[4])
                al_chromo = line[2]
                last_id = al_id
                if al_chromo in al_dict:
                    al_dict[al_chromo].append((al_id, al_start, al_end, al_chromo))
                else:
                    al_dict[al_chromo] = [(al_id, al_start, al_end, al_chromo)]
    try:
        f = open(args['dist_file'], 'w')
        f.close()
    except:
        vprint('Could not open distance file: Check your permissions')
    for chromo in sorted(al_dict.keys()):
        al_list = sort_al(al_dict[chromo])
        all_als += len(al_list)
        if args['excluded'] and (chromo in args['excluded']):
            removed_chromo += len(al_list)
            continue
        end = -dist - 1 #don't exclude first al
        last_chromo = ''
        with open(args['dist_file'], 'a') as dist_file:
            for al in al_list: #al = (id, start, end, chromo)
                if end + dist < al[1] or (args['parts'] and ('.1' not in al[0])):
                    if last_chromo == al[3]:
                        dist_file.write('{0}\t{1}\n'.format(al[0], al[1] - end))
                    else:
                        dist_file.write('Chromo: ' + al[3] + '\n')
                        dist_file.write('{0}\t{1}\n'.format(al[0], al[1] - end))
                        last_chromo = al[3]
                    filtered_als.append(al[0])
                    end = al[2]
    removed_als = all_als - len(filtered_als)#[al[0] for al in all_als if al[0] not in filtered_als]
    removed_distance = removed_als - removed_chromo
    vprint(str(removed_distance) + '/' + str(all_als) + ' ALs removed due to inter-distance filter.')
    vprint(str(removed_chromo) + '/' + str(all_als) + ' ALs removed due chromossome filter.')
    return filtered_als

def sort_al(al_list):
#    def c(x, y):
#        return cmp(x[1], y[1])
    sorted_al = sorted(al_list, key=lambda c: c[1])
    return sorted_al

def write_seqs(outpath, seq_dict, minseqs, minsize, vprint, chromo_sep=False):
    files = []
    for al in seq_dict:
        seqs = seq_dict[al]
        new_file = outpath + al + '.fasta'
        records = []
        filter = True
        if len(seqs) < minseqs:
            vprint(al, str(len(seqs)), ' too few sequences')
            filter = False
        for sp in sorted(seqs): 
            if len(seqs[sp]) < minsize:
                vprint(al, str(len(seqs[sp])), ' sequence size')
                filter = False
                break
        if filter:
            for n, sp in enumerate(sorted(seqs)):
                name = sp.split('.')[0].split('_')[0]
                records.append(SeqRecord(seqs[sp], id = name, description = al)) #ID is not n if there is someone missing!
            with open(new_file, 'w') as outfile:
                SeqIO.write(records, outfile, 'fasta')
            files.append(new_file.split('/')[-1])
    return files

def chunk_reader(fobj, chunk_size=1024):
    """Generator that reads a file in chunks of bytes"""
    while True:
        chunk = fobj.read(chunk_size)
        if not chunk:
            return
        yield chunk

def check_for_duplicates(files, path, vprint, hash=sha1):
    hashes = {}
    original_len = len(files)
    for f in files:
        hashobj = hash()
        for chunk in chunk_reader(open(f, 'rb')):
            hashobj.update(chunk)
        file_id = (hashobj.digest(), os.path.getsize(f))
        hashes[file_id] = f
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            full_path = os.path.join(dirpath, filename)
            hashobj = hash()
            for chunk in chunk_reader(open(full_path, 'rb')):
                hashobj.update(chunk)
            file_id = (hashobj.digest(), os.path.getsize(full_path))
            duplicate = hashes.get(file_id, None)
            if duplicate:    
                files.remove(duplicate)
                vprint("Duplicate found: %s and %s" % (full_path, duplicate))
    return files

def run_clustal(outpath, fasta_files):
    aligned_files = []
    clustal_log = outpath + 'clustal_log.txt'
    clean = open(clustal_log, 'w')
    clean.close()
    for fp in fasta_files:
        try:
            command = 'clustalo --force -i ' + outpath+fp +\
                      ' --outfmt=FASTA -o ' + outpath+fp.split('/')[-1].replace('.fasta', '.aln')
            with open(clustal_log, 'a') as log:
                log.write(fp + ' ' + command)
                a = Popen(shlex.split(command), stdout=log, stderr=log)
                a.wait()
            aligned_files.append(fp.split('/')[-1].replace('.fasta', '.aln'))
        except OSError:
            try:
                command = 'clustalw2 -INFILE=' + outpath+fp +\
                          ' -ALIGN -OUTPUT=FASTA -OUTFILE=' + outpath+fp.split('/')[-1].replace('.fasta', '.aln')
                with open(clustal_log, 'a') as log:
                    log.write(fp + ' ' + command)
                    a = Popen(shlex.split(command), stdout=log, stderr=log)
                    a.wait()
                aligned_files.append(fp.split('/')[-1].replace('.fasta', '.aln'))
            except OSError:
                print('Clustalo not found')
                raise
    return aligned_files
                
def make_nexus(path, aligned_files, min_align_size, pick, nogaps, vprint):
    alns = [al.split('/')[-1] for al in aligned_files]
    size_dict = {}
    old_alns = []
    removed_alns = []
    nexus_files = []
    if pick:
        if pick > len(alns):
            pick = len(alns)
        old_alns = alns
        alns = sample(alns, pick)
        for alignment in alns:
            old_alns.remove(alignment)
    else:
        pick = len(alns)
        assert len(old_alns) == 0
    for i in range(pick):
        infile = alns[i]
        splited_infile = infile.split('.')
        if len(splited_infile) > 2:
            chromo = splited_infile[1]
        outfile = infile.replace('.aln', '.nexus')
        vprint('making nexus:', path + infile, path + outfile)
        size, gaps = fastatonexus(path+infile, path+outfile, keep_gaps = not nogaps, vprint = vprint)
        if size - gaps < min_align_size:
            vprint(infile, 'removed due to small alignment size.')
            removed_alns.append(infile)
            if len(old_alns) > 0: 
                random_element = choice(old_alns)
                vprint('Added ', random_element, 'to replace it.')
                old_alns.remove(random_element)
                alns.append(random_element)
        else:
            size_dict[infile] = size
            nexus_files.append(outfile)
    for alignment in removed_alns:
        alns.remove(alignment) 
    vprint(str(len(alns)), 'filtered files.')
    sizes = []
    for f in sorted(alns):
        sizes.append(size_dict[f])
    return sizes, nexus_files

def join_fasta(path, aligned_files, outfile='all.fasta'):
    aligned_files.sort()
    with open(path + outfile, 'w') as fasta_out:
        for f in aligned_files:
            f = f.replace('.nexus', '.aln')
            f_id = f.split('.')[0]
            with open(path + f, 'r') as fasta:
                for l in fasta:
                    if l.startswith('>'):
                        l = '>' + f_id + '_' + l[1:]
                    fasta_out.write(l)

def concatenate_fasta(path, files, outfile='all.fasta'):
    with open(path + outfile, 'w') as fasta_out:
        for f in files:
            f_id = f.split('.')[0]
            with open(path + f, 'r') as fasta:
                for l in fasta:
                    if l.startswith('>'):
                        l = '>' + f_id + '_' + l[1:]
                    fasta_out.write(l)

def make_phylip(path, nexus_files, concat_outfile='all.phylip'):
    nexus_files.sort()
    with open(path + concat_outfile, 'w') as out:
        for f in nexus_files:
            f = f.replace('.nexus', '.aln')
            outfile = '.'.join(f.split('.')[:-1]) + '.phylip'
            seq = fasta2phylip(path+f, path+outfile)
            out.write(seq)

def fasta2phylip(infile, outfile): 
    al_id = infile.split('.')[0].split('/')[-1]
    seqs = []
    with open(infile, 'r') as fasta:
        seq = ''
        sp_ids = []
        for l in fasta:
            if l.startswith('>'):
                l = l.strip()[1:]
                sp_name = l.split()[0]
                al_index = l.split()[1]
                sp_ids.append('{}^{}'.format(al_index, sp_name))
                if seq:
                    seqs.append(seq)
                seq = ''
            else:
                seq += l.strip()
        seqs.append(seq)
    lines = []
    with open(outfile, 'w') as out:
        first_line ='\n {}     {}\n\n'.format(len(sp_ids), len(seqs[0])) 
        out.write(first_line)
        lines.append(first_line)
        for sp in zip(sp_ids, seqs):
            line = sp[0] +'\n'+ sp[1] + '\n'
            out.write(line)
            lines.append(line)
    return ''.join(lines)

def join_nexus_by_chromo(path, aligned_files, min_align_size, nogaps, vprint):
    chromo_dict = {}
    for af in aligned_files:
        chromo = af.split('.')[1]
        if chromo in chromo_dict:
            chromo_dict[chromo].append(af)
        else:
            chromo_dict[chromo] = [af]
    for chromo in chromo_dict:
        sizes, nexus_files = make_nexus(path, chromo_dict[chromo], min_align_size, None, nogaps, vprint)
        join_nexus(path, sizes, nexus_files, 'all.%s.nexus'%chromo)

def join_nexus(path, sizes, nexus_files, outfile='all.nexus', ntax=0):
    soma_old = 0
    soma = 0
    seqs = ''
    append = '''

;
end;

begin sets;
'''
    for n, s in enumerate(sizes):
        soma = soma_old + s
        append += 'charset "locus' + str(n + 1) + '" = ' + str(soma_old + 1) + '-' + str(soma) + ';\n'
        soma_old = soma
    append += 'end;'
    start = '#NEXUS\nbegin data;\n    dimensions ntax=' + str(ntax) + ' nchar=' + str(soma) + ';\n    format datatype=dna missing=? gap=- interleave;\nmatrix\n'''
    nexus_files.sort()
    with open(path + outfile, 'w') as nex_out:
        nex_out.write(start)
        for f in nexus_files:
            with open(path + f, 'r') as nex_in:
                seqs = ''
                for l in nex_in:
                    if l[:-1] == 'matrix':
                        l = next(nex_in)
                        while l[:-1] != ';':
                            seqs += l
                            l = next(nex_in)
                nex_out.write(seqs)
        nex_out.write(append)

def get_sizes(nexus_files):
    sizes = []
    for f in sorted(nexus_files):
        i = AlignIO.read(f, 'nexus')
        sizes.append(len(i[0]))
    return sizes

def fastatonexus(infile, outfile, format_in = 'fasta', format_out = 'nexus', protein = False, keep_gaps = False, keep_n = False, vprint = lambda: None):
    align = ''
    gaps = 0
    with open(infile, 'r') as handle:
        i = AlignIO.read(handle, format_in)
        columns = len(i[0])
        for col in range(columns):
            if keep_n or 'N' not in str(i[:,col:col+1]):
                if '-' in str(i[:,col:col+1]):  
                    gaps += 1
                    if keep_gaps:
                        if align:
                            align += i[:,col:col+1]
                        else:
                            align = i[:,col:col+1]
                else:
                    if align:
                        align += i[:,col:col+1]
                    else:
                        align = i[:,col:col+1]
    if not align:
        vprint('No sequences found in', infile)
    if keep_gaps:
        size = align.get_alignment_length()
    else:
        size = align.get_alignment_length() + gaps
    vprint(infile, 'alignment size', str(size))
    align.sort()
    with open(outfile, 'w') as out:
        try:
            AlignIO.write(align, out, format_out)
        except:
            vprint('Error while saving nexus', infile)
            raise
    return size, gaps

if __name__ == '__main__':
    args = argument_parser()
    main(args)   

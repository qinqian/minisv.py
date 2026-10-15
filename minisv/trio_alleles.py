"""Candidate hypotheses from TRGT allele lengths, not mean read differences."""
import csv
import json
from pathlib import Path

import pysam


def allele_hypotheses(row, minimum=50):
    if not row['denovo_status'].startswith('Y:'):
        return []
    child_index = int(row['index'])
    child = int(row['child_AL'].split(',')[child_index])
    origin = row['allele_origin']
    parents = []
    for tag, field in [('F', 'father_AL'), ('M', 'mother_AL')]:
        for index, value in enumerate(row[field].split(','), 1):
            if origin.startswith(('F:', 'M:')) and origin != f'{tag}:{index}':
                continue
            parents.append((tag, index, int(value)))
    if not parents:
        raise ValueError('parent allele index does not exist')
    return [dict(parent=tag, parent_index=index, child_index=child_index,
                 parent_length=length, child_length=child, delta=child-length,
                 origin_resolved=origin.startswith(('F:', 'M:')))
            for tag, index, length in parents if abs(child-length) >= minimum]


def build_candidates(signal, spanning_bam, output, audit, minimum=50, indexed=False):
    inputs = {Path(p).resolve() for p in (signal, spanning_bam)}
    if Path(output).resolve() in inputs or Path(audit).resolve() in inputs or Path(output).resolve() == Path(audit).resolve():
        raise ValueError('output and audit must be distinct from each other and inputs')
    candidates = {}
    with open(signal) as source:
        header=source.readline().rstrip().split('\t')
        status_index=header.index('denovo_status')
        for line in source:
            fields=line.rstrip().split('\t')
            if not fields[status_index].startswith('Y:'):continue
            row=dict(zip(header,fields))
            if not row['denovo_status'].startswith('Y:'):
                continue
            hypotheses = allele_hypotheses(row, minimum)
            if not hypotheses:
                continue
            key = (row['trid'], int(row['index']))
            if key in candidates:
                raise ValueError(f'duplicate allele: {key}')
            candidates[key] = (row, hypotheses, set())
    with pysam.AlignmentFile(spanning_bam) as bam:
        def selected_reads():
            if not indexed:
                yield from bam.fetch(until_eof=True)
                return
            loci={key[0]:value[0] for key,value in candidates.items()}
            for trid,row in loci.items():
                for read in bam.fetch(row['chrom'],max(0,int(row['start'])-1),int(row['end'])):
                    if read.has_tag('TR') and read.get_tag('TR')==trid:yield read
        for read in selected_reads():
            if read.is_unmapped or not read.has_tag('TR') or not read.has_tag('AL'):
                continue
            key = (read.get_tag('TR'), int(read.get_tag('AL')))
            if key in candidates:
                candidates[key][2].add(read.query_name)
    with open(output, 'w') as out, open(audit, 'w') as log:
        out.write('##fileformat=VCFv4.2\n##source=minisv_parent_delta_hypotheses\n')
        for name, number, typ, desc in [
                ('SVTYPE','1','String','Parent-relative hypothesis type'),
                ('SVLEN','1','Integer','Child minus parental allele length'),
                ('END','1','Integer','Catalog interval end'),
                ('TRID','1','String','Repeat catalog identifier'),
                ('RNAMES','.','String','Allele-specific molecule identifiers'),
                ('SUPPORT','1','Integer','Unique allele read count'),
                ('DELTA_BASIS','1','String','Coordinate semantics'),
                ('PARENT_ALLELE','1','String','Candidate parental allele'),
                ('ORIGIN_RESOLVED','1','Integer','Origin assigned by TRGT')]:
            out.write(f'##INFO=<ID={name},Number={number},Type={typ},Description="{desc}">\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        count = 0
        for key, (row, hypotheses, names) in sorted(candidates.items()):
            for hyp in hypotheses:
                count += 1
                sid = f'TRGT.ALDELTA.{count}'
                log.write(json.dumps(dict(id=sid, trid=key[0], **hyp,
                          reads=len(names), status='CANDIDATE' if names else 'NO_READS'))+'\n')
                if not names:
                    continue
                typ = 'INS' if hyp['delta'] > 0 else 'DEL'
                pos = int(row['start']) + 1
                end = pos if typ == 'INS' else int(row['end'])
                info = (f"SVTYPE={typ};SVLEN={hyp['delta']};END={end};TRID={key[0]};"
                        f"RNAMES={','.join(sorted(names))};SUPPORT={len(names)};"
                        f"DELTA_BASIS=PARENT_ALLELE;PARENT_ALLELE={hyp['parent']}:{hyp['parent_index']};"
                        f"ORIGIN_RESOLVED={int(hyp['origin_resolved'])}")
                out.write(f"{row['chrom']}\t{pos}\t{sid}\tN\t<{typ}>\t.\tPASS\t{info}\n")

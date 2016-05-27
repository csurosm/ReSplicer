/*
 * The MIT License
 *
 * Copyright 2016 Mikl&oacute;s Cs&#369;r&ouml;s.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package splice;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import splice.annotation.AnnotatedGenomes;
import splice.annotation.GenePred;
import splice.sequence.DNASequence;
import splice.sequence.FastaSequence;
import splice.sequence.ProteinSequence;

/**
 *
 * Writes out DNA sequences for given set of genes: the 
 * genomic sequences cover the genes with exons and introns, as well as some overhead upstream and downstream.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class writeGenomicContext 
{
    
    private static final int OVERHANG_LENGTH = 100; 
    
    private writeGenomicContext()
    {
        
    }
    
    private void uzsgyi(String[] args) throws Exception
    {
        if (args.length==0 || args[0].equals("-h")){
            System.err.println("Tests: call as "+getClass().getName()+" org1=g1.fa,org2=g2.fa annot.txt alignment.fa");
            System.err.println("Reads gene annotations from file.");
            System.exit(9);
        }
        int arg_idx=0;
        String genome_fa_list = args[arg_idx++];
        String annotation_file = args[arg_idx++];
        String alignment_file = args[arg_idx++];
        
        AnnotatedGenomes annotations = new AnnotatedGenomes();
        List<String> wanted_organisms = annotations.readMultipleGenomes(genome_fa_list);
        annotations.readAnnotations(annotation_file);

        ProteinSequence[] protein_alignment = ProteinSequence.readMultiFasta(alignment_file);
        
        Map<String, ProteinSequence> matched_sequences = new HashMap<>();
        
        for(ProteinSequence pseq: protein_alignment)
        {
            GenePred gp = annotations.getAnnotationForGene(pseq);
            if (gp!=null)
            {
                String organism = gp.getChrom().getOrganism();
                matched_sequences.put(organism, pseq);
            }
            //System.out.println("#*wGC.uzsgyi "+pseq.getIdent()+"\torg "+pseq.getTagValue("organism")+"\tgp "+
            //        (gp==null?gp:gp.shortDesc()));
        }
        
        
        
        for (String org: wanted_organisms)
        {
            ProteinSequence pseq = matched_sequences.get(org);
            GenePred gene = annotations.getAnnotationForGene(pseq);
            DNASequence chrom = gene.getChrom();
            
            int cds_start = gene.getCDSStart();
            int cds_end = gene.getCDSEnd();
            
            int start_pos = Math.max(0, cds_start-OVERHANG_LENGTH);
            int end_pos = Math.min(chrom.getLength(), cds_end+OVERHANG_LENGTH);
            
            //System.out.println("#*wGC.uzsgyi "+pseq.getIdent()+"\torg "+pseq.getTagValue("organism")+"\tcoord "+chrom.getIdent()+gene.getStrand()+":"+cds_start+".."+cds_end);
            assert (end_pos>start_pos); // not that it matters from getSubString
            String dnaseq = chrom.getSubstring(start_pos, end_pos);
            if (gene.isReverseStrand())
            {
                dnaseq = DNASequence.reverseComplement(dnaseq);
            }
            StringBuilder defline = new StringBuilder(pseq.getDefLine());
            defline.append(" /region=").append(chrom.getIdent()).append(':').append(1+start_pos).append("..").append(end_pos);
            defline.append(" /strand=").append(gene.getStrand());
            if (pseq.getTagValue(FastaSequence.DEFLINE_TAG_ORGANISM)==null)
            {
                defline.append(" /" + FastaSequence.DEFLINE_TAG_ORGANISM + "=").append(org);
            }
            
            String fasta_seq = FastaSequence.toFasta(50, defline.toString(), dnaseq);
            System.out.println(fasta_seq);
        }
        
        
    }
    
    public static void main(String[] args) throws Exception
    {
        (new writeGenomicContext()).uzsgyi(args);
    }
}
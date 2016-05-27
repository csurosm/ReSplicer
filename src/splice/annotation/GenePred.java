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

package splice.annotation;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import splice.sequence.AminoAcid;
import splice.sequence.Codon;
import splice.sequence.DNASequence;
import static splice.sequence.DNAStrand.REVERSE;
import splice.sequence.ProteinSequence;
import splice.sequence.TranslationTable;

/**
 * Gene as a set of exons, with CDS annotation. Same information as GenePred tracks in UCSC Genome browser.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class GenePred 
{    
    public GenePred(String name, DNASequence chrom, char strand, int txStart, int txEnd) 
    {
        this.name=name;
        this.chrom=chrom;
        this.strand=strand;
        this.txStart=txStart;
        this.txEnd=txEnd;
        this.exonCount=0;
    }
    
    public GenePred(String name, DNASequence chrom, char strand, int[] exonStarts, int[] exonEnds, int cdsStart, int cdsEnd)
    {
        this.name=name;
        this.chrom=chrom;
        this.strand=strand;
        
        if (exonStarts == null)
        {
            this.exonEnds = new int[1];
            this.exonStarts = new int[1];
            
            this.exonCount =1;
            exonStarts[0] = cdsStart;
            exonEnds[0] = cdsEnd;
        } else
        {
            assert(exonEnds != null);
            
            assert(exonStarts.length == exonEnds.length);

            // check exon ends / exon starts
            List<Integer> es = new ArrayList<>();
            List<Integer> ee = new ArrayList<>();
            
            int eidx = 0;
            if (exonStarts[eidx]>exonEnds[eidx])
                throw new IllegalArgumentException("Exon length is negative: "+exonStarts[eidx]+".."+exonEnds[eidx]+"; "+name);
            int current_end = exonEnds[eidx];
            es.add(exonStarts[eidx]);
            ee.add(current_end);
            
            eidx++;
            while (eidx<exonStarts.length)
            {
                int current_start = exonStarts[eidx];
                if (current_start == current_end) // intron with length 0
                {
                    System.err.println("#WARNING("+getClass().getName()+") intron "+eidx+" length 0, exons are fused: exon_starts "+Arrays.toString(exonStarts)+", exon_ends "+Arrays.toString(exonEnds)+"; "+name);
                    current_end = exonEnds[eidx];
                    ee.set(ee.size()-1, current_end);
                } else if (current_start<current_end)
                    throw new IllegalArgumentException("Intron "+eidx+" length is negative: exon_starts "+Arrays.toString(exonStarts)+", exon_ends "+Arrays.toString(exonEnds)+"; "+name);
                else 
                {
                    if (current_start<current_end+3)
                        System.err.println("#WARNING("+getClass().getName()+") intron "+eidx+" is very short: exon_starts "+Arrays.toString(exonStarts)+", exon_ends "+Arrays.toString(exonEnds)+"; "+name);
                    current_end = exonEnds[eidx];
                    es.add(current_start);
                    ee.add(current_end);                    
                }
                eidx++;
            }
            
            assert (es.size() == ee.size());
            int[] sanitized_exon_starts;
            int[] sanitized_exon_ends;
            if (es.size() == exonStarts.length) // nothing got corrected
            {
                sanitized_exon_starts = exonStarts;
                sanitized_exon_ends = exonEnds;
            } else
            {
                sanitized_exon_starts = new int[es.size()];
                sanitized_exon_ends = new int[ee.size()];
                for (eidx=0; eidx<es.size(); eidx++)
                {
                    sanitized_exon_starts[eidx] = es.get(eidx);
                    sanitized_exon_ends[eidx] = ee.get(eidx);
                }
            }            
            
            this.exonStarts = sanitized_exon_starts;
            this.exonEnds = sanitized_exon_ends;
            this.exonCount = sanitized_exon_starts.length;
            
        }
        this.txStart = this.exonStarts[0];
        this.txEnd = this.exonEnds[exonCount-1];

        this.cdsEnd = cdsEnd;
        this.cdsStart = cdsStart;
    }
    
    private final String name; // gene name
    private final DNASequence chrom; // chromosome name
    private final char strand;
    private int txStart; // transcription start; 0-based inclusive
    private int txEnd; // transcription end;  0-based exclusive, >= txStart
    private int cdsStart; // coding sequence start; 0-based inclusive
    private int cdsEnd; // coding sequence end; 0-based exclusive
    
    private int exonCount;
    private int[] exonStarts;
    private int[] exonEnds;
    
    public String getName(){return name;}
    public DNASequence getChrom(){return chrom;}
    public char getStrand(){return strand;}
    public int getTxStart(){return txStart;}
    public int getTxEnd(){return txEnd;}
    public int getCDSStart(){return cdsStart;}
    public int getCDSEnd(){return cdsEnd;}
    public int getExonStart(int eidx){return exonStarts[eidx];}
    public int getExonEnd(int eidx){return exonEnds[eidx];}
    public int getExonCount(){return exonCount;}
    
    /**
     * Sets the column index for the chromosome column.
     * 
     * @param col 0-based column index
     */
    public static void setChromoColumn(int col)
    {
        CHROMO_COL = col;
    }
    
    /**
     * Sets the column index for the gene identifier.
     * 
     * @param col 0-based column index
     */
    public static void setAccnoColumn(int col)
    {
        ACCNO_COL = col;
    }
    
    /**
     * Sets the column index for the strand column.
     * 
     * @param col 0-based column index
     */
    public static void setStrandColumn(int col)
    {
        STRANDS_COL = col;
    }
    
    
    /**
     * Sets the column indice for the exon coordinates. Both exon start and exon ends 
     * are comma-separated positions along the given chromosome (in increasing order). 
     * Exons are encoded as half-open intervals, as in the UCSC genome browser's 
     * gene prediction tracks.  
     * 
     * @param start_col 0-based column index for exon start positions
     * @param end_col 0-based column index for exon end positions
     */
    public static void setExonColumns(int start_col, int end_col)
    {
        EXON_START_COL = start_col;
        EXON_END_COL = end_col;
    }
    
    /**
     * Sets the column index for the organism column.
     * 
     * @param col 0-based column index
     */
    public static void setOrganismColumn(int col)
    {
        ORGANISM_COL = col;
    }
    
    private static final String OPTION_CHROMO_COL = "-colidx_chromo";
    private static final String OPTION_ACCNO_COL = "-colidx_gene";
    private static final String OPTION_STRANDS_COL = "-colidx_strand";
    private static final String OPTION_EXON_COL = "-columns_exon";
    private static final String OPTION_ORGANISM_COL = "-colidx_organism";

    private static int CHROMO_COL    = 0;
    private static final int GENE_COL      = 1;
    private static int ACCNO_COL     = 2;
    private static final int CDS_COL       = 3;
    private static int STRANDS_COL   = 4;
    private static final int TRANSCRIPT_LENGTH_COL = 5;
    private static final int INTRON_POSITIONS_COL = 6;
    private static final int INTRON_LENGTH_COL = 7;
    private static final int GI_COL        = 8;
    private static int ORGANISM_COL  = 9;
    private static final int LEFTMOST_POS_COL = 10;
    private static int EXON_START_COL = 11;
    private static int EXON_END_COL = 12;
    
    /**
     * Parses one line of gene prediction for multiple genomes.
     * 
     * @param line input line 
     * @param multiple_genomes mapping between organism codes and stored genomes 
     * @return corresponding genomic exon-intron structure or null if parsing fails
     */
    public static GenePred parseTabbed(String line, Map<String,DNASequence.Genome> multiple_genomes)
    {
        String[] fields = line.split("\\t"); //
        String organism = fields[ORGANISM_COL];
        DNASequence.Genome g = multiple_genomes.get(organism);
        if (g==null) // unknown organism
        {
            return null;
        } else
            return parseTabbed(line, g);
    }
    
    /**
     * Parses one line of gene prediction for a single genomes.
     * 
     * @param line input line 
     * @param genome the parsed genome
     * @return corresponding genomic exon-intron structure or null if parsing fails
     */
    public static GenePred parseTabbed(String line, DNASequence.Genome genome)
    {
        String[] fields = line.split("\\t"); //
        // # seq   gene    accno   CDS     strand  length  ipos    ilen    gi      organism        leftmost
        
        DNASequence chromo = genome.getChromosome(fields[CHROMO_COL]);
        //System.out.println("#**GP.pT "+chromo+"\t"+fields[CHROMO_COL]);
        
        //String gene_name = fields[GI_COL];
        String gene_id = fields[ACCNO_COL];
        char strand = fields[STRANDS_COL].charAt(0);
        
        String[] xs = fields[EXON_START_COL].split(",");
        String[] xe = fields[EXON_END_COL].split(",");
        
        int num_exons = xs.length;
        while (num_exons>0 && xs[num_exons-1].equals("")) // final , 
            num_exons--;
        int[] exonStarts = new int[num_exons]; 
        int[] exonEnds = new int[num_exons];
        for (int eidx=0; eidx<num_exons; eidx++)
        {
            exonStarts[eidx]=Integer.parseInt(xs[eidx]);
            exonEnds[eidx] = Integer.parseInt(xe[eidx]);
        }
        int txStart = exonStarts[0];
        int txEnd = exonEnds[num_exons-1];
        GenePred gene = new GenePred(gene_id, chromo, strand, exonStarts, exonEnds, txStart, txEnd);
        gene.exonCount = num_exons;
//        gene.exonStarts = exonStarts;
//        gene.exonEnds = exonEnds;
        
        gene.txStart = txStart;
        gene.txEnd = txEnd;
//        gene.cdsStart = txStart;
//        gene.cdsEnd = txEnd;
        
        return gene;
    }
    
    public boolean isReverseStrand()
    {
        return getStrand() == REVERSE;
    }
    
    /**
     * Coding sequence in a region, taking strandedness into account.
     * 
     * @param region_start minimum coordinate (inclusive)
     * @param region_end maximum coordinate (exclusive)
     * @return sequence (complemented if necessary) within the regions, introns are skipped 
     */
    public String getSequence(int region_start, int region_end)
    {
        assert (region_start<region_end);
        //DNASequence ref = G.getChromosome(getChrom());
        
        StringBuilder sb = new StringBuilder();
        
        int eidx = 0;
        do
        {
            int xs = getExonStart(eidx);
            int xe = getExonEnd(eidx);
            if (xe>region_start && xs<region_end)
            {
                int pos0 = Math.max(xs, region_start);
                int pos1 = Math.min(xe, region_end);
                if (pos1 > pos0)
                {
                    String exon_seq = getChrom().getSubstring(pos0, pos1);
                    sb.append(exon_seq);
                }
            }
            ++eidx;
        } while (eidx<getExonCount());
        if (isReverseStrand())
            return DNASequence.reverseComplement(sb.toString());
        else
            return sb.toString();
    }
    
    public Codon getCDSCodonAt(int genome_pos)
    {
        int eidx = 0;
        while (eidx < getExonCount() && getExonEnd(eidx)<=genome_pos) eidx++;
        if (eidx == getExonCount()) // too far
            return null;
        if (getExonStart(eidx)>genome_pos) // intronic position or too early
            return null;
        DNASequence ref = getChrom();
        if (isReverseStrand())
        {
            byte nuc1 = DNASequence.complement(ref.getNucleotideAt(genome_pos));
            if (DNASequence.isAmbiguous(nuc1))
                return null;
            if (genome_pos == getExonStart(eidx))
            {
                --eidx;
                genome_pos = getExonEnd(eidx);
            }
            --genome_pos;
            byte nuc2 = DNASequence.complement(ref.getNucleotideAt(genome_pos));
            if (DNASequence.isAmbiguous(nuc2))
                return null;
            if (genome_pos == getExonStart(eidx))
            {
                --eidx;
                genome_pos = getExonEnd(eidx);
            }
            --genome_pos;
            byte nuc3 = DNASequence.complement(ref.getNucleotideAt(genome_pos));
            if (DNASequence.isAmbiguous(nuc3))
                return null;
            return new Codon(DNASequence.toChar(nuc1), DNASequence.toChar(nuc2), DNASequence.toChar(nuc3));
        } else
        {
            byte nuc1 = ref.getNucleotideAt(genome_pos);
            if (DNASequence.isAmbiguous(nuc1))
                return null;
            genome_pos++;
            if (genome_pos == getExonEnd(eidx))
            {
                ++eidx;
                genome_pos = getExonStart(eidx);
            }
            byte nuc2 = ref.getNucleotideAt(genome_pos);
            if (DNASequence.isAmbiguous(nuc2))
                return null;
            genome_pos++;
            if (genome_pos == getExonEnd(eidx))
            {
                ++eidx;
                genome_pos = getExonStart(eidx);
            }
            byte nuc3 = ref.getNucleotideAt(genome_pos);
            if (DNASequence.isAmbiguous(nuc3))
                return null;
            return new Codon(DNASequence.toChar(nuc1), DNASequence.toChar(nuc2), DNASequence.toChar(nuc3));
        }
        
    }
    
    public String getSequenceCDS()
    {
        return getSequence(getCDSStart(), getCDSEnd());
    }
    
    private int[] getExonPhases()
    {
        int num_exons = getExonCount();
        int[] exon_phase = new int[num_exons];
        if (isReverseStrand())
        {
            int eidx = num_exons;
            int cds_len = 0;
            do
            {
                eidx--;
                exon_phase[eidx] = (3-(cds_len % 3)) % 3;
                cds_len += getExonEnd(eidx)-getExonStart(eidx);
            } while (eidx>0);
        } else
        {
            int eidx=0;
            int cds_len = 0 ;
            while (eidx<num_exons)
            {
                exon_phase[eidx] = (3-(cds_len % 3)) % 3;
                cds_len += getExonEnd(eidx)-getExonStart(eidx);
                eidx++;
            }
        }
        return exon_phase;
    }
    
    private String intronInfoString()
    {
        StringBuilder sb_ipos = new StringBuilder("{i ");
        StringBuilder sb_ilen = new StringBuilder("/intron-length=");
        int num_exons = getExonCount();
        if (num_exons==1)
        {
            sb_ipos.append("00");
            sb_ilen.append("00");
        }
        
        if (isReverseStrand())
        {
            //
            //  <<<<<<<<A..D<<<<<<<<A..D<<<<<<<<
            //    exon 0     exon 1      exon 2
            //
            int eidx = num_exons-1;
            int cds_len = 0;
            int current_exon_last = getExonEnd(eidx);
            while (eidx>0)
            {
                int current_exon_first = getExonStart(eidx);
                int exon_len = current_exon_last-current_exon_first;
                cds_len += exon_len;
                int intron_start_cds = cds_len;
                sb_ipos.append(intron_start_cds);
                eidx--;
                int intron_end_pos=getExonEnd(eidx);
                int intron_length = current_exon_first - intron_end_pos;
                sb_ilen.append(intron_length);
                if (eidx>0)
                {
                    current_exon_last = intron_end_pos;
                    sb_ipos.append(",");
                    sb_ilen.append(",");
                }
            }
        } else
        {
            //
            //  >>>>>>>D..A>>>>>>>>D..A>>>>>>>>
            //  exon 0     exon 1      exon 2
            //
            int eidx = 0;
            int cds_len = 0;
            int current_exon_first = getExonStart(eidx);
            while (eidx+1<num_exons)
            {
                int current_exon_last = getExonEnd(eidx);
                int exon_len = current_exon_last - current_exon_first;
                cds_len += exon_len;
                sb_ipos.append(cds_len);
                eidx++;
                int intron_end_pos = getExonStart(eidx);
                int intron_length = intron_end_pos-current_exon_last;
                sb_ilen.append(intron_length);
                if (eidx+1<num_exons)
                {
                    current_exon_first = intron_end_pos;
                    sb_ipos.append(",");
                    sb_ilen.append(",");
                }
            }
        }
        sb_ipos.append(" i}");
        sb_ipos.append(" ");
        sb_ipos.append(sb_ilen);
        return sb_ipos.toString();
    }
    
    public ProteinSequence getTranslation(String defline, TranslationTable tbl)
    {
        ProteinSequence pseq = new ProteinSequence(defline+" "+intronInfoString());
        String cds = getSequenceCDS();
        
        if (cds.length() % 3 !=0)
        {
            System.err.println("#WARNING("+getClass().getName()+"): coding sequence length is not a multiple of 3, sequence will be truncated and frame may be wrong ("+defline+")");
        }
        //assert ((cds.length() % 3)==0);
        
        int pos = 0;
        while (pos+3<=cds.length())
        {
            char c1 = cds.charAt(pos++);
            char c2 = cds.charAt(pos++);
            char c3 = cds.charAt(pos++);
            if (DNASequence.isAmbiguous(c1)||DNASequence.isAmbiguous(c2)||DNASequence.isAmbiguous(c3))
            {
                pseq.append(AminoAcid.X);
            } else
            {
                Codon cod_cds = new Codon(c1, c2, c3);
                pseq.append(tbl.toAA(cod_cds));
            }
        }
        return pseq;
    }
    
    /**
     * DNA sequence for an intron
     * @param intron_idx 0-based order of the intron on the reference strand
     * @return intron sequence
     */
    public String getSequenceIntron(int intron_idx)
    {
        int is = getExonEnd(intron_idx);
        int ie = getExonStart(intron_idx+1);
        //DNASequence ref = G.getChromosome(getChrom());
        String retval
                = (isReverseStrand())
                ? getChrom().getSubstring(ie, is)
                : getChrom().getSubstring(is, ie);
        
        //System.out.println("#**GP.gSI "+is+".."+ie+"/"+getStrand()+"\tidx "+intron_idx+"\tlen "+retval.length()+"\t"+retval.substring(0,6)+".."+retval.substring(retval.length()-6,retval.length())+"\t// "+toString());
        return retval;
    }
   
    /**
     * Produces a copy where donor site is shifted
     * @param donor_pos new genomic position for donor site
     * @return an alternative transcript (with the same name)
     */
    public GenePred editDonorShift(int donor_pos)
    {
        return editSiteShift(donor_pos, true);
    }
    
    /**
     * Produces a copy where acceptor site is shifted
     * @param acceptor_pos new genomic position for acceptor site
     * @return an alternative transcript (with the same name)
     */
    public GenePred editAcceptorShift(int acceptor_pos)
    {
        return editSiteShift(acceptor_pos, false);
    }
    
    /**
     * Common code for {@link #editDonorShift(int)} and {@link #editAcceptorShift(int) }
     * @param pos new position for splice site
     * @param is_donor donor or acceptor site
     * @return modified transcript with the same name
     */
    private GenePred editSiteShift(int pos, boolean is_donor)
    {
        int num_exons=getExonCount();
        int[] new_exon_start = new int[num_exons];
        int[] new_exon_end = new int[num_exons];
        if (isReverseStrand() == is_donor)
        {
            // <<<Aa......dD<<<
            int eidx=0;

            while (getExonEnd(eidx)<=pos)
            {
                new_exon_start[eidx] = getExonStart(eidx);
                new_exon_end[eidx] = getExonEnd(eidx);
                eidx++;
            }
            new_exon_start[eidx] = pos+1;
            new_exon_end[eidx] = getExonEnd(eidx);
            eidx++;
            while (eidx<num_exons)
            {
                new_exon_start[eidx] = getExonStart(eidx);
                new_exon_end[eidx] = getExonEnd(eidx);
                eidx++;
            }
        } else
        {
            // >>>Dd...aA>>>
            int eidx=num_exons-1;
            while (getExonStart(eidx)>pos)
            {
                new_exon_start[eidx] = getExonStart(eidx);
                new_exon_end[eidx] = getExonEnd(eidx);
                eidx--;
            }
            new_exon_start[eidx] = getExonStart(eidx);
            new_exon_end[eidx] = pos;
            eidx--;
            while (eidx>=0)
            {
                new_exon_start[eidx] = getExonStart(eidx);
                new_exon_end[eidx] = getExonEnd(eidx);
                eidx--;
            }
        }
        GenePred modified_transcript = new GenePred(getName(), getChrom(), getStrand(), new_exon_start, new_exon_end, new_exon_start[0], new_exon_end[new_exon_end.length-1]);
        return modified_transcript;        
    }
    
    /**
     * Produces a copy where an intron is turned into coding sequence
     * @param donor_pos donor site position for intron start
     * @param acceptor_pos acceptor site posiiton for intronn end
     * @return modified transcript (same name)
     */
    public GenePred editExonify(int donor_pos, int acceptor_pos)
    {
        int exon_start_pos, exon_end_pos;
        if (isReverseStrand())
        {
            // <<<<Aa...dD<<<
            if (donor_pos <= acceptor_pos)
            {
                throw new IllegalArgumentException("Bad order for splice sites: donor="+donor_pos+", acceptor="+acceptor_pos+"; "+shortDesc(true));
            }
            
            assert (donor_pos>acceptor_pos);

            exon_end_pos = acceptor_pos;
            exon_start_pos = donor_pos+1;
        } else
        {
            // >>>>Dd....Aa>>>>
            exon_start_pos = acceptor_pos+1;
            exon_end_pos = donor_pos;
        }
        
        assert (exon_end_pos<exon_start_pos);
        // now exon_end_pos is the end of the exon that will be fused with the one startong at exon_start_pos
        List<Integer> es = new ArrayList<>();
        List<Integer> ee = new ArrayList<>();

        int eidx =0;
        while (exon_end_pos > getExonEnd(eidx))
        {
            es.add(getExonStart(eidx));
            ee.add(getExonEnd(eidx));
            eidx++;
        }
        assert(getExonEnd(eidx)==exon_end_pos);
        es.add(getExonStart(eidx));
        eidx++;
        while (exon_start_pos > getExonStart(eidx))
        {
            eidx++;
        }
        assert(getExonStart(eidx)==exon_start_pos);
        ee.add(getExonEnd(eidx));
        eidx++;
        while (eidx<getExonCount())
        {
            es.add(getExonStart(eidx));
            ee.add(getExonEnd(eidx));
            eidx++;
        }
        
        assert (es.size() == ee.size());
        int[] new_exon_start = new int[es.size()];
        int[] new_exon_end = new int[ee.size()];
        
        for (eidx=0; eidx<new_exon_start.length; eidx++)
        {
            new_exon_start[eidx] = es.get(eidx);
            new_exon_end[eidx] = ee.get(eidx);
        }
        
        GenePred modified_transcript = new GenePred(getName(), getChrom(), getStrand(), new_exon_start, new_exon_end, new_exon_start[0], new_exon_end[new_exon_end.length-1]);

//        System.out.println("#*GP.eEx don "+donor_pos+" acc "+acceptor_pos+"; was "+this.shortDesc(true)+" -> "+modified_transcript.shortDesc(true));
//        System.out.println("#*GP.eI oldseq "+getTranslation("old",TranslationTable.STANDARD_TABLE).toFasta(60));
//        System.out.println("#*GP.eI newseq "+modified_transcript.getTranslation("new",TranslationTable.STANDARD_TABLE).toFasta(60));
        return modified_transcript;        
    }
    
    /**
     * Produces a copy where an intron is inserted
     * 
     * @param donor_pos intron's start position
     * @param acceptor_pos intron's end position
     * @return modified transcript (same name)
     */
    public GenePred editIntronify(int donor_pos, int acceptor_pos)
    {
//        System.out.println("#*GP.eI don "+donor_pos+" acc "+acceptor_pos); System.out.flush();
        
        int exon_start_pos, exon_end_pos;
        if (isReverseStrand())
        {
            // Aa<<<<<dD
            assert (donor_pos>acceptor_pos);
            exon_end_pos = acceptor_pos;
            exon_start_pos = donor_pos + 1;
        } else
        {
            // Dd>>>aA
            assert (acceptor_pos>donor_pos);
            exon_end_pos = donor_pos;
            exon_start_pos = acceptor_pos + 1;
        }
        
        assert (exon_end_pos<exon_start_pos);
        
        // now exon_end_pos is the new end position for an exon
        // exon_start_pos is the new start position for an exon 
        List<Integer> es = new ArrayList<>();
        List<Integer> ee = new ArrayList<>();

        int eidx =0;
        while (exon_end_pos > getExonEnd(eidx))
        {
            es.add(getExonStart(eidx));
            ee.add(getExonEnd(eidx));
            eidx++;
        }
        
        assert (exon_end_pos>=getExonStart(eidx));
        es.add(getExonStart(eidx));
        ee.add(exon_end_pos);
        
        while (exon_start_pos>getExonEnd(eidx))
        {
            eidx++;
        }
        
        es.add(exon_start_pos);
        ee.add(getExonEnd(eidx));
        eidx++;
        
        while (eidx<getExonCount())
        {
            es.add(getExonStart(eidx));
            ee.add(getExonEnd(eidx));
            eidx++;
        }
        
        assert (es.size() == ee.size());
        int[] new_exon_start = new int[es.size()];
        int[] new_exon_end = new int[ee.size()];
        
        for (eidx=0; eidx<new_exon_start.length; eidx++)
        {
            new_exon_start[eidx] = es.get(eidx);
            new_exon_end[eidx] = ee.get(eidx);
        }
        
        GenePred modified_transcript = new GenePred(getName(), getChrom(), getStrand(), new_exon_start, new_exon_end, new_exon_start[0], new_exon_end[new_exon_end.length-1]);

//        System.out.println("#*GP.eIn don "+donor_pos+" acc "+acceptor_pos+"; was "+this.shortDesc(true)+" -> "+modified_transcript.shortDesc(true));
//        System.out.println("#*GP.eI oldseq "+getTranslation("old",TranslationTable.STANDARD_TABLE).toFasta(60));
//        System.out.println("#*GP.eI newseq "+modified_transcript.getTranslation("new",TranslationTable.STANDARD_TABLE).toFasta(60));
        return modified_transcript;
    }
    
    /**
     * Description for log and error messages.
     * 
     * @param want_exons whether exon start-end positions are needed
     * @return a one-line string about this gene
     */
    public String shortDesc(boolean want_exons)
    {
        StringBuilder sb = new StringBuilder();
        sb.append("GP['");
        sb.append(getName());
        sb.append("' @ ");
        sb.append(txCoordString());
        int num_exons = getExonCount();
        sb.append(" exons[");
        sb.append(num_exons);
        sb.append("]");
        if (want_exons)
        {
            sb.append("={");
            int eidx=0;
            while (eidx<num_exons)
            {
                sb.append(getExonStart(eidx));
                sb.append("..");
                sb.append(getExonEnd(eidx));
                eidx++;
                if (eidx<num_exons)
                    sb.append(",");
            } 
            sb.append("}");
        } 
        sb.append("]");
        return sb.toString();
    }
    
    /**
     * Short description of the gene.
     * @return a one-line string about this gene
     */
    public String shortDesc()
    {
        return shortDesc(false);
    }
    
    
    private static final String GFF_FEATURE_TYPE = "CDS";
    public String[] toGFF3(String source)
    {
        int num_exons = getExonCount();
        String[] gff_lines = new String[num_exons];
        int[] exon_phases = getExonPhases();
        
        for (int eidx=0; eidx<num_exons; eidx++)
        {
            StringBuilder sb = new StringBuilder();
            sb.append(getChrom().getIdent()).append("\t");
            sb.append(source).append("\t");
            sb.append(GFF_FEATURE_TYPE).append("\t");
            sb.append(getExonStart(eidx)+1).append("\t");
            sb.append(getExonEnd(eidx)).append("\t");
            sb.append(".\t"); // score
            sb.append(getStrand()).append("\t");
            sb.append(exon_phases[eidx]).append("\t");
            // attributes
            sb.append("ID=").append(getName());
            gff_lines[eidx] = sb.toString();
        }
        return gff_lines;
    }
    
    
    /**
     * @return a line in the syntax of the knownGene annotation track
     */
    @Override
    public String toString()
    {
        StringBuilder sb=new StringBuilder(name);
        sb.append('\t'); sb.append(chrom);
        sb.append('\t'); sb.append(strand);
        sb.append('\t'); sb.append(txStart);
        sb.append('\t'); sb.append(txEnd);
        sb.append('\t'); sb.append(cdsStart);
        sb.append('\t'); sb.append(cdsEnd);
        sb.append('\t'); sb.append(exonCount);
        sb.append('\t'); 
        for (int i=0; i<exonCount;i++){
            sb.append(exonStarts[i]);
            sb.append(',');
        }
        sb.append('\t');
        for (int i=0; i<exonCount; i++){
            sb.append(exonEnds[i]);
            sb.append(',');
        }
        return sb.toString();
    }
    
    /** 
     * @return a line in our own annotation syntax
     */
    public String toAnnotationString()
    {
        StringBuilder sb = new StringBuilder();
        sb.append(getChrom().getIdent()).append("\t");
        sb.append(getName()).append("\t").append(getName()).append("\t"); // name and accno
        sb.append("full\t").append(getStrand()).append("\t");
        
        int cds_len = 0;
        int num_exons = getExonCount();
        List<Integer> intron_offsets = new ArrayList<>();
        List<Integer> intron_lengths = new ArrayList<>();
        
        if (isReverseStrand())
        {
            int eidx = getExonCount()-1;
            int current_exon_last = getExonEnd(eidx);
            while (eidx>=0)
            {
                int current_exon_first = getExonStart(eidx);
                int exon_len = current_exon_last-current_exon_first;
                cds_len += exon_len;
                eidx--;
                if (eidx>=0)
                {
                    intron_offsets.add(cds_len);
                    int intron_end_pos=getExonEnd(eidx);
                    int intron_length = current_exon_first - intron_end_pos;
                    intron_lengths.add(intron_length);
                    current_exon_last = intron_end_pos;
                }
            }
            
        } else
        {
            //
            //  >>>>>>>D..A>>>>>>>>D..A>>>>>>>>
            //  exon 0     exon 1      exon 2
            //
            int eidx = 0;
            int current_exon_first = getExonStart(eidx);
            while (eidx<num_exons)
            {
                int current_exon_last = getExonEnd(eidx);
                int exon_len = current_exon_last - current_exon_first;
                cds_len += exon_len;
                
                eidx++;
                if (eidx<num_exons)
                {
                    intron_offsets.add(cds_len);
                    int intron_end_pos = getExonStart(eidx);
                    int intron_length = intron_end_pos-current_exon_last;
                    intron_lengths.add(intron_length);
                    current_exon_first = intron_end_pos;
                }
            }
        }
        
        // aa length
        if (cds_len % 3 == 0)
            sb.append(cds_len/3).append("\t");
        else
            sb.append(cds_len/3.0).append("\t");
        if (num_exons==1)
        {
            sb.append("00\t00\t");
        } else
        {
            int eidx=0;
            sb.append(intron_offsets.get(eidx));
            eidx++;
            while (eidx+1<num_exons)
            {
                sb.append(",").append(intron_offsets.get(eidx));
                eidx++;
            }
            sb.append("\t");
            eidx=0;
            sb.append(intron_lengths.get(eidx));
            eidx++;
            while (eidx+1<num_exons)
            {
                sb.append(",").append(intron_lengths.get(eidx));
                eidx++;
            }
            sb.append("\t");
        }
        sb.append(getName()).append("\t"); // gi
        
        sb.append(getChrom().getOrganism()).append("\t");
        sb.append(1+getExonStart(0)).append("\t"); // 1-based CDS start
        
        {
            int eidx=0;
            sb.append(getExonStart(eidx));
            eidx++;
            while (eidx<num_exons)
            {
                sb.append(",").append(getExonStart(eidx));
                eidx++;
            }
            eidx=0;
            sb.append("\t").append(getExonEnd(eidx));
            eidx++;
            while(eidx<num_exons)
            {
                sb.append(",").append(getExonEnd(eidx));
                eidx++;
            }
        }
        return sb.toString();
    }
    
    public String txCoordString()
    {
        StringBuilder sb=new StringBuilder();
        
//        if (!chrom.startsWith("chr"))
//            sb.append("chr");
        sb.append(chrom.getIdent());
        sb.append(".");
        sb.append(strand);
        sb.append(':');
        sb.append(strand);
        if (isReverseStrand())
        {
            sb.append(txEnd-1);
            sb.append("..");
            sb.append(txStart);
        } else
        {
            sb.append(txStart);
            sb.append("..");
            sb.append(txEnd-1);
        } 
        return sb.toString();
    }
    
    /**
     * Checks if upcoming arguments are setting parameters in this class.
     * 
     * @param args Command-line arguments
     * @param arg_idx current argument index
     * @return new argument index (larger if we processed some values)
     */
    public static int parseCommandLineOptions(String[] args, int arg_idx)
    {
        while (arg_idx<args.length && args[arg_idx].startsWith("-"))
        {
            String option = args[arg_idx++];
            String value = args[arg_idx++];
            if (OPTION_CHROMO_COL.equals(option))
            {
                setChromoColumn(Integer.parseInt(value));
            } else if (OPTION_ACCNO_COL.equals(option))
            {
                setAccnoColumn(Integer.parseInt(value));
            } else if (OPTION_STRANDS_COL.equals(option))
            {
                setStrandColumn(Integer.parseInt(value));
            } else if (OPTION_ORGANISM_COL.equals(option))
            {
                setOrganismColumn(Integer.parseInt(value));
            } else if (OPTION_EXON_COL.equals(option))
            {
                int comma_at = value.indexOf(',');
                setExonColumns(Integer.parseInt(value.substring(0, comma_at)), Integer.parseInt(value.substring(comma_at+1)));
            } else // uknnown
            {
                arg_idx -= 2;
                break;
            }
        }
        
        return arg_idx;
    }
    
       
    /** 
     * Test: reads genome sequence and gene annotations, and writes out the 
     * coding sequences.
     * @param args command-line arguments
     * @throws Exception if something goes wrong
     */
    public static void main(String[] args) throws Exception 
    {
        if (args.length==0 || args[0].equals("-h")){
            System.err.println("Tests: call as $0 <Fasta> <file>");
            System.err.println("Reads gene annotations from file for a single genome.");
            System.exit(9);
        }
        DNASequence[] refseq = DNASequence.readMultiFasta(args[0]);
        DNASequence.Genome G = new DNASequence.Genome(refseq);
        BufferedReader R=new BufferedReader(new java.io.FileReader(args[1]));
        String line ;
        do
        {
            line = R.readLine();
            if (line != null && !line.startsWith("#"))
            {
                line = line.trim();
                GenePred gene = GenePred.parseTabbed(line, G);
                String gs = gene.getSequenceCDS();
                System.out.println(">"+gene+"\n\t"+gs);
            }
        } while (line != null);
        
        R.close();
    }
}

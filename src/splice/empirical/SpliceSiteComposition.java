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

package splice.empirical;

import java.util.List;
import splice.annotation.AnnotatedGenomes;
import splice.annotation.GenePred;
import splice.model.LogScaling;
import splice.sequence.DNASequence;

/**
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class SpliceSiteComposition 
{
    public static boolean NONCANONICAL_SITES_ALLOWED = true;
    
//    public static int MOTIF_INTRONIC_LENGTH=4; //4;
//    public static int MOTIF_BOUNDARY_LENGTH=6;
    
    // >>>MMMMbbbbbb 
    // 5' I
    //                 I  5'
    //        bbbbbbMMMM<<<<
    //                 

    //                 I  3'
    //        bbbbbbMMMM>>>>
    //
    // <<<MMMMbbbbbb
    // 3' I
    
    private final Counter<Integer> motif_count;
    private final int[][] boundary_count;
    private final boolean is_donor_site;
    private final int motif_intronic_length;
    private final int motif_boundary_length;
    
    
    public SpliceSiteComposition(int motif_intronic_length, int motif_boundary_length, boolean donor_site)
    {
        this.motif_intronic_length = motif_intronic_length;
        this.motif_boundary_length = motif_boundary_length;
        
        this.motif_count = new Counter<>();
        this.boundary_count = new int[motif_boundary_length][4];
        this.is_donor_site = donor_site;
    }
    
    /**
     * Generic donor site for 98% GT, 1% GC, 1% AT.
     * 
     * @return 
     */
    public static SpliceSiteComposition.Scoring genericDonorScoring(LogScaling scaling_score)
    {
        SpliceSiteComposition ss5 = new SpliceSiteComposition(2,0,true);
        for (int rep=0; rep<98*2; rep++)
            ss5.countSite("GT", 0);
        if (NONCANONICAL_SITES_ALLOWED)
        {
            ss5.countSite("GC",0);
            ss5.countSite("AT", 0);
        }
//        ss5.countSite("GG",0);
//        ss5.countSite("GA",0);
        
        SpliceSiteComposition null_ss = genericNullSite(true);
        return ss5.getScoring(scaling_score, null_ss);
    }
    
    /**
     * Generic acceptor site for 98% GT, 1% AC, 1% TG.
     * 
     * @return 
     */
    public static SpliceSiteComposition.Scoring genericAcceptorScoring(LogScaling scaling_score)
    {
        SpliceSiteComposition ss3 = new SpliceSiteComposition(2,0,false);
        for (int rep=0; rep<98; rep++)
            ss3.countSite(DNASequence.reverseComplement("AG"), 0); // GT-AG
        if (NONCANONICAL_SITES_ALLOWED)
        {
            ss3.countSite(DNASequence.reverseComplement("AC"), 0); // AT-AC
            ss3.countSite(DNASequence.reverseComplement("TG"),0);  
        }
//        ss3.countSite("AA",0);
//        ss3.countSite("GG",0);
        
        
        
        SpliceSiteComposition null_ss = genericNullSite(false);
        SpliceSiteComposition.Scoring sc = ss3.getScoring(scaling_score, null_ss);

//        System.out.println("#*SSC.gAS\t"+Arrays.toString(sc.mot_score));
//        ss3.writeData(System.out);
//        null_ss.writeData(System.out);
//        
        
        return sc;
    }
    
    /**
     * Generic null site with uniform motif distribution
     * 
     * @param is_donor whether a donor site is needed
     * @return 
     */
    public static SpliceSiteComposition genericNullSite(boolean is_donor)
    {
        String deBruijn2 = "GGAACCTTGCATAGTCG"; // all length-2 DNA sequences
        SpliceSiteComposition null_site = new SpliceSiteComposition(2,0,is_donor);
        for (int pos=0; pos<16; pos++)
            null_site.countSite(deBruijn2, pos);
        return null_site;
    }
    
    public int[] getAllMotifCounts()
    {
        int num_motifs = 1<<(2*motif_intronic_length);
        int[] cnt = new int[num_motifs];
        
        for (int motif_idx: motif_count.keySet())
        {
            if (motif_idx>=0)
                cnt[motif_idx] = motif_count.getCount(motif_idx);
        }
        return cnt;
    }
    
    
    public Scoring getScoring(LogScaling scaling_score, SpliceSiteComposition ss_null)
    {
        assert (ss_null.motif_intronic_length==this.motif_intronic_length
                && ss_null.motif_boundary_length==this.motif_boundary_length
                && ss_null.is_donor_site == this.is_donor_site);
        return new Scoring(scaling_score, ss_null);
    }
    
    public class Scoring
    {
        private Scoring(LogScaling scoring_scale, SpliceSiteComposition ss_null)
        {
            int[] ss_mot1 = getAllMotifCounts();
            int[] ss_mot0 = ss_null.getAllMotifCounts();
            mot_score = scoring_scale.getScaledScore(ss_mot1, ss_mot0);
            //System.out.println("#*SSC.S m1 "+Arrays.toString(ss_mot1)+"\t"+Arrays.toString(ss_mot0)+"\tsc "+Arrays.toString(mot_score));
            bnd_score = new int[motif_boundary_length][];
            for (int j=0; j<motif_boundary_length; j++)
                bnd_score[j] = scoring_scale.getScaledScore(boundary_count[j], ss_null.boundary_count[j]);
        }
        
        private final int[] mot_score;
        private final int[][] bnd_score;
        
        public int getSiteLength()
        {
            return motif_intronic_length+motif_boundary_length;
        }
        
        public int getScore(DNASequence seq, int pos, boolean is_forward_strand)
        {
            if (is_forward_strand == is_donor_site)
            {
                if (pos+motif_intronic_length+motif_boundary_length > seq.getLength())
                    return Integer.MIN_VALUE;
                int motif_idx = encodeMotif(seq, pos, true);
                int sc = (motif_idx == -1?Integer.MIN_VALUE:mot_score[motif_idx]);
                int bpos = pos+motif_intronic_length;
                for (int bidx=0; bidx<motif_boundary_length && sc!=Integer.MIN_VALUE; bidx++, bpos++)
                {
                    byte nuc_fw = seq.getNucleotideAt(bpos);
                    if (DNASequence.isAmbiguous(nuc_fw))
                        sc = Integer.MIN_VALUE;
                    else
                    {
                        int nuc_idx = DNASequence.toIndex(nuc_fw);
                        sc += bnd_score[bidx][nuc_idx];
                    }
                }
                return sc;
            } else
            {
                if (pos-motif_intronic_length-motif_boundary_length < -1)
                    return Integer.MIN_VALUE;
                int motif_idx = encodeMotif(seq, pos, false);
                int sc = (motif_idx == -1?Integer.MIN_VALUE:mot_score[motif_idx]);
                int bpos = pos-motif_intronic_length;
                for (int bidx=0; bidx<motif_boundary_length && sc!=Integer.MIN_VALUE; bidx++, bpos--)
                {
                    byte nuc_fw = seq.getNucleotideAt(bpos);
                    if (DNASequence.isAmbiguous(nuc_fw))
                        sc = Integer.MIN_VALUE;
                    else
                    {
                        int nuc_idx = DNASequence.toIndex(DNASequence.complement(nuc_fw));
                        sc += bnd_score[bidx][nuc_idx];
                    }
                }
                return sc;
            }
        }
        
        public int getScore(String seq, int pos, boolean is_forward_strand)
        {
            if (is_donor_site == is_forward_strand)
            {
                if (pos+motif_intronic_length+motif_boundary_length > seq.length())
                    return Integer.MIN_VALUE;
                int motif_idx = encodeMotif(seq, pos, true);
                int sc = (motif_idx == -1?Integer.MIN_VALUE:mot_score[motif_idx]);
                int bpos = pos+motif_intronic_length;
                for (int bidx=0; bidx<motif_boundary_length && sc!=Integer.MIN_VALUE; bidx++, bpos++)
                {
                    char nuc_fw = seq.charAt(bpos);
                    if (DNASequence.isAmbiguous(nuc_fw))
                        sc = Integer.MIN_VALUE;
                    else
                    {
                        int nuc_idx = DNASequence.toIndex(nuc_fw);
                        sc += bnd_score[bidx][nuc_idx];
                    }
                }
                return sc;
            } else
            {
                if (pos-motif_intronic_length-motif_boundary_length < -1)
                    return Integer.MIN_VALUE;
                int motif_idx = encodeMotif(seq, pos, false);
                int sc = (motif_idx == -1?Integer.MIN_VALUE:mot_score[motif_idx]);
                int bpos = pos-motif_intronic_length;
                for (int bidx=0; bidx<motif_boundary_length && sc!=Integer.MIN_VALUE; bidx++, bpos--)
                {
                    char nuc_fw = seq.charAt(bpos);
                    if (DNASequence.isAmbiguous(nuc_fw))
                        sc = Integer.MIN_VALUE;
                    else
                    {
                        int nuc_idx = DNASequence.toIndex(DNASequence.complement(nuc_fw));
                        sc += bnd_score[bidx][nuc_idx];
                    }
                }
                return sc;
            }
        }        
    }
    
    public int[][] getBoundaryCounts()
    {
//        int[][] cnt = new int[motif_boundary_length][4];
//        if (is_donor_site)
//        {
//            for (int i=0; i<motif_boundary_length; i++)
//            {
//                cnt[i][NucleotideEncoding.IDX_A] = boundary_count[i][NucleotideEncoding.IDX_A];
//                cnt[i][NucleotideEncoding.IDX_C] = boundary_count[i][NucleotideEncoding.IDX_C];
//                cnt[i][NucleotideEncoding.IDX_G] = boundary_count[i][NucleotideEncoding.IDX_G];
//                cnt[i][NucleotideEncoding.IDX_T] = boundary_count[i][NucleotideEncoding.IDX_T];
//                
//            }
//        } else
//        {
//            for (int i=0; i<motif_boundary_length; i++)
//            {
//                cnt[i][NucleotideEncoding.IDX_A] = boundary_count[i][NucleotideEncoding.IDX_T];
//                cnt[i][NucleotideEncoding.IDX_C] = boundary_count[i][NucleotideEncoding.IDX_G];
//                cnt[i][NucleotideEncoding.IDX_G] = boundary_count[i][NucleotideEncoding.IDX_C];
//                cnt[i][NucleotideEncoding.IDX_T] = boundary_count[i][NucleotideEncoding.IDX_A];
//            }
//        }
//        return cnt;
        return boundary_count;
    }
    
    public void countIntronSites(GenePred gene)
    {
        int num_exons = gene.getExonCount();
        DNASequence ref = gene.getChrom();
        if (gene.isReverseStrand())
        {
            for (int eidx=1; eidx<num_exons; eidx++)
            {
                int site_pos = (is_donor_site?gene.getExonStart(eidx)-1:gene.getExonEnd(eidx-1));
                countSite(ref, site_pos, false);
            }
        } else // gene on forward strand
        {
            for (int eidx=1; eidx<num_exons; eidx++)
            {
                int site_pos = (is_donor_site?gene.getExonEnd(eidx-1):gene.getExonStart(eidx)-1);
                countSite(ref, site_pos, true);
            }
        }
    }
    
    /**
     * Collects statistics about motifs in coding sequences (i.e., *not* splice sites) for 
     * null distribution.
     * 
     * @param gene
     * @param include_introns whether intron sequences are to be included (probably not)
     */
    public void countNonSites(GenePred gene, boolean include_introns)
    {
        String dna_seq;
        if (include_introns)
        {
            int cds_start = gene.getCDSStart();
            int cds_end = gene.getCDSEnd();

            dna_seq = gene.getChrom().getSubstring(cds_start, cds_end);
            if (gene.isReverseStrand() == is_donor_site)
            {
                dna_seq = DNASequence.reverseComplement(dna_seq);
            }
        } else
        {
            dna_seq = gene.getSequenceCDS();
            if (!is_donor_site)
                dna_seq = DNASequence.reverseComplement(dna_seq);
        }

        int max_pos = dna_seq.length()-motif_intronic_length-motif_boundary_length;
        for (int pos=0; pos<max_pos; pos++)
            countSite(dna_seq, pos);
        
//        dna_seq = DNASequence.reverseComplement(dna_seq);
//        for (int pos=0; pos<max_pos; pos++)
//            countSite(dna_seq, pos);
    }

    /**
     * 
     * 
     * @param seq will not be reverse complemented for acceptor site!
     * @param pos 
     */
    private void countSite(String seq, int pos)
    {
        int motif_idx = encodeMotif(seq, pos, true);
        motif_count.increment(motif_idx);
        
        int bpos = pos+motif_intronic_length;
        assert (bpos+motif_boundary_length<=seq.length());

        for (int bidx=0; bidx<motif_boundary_length; bidx++, bpos++)
        {
            char nuc = seq.charAt(bpos);
            if (!DNASequence.isAmbiguous(nuc))
                boundary_count[bidx][DNASequence.toIndex(nuc)]++;
        }
    }
    
    private void countSite(DNASequence seq, int pos, boolean is_forward_strand)
    {
        int motif_idx = encodeMotif(seq, pos, is_forward_strand == is_donor_site);
        motif_count.increment(motif_idx);
    
        int bpos = (is_forward_strand==is_donor_site?pos+motif_intronic_length:pos-motif_intronic_length);
        
        for (int bidx=0; bidx<motif_boundary_length; bidx++)
        {
            byte nuc_fw = seq.getNucleotideAt(bpos);
            if (!DNASequence.isAmbiguous(nuc_fw))
            {
                int nuc_idx = DNASequence.toIndex(is_forward_strand==is_donor_site?nuc_fw:DNASequence.complement(nuc_fw));
                boundary_count[bidx][nuc_idx]++;
            }
            if (is_forward_strand==is_donor_site)
            {
                bpos++;
            } else
            {
                bpos--;
            }
        }
    }
    
    public int encodeMotif(String seq, int first_pos, boolean is_forward)
    {
        int mot = 0;
        int pos = first_pos;
        
        if (is_forward)
        {
            if (pos+motif_intronic_length<=seq.length())
            {
                for (int offset=0;offset<motif_intronic_length; offset++, pos++)
                {
                    char nuc = seq.charAt(pos);
                    if (DNASequence.isAmbiguous(nuc))
                    {
                        mot = -1; 
                        break;
                    } else
                    {
                        mot = (mot << 2)|DNASequence.toIndex(nuc);
                    }
                }
            } else
                mot = -1;
        } else 
        {
            if (pos>=motif_intronic_length-1)
            {
                for (int offset=0; offset<motif_intronic_length ; offset++)
                {
                    char nuc_fw = seq.charAt(pos);
                    if (DNASequence.isAmbiguous(nuc_fw))
                    {
                        mot = -1;
                        break;
                    } else
                    {
                        mot = (mot << 2) | DNASequence.toIndex(DNASequence.complement(nuc_fw));
                    }
                    pos--;
                }
            }
            else
                mot = -1;
        }
        return mot;
    }
    
    
    /**
     * Sequence motif at a given position on forward or reverse strand. 
     * 
     * @param seq
     * @param first_pos
     * @param is_forward_strand
     * @return -1 if too close to sequence boundary or if an ambiguous nucleotide is encountered; otherwise a 2*length-bit integer
     */
    public int encodeMotif(DNASequence seq, int first_pos, boolean is_forward_strand)
    {
        int mot = 0;
        int pos = first_pos;
        if (is_forward_strand)
        {
            if (pos+motif_intronic_length<=seq.getLength()) // or else does not fit
                for (int offset=0; offset<motif_intronic_length ; offset++, pos++)
                {
                    byte nuc_fw = seq.getNucleotideAt(pos);
                    if (DNASequence.isAmbiguous(nuc_fw))
                    {
                        mot = -1;
                        break;
                    } else
                    {
                        mot = (mot << 2) | DNASequence.toIndex(nuc_fw);
                    }
                }
            else
                mot = -1;
        } else // reverse strand
        {
            if (pos>=motif_intronic_length-1)
                for (int offset=0; offset<motif_intronic_length ; offset++)
                {
                    byte nuc_fw = seq.getNucleotideAt(pos);
                    if (DNASequence.isAmbiguous(nuc_fw))
                    {
                        mot = -1;
                        break;
                    } else
                    {
                        mot = (mot << 2) | DNASequence.toIndex(DNASequence.complement(nuc_fw));
                    }
                    pos--;
                }
            else
                mot = -1;
        }
        return mot;
    }
        
    public String decodeMotif(int motif_idx)
    {
        if (motif_idx<0)
            return null;
        
        char[] motc = new char[motif_intronic_length];
        int cidx=motc.length;
        do
        {
            --cidx;
            int nuc_idx = motif_idx & 3;
            motc[cidx] = DNASequence.toChar(nuc_idx);
            motif_idx = motif_idx >>> 2;
        } while (cidx>0);
        return new String(motc);
    }
    
    public void reportStatistics(java.io.PrintStream out)
    {
        out.println("#SSC motif length "+motif_intronic_length+"\tboundary "+motif_boundary_length+"\t"
        +(is_donor_site?"5'/donor":"3'/acceptor"));
        
        Integer[] sorted_motifs = motif_count.toArray(new Integer[0]);
        
        if (is_donor_site)
        {
            for (int mot_idx: sorted_motifs)
            {
                String mot = decodeMotif(mot_idx);
                out.println("#SSC motif ["+(mot==null?"null":mot.toUpperCase())+"\t"+motif_count.getCount(mot_idx));
            }
            {
                int bidx=0; 
                while (bidx<motif_boundary_length)
                {
                    int offset = bidx+motif_intronic_length;
                    out.print("#SSC bndry "+offset);
                    double tot = 0.0;
                    for (int nuc_idx=0; nuc_idx<4; nuc_idx++)
                    {
                        int cnt = boundary_count[bidx][nuc_idx];
                        out.print("\t"+DNASequence.toChar(nuc_idx)+" "+cnt);
                        tot += cnt;
                    }
                    for (int nuc_idx=0; nuc_idx<4; nuc_idx++)
                    {
                        double p = boundary_count[bidx][nuc_idx]/tot;
                        out.printf("\t%.3f", p);
                    }
                    out.println();
                    bidx++;
                }
            }
            
        } else
        {
            {
                int bidx=motif_boundary_length;
                do
                {
                    --bidx;
                    int offset = -bidx-motif_intronic_length;
                    out.print("#SSC bndry "+offset);
                    double tot = 0.0;
                    for (int nuc_idx=0; nuc_idx<4; nuc_idx++)
                    {
                        int cnt = boundary_count[bidx][nuc_idx];
                        out.print("\t"+DNASequence.complement(DNASequence.toChar(nuc_idx))+" "+cnt);
                        tot += cnt;
                    }
                    for (int nuc_idx=0; nuc_idx<4; nuc_idx++)
                    {
                        double p = boundary_count[bidx][nuc_idx]/tot;
                        out.printf("\t%.3f", p);
                    }
                    
                    out.println();
                } while (bidx>0);
            }
            for (int mot_idx: sorted_motifs)
            {
                String mot = decodeMotif(mot_idx);
                out.println("#SSC motif "+(mot==null?"null":DNASequence.reverseComplement(mot).toUpperCase())+"]\t"+motif_count.getCount(mot_idx));
            }
        }
    }
    
    public void writeData(java.io.PrintStream out)
    {
        out.println("#SPLICE");
        out.println(motif_intronic_length+"\t"+motif_boundary_length+"\t"+is_donor_site
            +"\t# motif-length boundary-length is_donor");
        Integer[] motif_indexes = motif_count.toArray(new Integer[0]);
        out.println(motif_indexes.length+"\t# number of motifs");
        for (int i=0; i<motif_indexes.length; i++)
        {
            int mot_idx = motif_indexes[i];
            out.println(mot_idx+"\t"+motif_count.getCount(mot_idx)+"\t# motif "+decodeMotif(mot_idx));
        }
        String nucleotide_order_info = 
                "("+DNASequence.toChar(0)+DNASequence.toChar(1)+DNASequence.toChar(2)+DNASequence.toChar(3)+")";
        for (int bidx=0; bidx<boundary_count.length; bidx++)
        {
            for (int nuc_idx=0; nuc_idx<4; nuc_idx++)
            {
                out.print(boundary_count[bidx][nuc_idx]);
                out.print('\t');
            }
            out.println("# boundary "+bidx+" "+nucleotide_order_info);
        }
    }
    
            
    private static String[] skipCommentLines(java.io.BufferedReader in) throws java.io.IOException
    {
        String line;
        do
        {
            line = in.readLine();
        } while (line != null && line.startsWith("#"));
        return (line==null?null:line.split("\\t"));
    }
    
    public static SpliceSiteComposition readData(java.io.BufferedReader in) throws java.io.IOException
    {
        String[] fields = skipCommentLines(in);
        int mlen = Integer.parseInt(fields[0]);
        int blen = Integer.parseInt(fields[1]);
        boolean isdon = fields[2].equals("true");
        
        SpliceSiteComposition ss = new SpliceSiteComposition(mlen, blen, isdon);
        fields = skipCommentLines(in);
        int nmot = Integer.parseInt(fields[0]);
        for (int j=0; j<nmot; j++)
        {
            fields = skipCommentLines(in);
            int mot_idx = Integer.parseInt(fields[0]);
            int cnt = Integer.parseInt(fields[1]);
            ss.motif_count.set(mot_idx, cnt);
        }
        for (int bidx=0; bidx<blen; bidx++)
        {
            fields = skipCommentLines(in);
            for (int nuc_idx=0; nuc_idx<4; nuc_idx++)
            {
                ss.boundary_count[bidx][nuc_idx] = Integer.parseInt(fields[nuc_idx]);
            }
        }
        return ss;
    }

    //
    // 5ss-motif 5ss-boundary 3ss-motif 3ss-boundary: 5 4 4 8 
    public static void main(String[] args) throws Exception
    {
        if (args.length==0 || args[0].equals("-h")){
            System.err.println("Tests: call as $0 org1=g1.fa,org2=g2.fa annot.txt 5ss-motif 5ss-boundary 3ss-motif 3ss-boundary");
            System.err.println("Outputs splice site compositions.");
            System.exit(9);
        }
        
        int arg_idx=0;
        String genome_fa_list = args[arg_idx++];
        String annotation_file = args[arg_idx++];
        int donor_motif_length = Integer.parseInt(args[arg_idx++]);
        int donor_boundary_length = Integer.parseInt(args[arg_idx++]);
        int acceptor_motif_length = Integer.parseInt(args[arg_idx++]);
        int acceptor_boundary_length = Integer.parseInt(args[arg_idx++]);

        AnnotatedGenomes annotations = new AnnotatedGenomes();
        List<String> wanted_organisms = annotations.readMultipleGenomes(genome_fa_list);
        annotations.readAnnotations(annotation_file);    
        
        SpliceSiteComposition ss5[] = new SpliceSiteComposition[wanted_organisms.size()];
        SpliceSiteComposition ss3[] = new SpliceSiteComposition[wanted_organisms.size()];
        
        for (int org_idx=0; org_idx<wanted_organisms.size(); org_idx++)
        {
            String org = wanted_organisms.get(org_idx);
            System.out.println("#SSC organism "+org);
            
            SpliceSiteComposition don = ss5[org_idx] = new SpliceSiteComposition(donor_motif_length, donor_boundary_length, true);
            SpliceSiteComposition acc = ss3[org_idx] = new SpliceSiteComposition(acceptor_motif_length, acceptor_boundary_length, false);

            for (GenePred gene: annotations.getAllAnnotations(org))
            {
                don.countIntronSites(gene);
                acc.countIntronSites(gene);

            } // for gene 
            
            //don.reportStatistics(System.out);
            //acc.reportStatistics(System.out);
            don.writeData(System.out);
            acc.writeData(System.out);
            
        } // for org
    }
}

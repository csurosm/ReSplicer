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

import splice.model.SpliceSiteHomology;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import splice.annotation.AnnotatedGenomes;
import splice.annotation.GenePred;
import splice.empirical.Aggregator;
import splice.empirical.Counter;
import splice.sequence.ProteinSequence;

/**
 *
 * Counts how often introns are near gaps in alignments.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class countGaps 
{
    
    private void initDataStructures(List<String> wanted_organisms)
    {
        int norg = wanted_organisms.size();
        this.organism_index = new HashMap<>();
        pairwise_counters = new GapCounter[norg][norg];
        for (int i=0; i<norg; i++)
        {
            String org =wanted_organisms.get(i);
            organism_index.put(org, i);
            for (int j=0; j<norg; j++)
                pairwise_counters[i][j]=new GapCounter();
        }
    }
    
    private GapCounter[][] pairwise_counters;
    private Map<String, Integer> organism_index;

    private static void printTally(java.io.PrintStream out, String prefix, Counter<Integer> cnt)
    {
        Integer[] lengths = cnt.keySet().toArray(new Integer[0]);
        Arrays.sort(lengths);
        for (int len: lengths)
        {
            out.println(prefix+"\t"+len+"\t"+cnt.getCount(len));
        }
        out.println(prefix+"\ttotal\t"+cnt.totalCounts());
    }
    
    private static class GapCounter
    {
        final Counter<Integer> insert_length;
        final Counter<Integer> delete_length;
        final Counter<Integer> match_length;
        
        private enum State {MATCH, INSERT, DELETE};
        
        State current_state;
        int segment_length;
        
        final List<Set<Integer>> segment_introns;
        final List<Counter<Integer>> upstream_gap_distance; 
        final List<Counter<Integer>> downstream_gap_distance;
        
        /**
         * Class for collecting statistics about 
         * gap lengths, and the intron placements relative to the gaps.
         * 
         * Phase-1 and -2 introns in the i-th codon of a length-n run of 
         * match columns counted for i upstream distance, and n-i+1 downstream distance.
         * Phase-1 and -2 introns within indels counted for distances 0.
         * Phase-3 intron after the i-th codon in a length-n match segment 
         * is counted for i upstream distance and n-i downstream distance. 
         * Phase-3 intron preceded by indel is counted for 0 upstream; phase-3
         * succeded by indel is counted for 0 downstream.
         */
        GapCounter()
        {
            insert_length = new Counter<>();
            delete_length = new Counter<>();
            match_length = new Counter<>();
            
            this.segment_introns = new ArrayList<>();
            this.upstream_gap_distance = new ArrayList<>();
            this.downstream_gap_distance = new ArrayList<>();
            for (int i=0; i<6; i++) 
            {
                segment_introns.add(new HashSet<>());
                upstream_gap_distance.add(new Counter<>());
                downstream_gap_distance.add(new Counter<>());
            }

            current_state = State.MATCH;
            segment_length=0;
        }
        
        /**
         * Reaching a gap: all introns in match segment are now counted 
         * for downstream distances.
         */
        private void storeDownstreamGapsAll()
        {
            for (int i=0; i<6; i++)
            {
                for (int offset: segment_introns.get(i))
                {
                    downstream_gap_distance.get(i).increment(segment_length-offset);
                }
                segment_introns.get(i).clear();
            }
        }
        
        /**
         * Phase-3 intron in preceding column is counted for downstream distance 0.
         * 
         * @param sidx 1 or 2
         */
        private void storeDownstreamGaps0(int sidx)
        {
            int iph = getIntronIndex(sidx,3);
            if (!segment_introns.get(iph).isEmpty())
            {
                downstream_gap_distance.get(iph).increment(0);
                segment_introns.get(iph).clear();
            }
        }
        
        /**
         * Indexing in {@link #downstream_gap_distance} and {@link #segment_introns} lists.
         * 
         * @param seq_idx 1 or 2 
         * @param phase 1, 2 or 3
         * @return 0..5
         */
        private int getIntronIndex(int seq_idx, int phase)
        {
            return 3*seq_idx-phase;
        }
        
        /**
         * Counting downstream distances are the end of the sequence.
         */
        void endAlignment()
        {
            if (current_state==State.MATCH)
            {
                match_length.increment(segment_length);
                storeDownstreamGapsAll();
            } else if (current_state == State.DELETE)
            {
                delete_length.increment(segment_length);
                // phase-3 intron in last position is ignored
            } else if (current_state==State.INSERT)
            {
                insert_length.increment(segment_length);
                // phase-3 intron in last position is ignored
            }
            segment_length = 0;
        }
        
        void reportStatistics(java.io.PrintStream out)
        {
            printTally(out, "#INSERT", insert_length);
            printTally(out, "#DELETE", delete_length);
            printTally(out, "#MATCH", match_length);
            
            for (int sidx=1; sidx<=2; sidx++)
            {
                for (int ph=1; ph<=3; ph++)
                {
                    int idx = getIntronIndex(sidx, ph);
                    printTally(out, "#UP  \t"+sidx+".ph"+ph, upstream_gap_distance.get(idx));
                    printTally(out, "#DOWN\t"+sidx+".ph"+ph, downstream_gap_distance.get(idx));
                }
            }
        }
        
        void addPair(boolean is_gap1, int intron_phase1, boolean is_gap2, int intron_phase2)
        {
            if (is_gap1)
            {
                if (is_gap2)
                {
                    // nothing to do
                } else // insert
                {
                    if (current_state==State.MATCH) // M->I
                    {
                        if(segment_length>0) // except for the very first call right after initialization with segment_length=0
                        {
                            match_length.increment(segment_length);
                            storeDownstreamGapsAll();
                        }
                        segment_length=0;
                    } else if (current_state == State.DELETE) // D-> I
                    {
                        delete_length.increment(segment_length);
                        storeDownstreamGaps0(1);
                        segment_length = 0;
                    } else if (current_state == State.INSERT) // I->I
                    {
                        storeDownstreamGaps0(2);
                    }
                    segment_length++;
                    if (intron_phase2 != SpliceSiteHomology.NO_INTRON)
                    {
                        int iph = getIntronIndex(2, intron_phase2);
                        upstream_gap_distance.get(iph).increment(0);
                        if (intron_phase2==1 || intron_phase2==2)
                            downstream_gap_distance.get(iph).increment(0);
                        else if (intron_phase2==3) // wait to see next position
                            segment_introns.get(iph).add(0);
                    }
                    current_state = State.INSERT;
                }
            } else // no gap 1
                        {
                if (is_gap2) // delete
                {
                    if (current_state==State.MATCH) // M->D
                    {
                        if (segment_length>0)
                        {
                            match_length.increment(segment_length);
                            storeDownstreamGapsAll();
                        }
                        segment_length = 0;
                    } else if (current_state == State.DELETE) // D->D
                    {
                        storeDownstreamGaps0(1);
                    } else if (current_state == State.INSERT) // I->D
                    {
                        insert_length.increment(segment_length);
                        storeDownstreamGaps0(2);
                        segment_length = 0;
                    }
                    segment_length++;
                    
                    if (intron_phase1 != SpliceSiteHomology.NO_INTRON)
                    {
                        int iph = getIntronIndex(1, intron_phase1);
                        upstream_gap_distance.get(iph).increment(0);
                        if (intron_phase1==1 || intron_phase1==2)
                            downstream_gap_distance.get(iph).increment(0);
                        else if (intron_phase1==3)
                            segment_introns.get(iph).add(0);
                    }
                    this.current_state = State.DELETE;
                } else // match
                {
                    if (current_state == State.MATCH)
                    {
                        // M -> M
                    } else if (current_state == State.DELETE)
                    {
                        // D-> M
                        delete_length.increment(segment_length); // keep phase-3 intron at boundary in segment_introns
                        segment_length=0;
                    } else
                    {
                        insert_length.increment(segment_length);
                        segment_length=0;
                    }
                    segment_length++;
                    if (intron_phase1!=SpliceSiteHomology.NO_INTRON)
                    {
                        int iph = getIntronIndex(1, intron_phase1);
                        upstream_gap_distance.get(iph).increment(segment_length);
                        if (intron_phase1==3)
                            segment_introns.get(iph).add(segment_length);
                        else
                            segment_introns.get(iph).add(segment_length-1);
                    } 
                    if (intron_phase2!=SpliceSiteHomology.NO_INTRON)
                    {
                        int iph = getIntronIndex(2, intron_phase2);
                        upstream_gap_distance.get(iph).increment(segment_length);
                        if (intron_phase2==3)
                            segment_introns.get(iph).add(segment_length);
                        else
                            segment_introns.get(iph).add(segment_length-1);
                    }
                    current_state = State.MATCH;
                } // gap2
            } // gap2
        } // addPair()
    }
    
    
    private void reportStatistics(java.io.PrintStream out)
    {
        String[] wanted_organisms = organism_index.keySet().toArray(new String[0]);
        Arrays.sort(wanted_organisms, new Comparator<String>(){
            @Override
            public int compare(String o1, String o2) 
            {
                return Integer.compare(organism_index.get(o1), organism_index.get(o2));
            }
        });
        
        for (int i=0; i<wanted_organisms.length; i++)
        {
            Counter<Integer> ph_cnt = phase_counters[i].phases;
            int tot = ph_cnt.getCount(1)+ph_cnt.getCount(2)+ph_cnt.getCount(3);
            Counter<Integer> intr_cnt =
                phase_counters[i].introns_per_gene;
            int ngenes = intr_cnt.totalCounts();
            int tot_intr=0;
            for (int j: intr_cnt.keySet())
                tot_intr += j*intr_cnt.getCount(j);
            double avg_introns_per_gene = tot_intr / ((double)ngenes);
            out.println("#TALLY\t"+wanted_organisms[i]
                    +"\t"+ph_cnt.getCount(1)
                    +"\t"+ph_cnt.getCount(2)
                    +"\t"+ph_cnt.getCount(3)
                    +"\t"+tot
                    +"\t// ngenes "+ngenes+", totintr "+tot_intr+", avgpergene "+avg_introns_per_gene);
            
        }
        
        
        
        for (int i=0; i<wanted_organisms.length; i++)
            for (int j=0; j<i; j++)
            {
                out.println("#STATISTICS\t"+wanted_organisms[i]+"/"+wanted_organisms[j]);
                pairwise_counters[i][j].reportStatistics(out);
            }
    }
    
    private class PhaseCounter
    {
        private Counter<Integer> phases = new Counter<>();
        private Counter<Integer> introns_per_gene = new Counter<>();
    }
    
    private PhaseCounter[] phase_counters = null;
    
    private void processAlignment(SpliceSiteHomology aligned_sites)
    {
        if (phase_counters == null)
        {
            phase_counters = new PhaseCounter[aligned_sites.getNumSequences()];
            for (int ali_i=0; ali_i<aligned_sites.getNumSequences(); ali_i++)
            {
                phase_counters[ali_i] = new PhaseCounter();
            }
        }
        for (int ali_i=0; ali_i<aligned_sites.getNumSequences(); ali_i++)
        {
            GenePred gene_i = aligned_sites.getGene(ali_i);
            String org_i = aligned_sites.getOrganism(ali_i);
            int i = organism_index.get(org_i);
            phase_counters[i].introns_per_gene.increment(gene_i.getExonCount()-1);
        }
        
        for (int col=0; col<aligned_sites.getNumAlignmentColumns(); col++)
        {
            for (int ali_i=0; ali_i<aligned_sites.getNumSequences(); ali_i++)
            {
                GenePred gene_i = aligned_sites.getGene(ali_i);
                String org_i = aligned_sites.getOrganism(ali_i);
                int i = organism_index.get(org_i);

                boolean is_gap_i = aligned_sites.getAlignedGenomicPosition(ali_i, col)==SpliceSiteHomology.ALIGNMENT_GAP;
                int phase_i = aligned_sites.getIntronPhase(gene_i, col);
                phase_counters[i].phases.increment(phase_i);

                for (int ali_j=0; ali_j<aligned_sites.getNumSequences(); ali_j++)
                {
                    GenePred gene_j = aligned_sites.getGene(ali_j);
                    String org_j = aligned_sites.getOrganism(ali_j);
                    int j = organism_index.get(org_j);
                    if (i>j)
                    {

                        boolean is_gap_j = aligned_sites.getAlignedGenomicPosition(ali_j, col)==SpliceSiteHomology.ALIGNMENT_GAP;
                        int phase_j = aligned_sites.getIntronPhase( gene_j, col);

                        GapCounter cnt = pairwise_counters[i][j];
                        cnt.addPair(is_gap_i, phase_i, is_gap_j, phase_j);
                    }
                }
            }
        }

        for (int i=0; i<aligned_sites.getNumSequences(); i++)
        {
            for(int j=0; j<aligned_sites.getNumSequences(); j++)
                if (i!=j)
                {
                    GapCounter cnt = pairwise_counters[i][j];
                    cnt.endAlignment();
                }
        }        
    }
    
    
    private void huzzah(String[] args) throws Exception
    {
        if (args.length==0 || args[0].equals("-h")){
            System.err.println("Call as "+getClass().getName()+" org1=g1.fa,org2=g2.fa annot.txt agg.data ali1.fa ali2.fa ...");
            System.err.println("Calculates gap statistics.");
            System.exit(9);
        } 
        
        int arg_idx=0;
        String genome_fa_list = args[arg_idx++];
        String annotation_file = args[arg_idx++];
        String scoring_file = args[arg_idx++];
        
        AnnotatedGenomes annotations = new AnnotatedGenomes();
        List<String> wanted_organisms = annotations.readMultipleGenomes(genome_fa_list);
        annotations.readAnnotations(annotation_file); 
        
        Aggregator all_scoring = Aggregator.readData(new GeneralizedFileReader(scoring_file));
        this.initDataStructures(wanted_organisms);
    
        int nproc = 0;
        while (arg_idx<args.length)
        {
            String alignment_file = args[arg_idx];
            ProteinSequence[] alignment= ProteinSequence.readMultiFasta(alignment_file);
            SpliceSiteHomology aligned_sites = new SpliceSiteHomology(all_scoring);
            if (aligned_sites.computeAlignmentColumns(annotations, alignment))        
            {
//                aligned_sites.printAlignment(System.out);
                processAlignment(aligned_sites);
                nproc++;
            } else
            {
                System.out.println("#SKIP "+alignment_file);
            }
            arg_idx++;
        }
        
        reportStatistics(System.out);
    }
    
    public static void main(String[] args) throws Exception
    {
        (new countGaps()).huzzah(args);
        
    }
    
}

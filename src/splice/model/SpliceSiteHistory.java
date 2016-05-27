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

package splice.model;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import splice.sequence.DNASequence;
import splice.sequence.DNAStrand;
import splice.sequence.MultiHomology;

/**
 *
 * @author Mikl&oacute;don Cs&#369;r&ouml;don
 */
public class SpliceSiteHistory 
{
    private final MultiHomology donor_homologies;
    private final MultiHomology acceptor_homologies;
    
    private final Map<DNAStrand, Map<SpliceSite, Boolean>> putative_sites;
    private final Map<String, DNAStrand> gene_strands;
    
    public SpliceSiteHistory(MultiHomology donor_homologies, MultiHomology acceptor_homologies)
    {
        this.donor_homologies = donor_homologies;
        this.acceptor_homologies = acceptor_homologies;
        this.putative_sites = new HashMap<>();
        for (DNAStrand seq: donor_homologies.getSequences())
        {
            putative_sites.put(seq, new HashMap<>());
        }
        SpliceSite no_site = new SpliceSite(SiteType.NOT_A_SITE, 0, -1);
        SpliceSite missing = new SpliceSite(SiteType.NOT_A_SITE, -1, -1);

        NO_INTRON = new IntronPlacement(no_site, no_site);
        MISSING_VALUE = new IntronPlacement(missing, missing);
        
        gene_strands = new HashMap<>();
        for (DNAStrand seq: donor_homologies.getSequences())
        {
            String org = seq.getDNA().getOrganism();
            gene_strands.put(org, seq);
        }
    }

    public DNAStrand getSequence(String organism)
    {
        return gene_strands.get(organism);
    }
    
    private final IntronPlacement NO_INTRON;
    private final IntronPlacement MISSING_VALUE;
    
    public void addDonorSite(DNAStrand seq, int pos, int phase, boolean is_annotated)
    {
        int mapped_pos = donor_homologies.getAbsolutePosition(seq, pos);
        SpliceSite donor = new SpliceSite(SiteType.DONOR, mapped_pos, phase);
        Map<SpliceSite, Boolean> known_sites = putative_sites.get(seq);
        boolean was_annotated = known_sites.getOrDefault(donor, false);
        known_sites.put(donor, was_annotated || is_annotated);
    }
    
    public void addAcceptorSite(DNAStrand seq, int pos, int phase, boolean is_annotated)
    {
        int mapped_pos = acceptor_homologies.getAbsolutePosition(seq, pos);
        SpliceSite acceptor = new SpliceSite(SiteType.ACCEPTOR, mapped_pos, phase);
        Map<SpliceSite, Boolean> known_sites = putative_sites.get(seq);
        boolean was_annotated = known_sites.getOrDefault(acceptor, false);
        known_sites.put(acceptor, was_annotated || is_annotated);
    }
    
    public void addNoIntron(DNAStrand seq, boolean is_annotated)
    {
        SpliceSite nointron = new SpliceSite(SiteType.NOT_A_SITE, 0, -1);
        Map<SpliceSite, Boolean> known_sites = putative_sites.get(seq);
        boolean was_annotated = known_sites.getOrDefault(nointron, false);
        known_sites.put(nointron, was_annotated || is_annotated);
    }
    
    public IntronPlacement getAnnotation(DNAStrand seq)
    {
        if (getAnnotatedSite(seq, SiteType.NOT_A_SITE)==null)
        {
            SpliceSite donor_site = getAnnotatedSite(seq, SiteType.DONOR);
            SpliceSite acceptor_site = getAnnotatedSite(seq, SiteType.ACCEPTOR);
            if (donor_site==null || acceptor_site == null)
                return null;
            else
                return new IntronPlacement(donor_site, acceptor_site);
        } else
        {
            return NO_INTRON;
        }
    }
    
    /**
     * The annotated donor/ acceptor / nointron site for a sequence.
     * 
     * @param seq
     * @param wanted_type
     * @return null if no annotated site of given type
     */
    private SpliceSite getAnnotatedSite(DNAStrand seq, SiteType wanted_type)
    {
        SpliceSite annotated_site = null;
        Map<SpliceSite, Boolean> known_sites = putative_sites.get(seq);
        if (known_sites!=null)
        {
            for (SpliceSite S: known_sites.keySet())
            {
                if (known_sites.get(S))
                {
                    if (S.type==wanted_type)
                    {
                        annotated_site=S;
                        break;
                    }
                }
            }
        }
        return annotated_site;
    }
    
    private String getSiteDescription(IntronPlacement inferred, DNAStrand seq)
    {
        Map<SpliceSite, Boolean> known_sites = putative_sites.get(seq);
        StringBuilder sb = new StringBuilder();
        
        if (known_sites==null)
        {
            assert (!inferred.hasIntron());
            sb.append("{no sites}");
        } else
        {
            SpliceSite annotated_donor = null;
            SpliceSite annotated_acceptor = null;
            for (SpliceSite S: known_sites.keySet())
            {
                if (known_sites.get(S))
                {
                    if (S.type==SiteType.DONOR)
                        annotated_donor=S;
                    else
                        annotated_acceptor=S;
                }
            }
            if (annotated_donor==null)
            {
                if (!inferred.hasIntron())
                {
                    // no annotated sites, no inferred sites
                } else
                {
                    // no annotated sites but inferred intron presence
                    sb.append("ss5(+").append(inferred.donor_site.getDinucleotideMotif(seq)).append(")");
                    sb.append(",");
                    sb.append("ss3(+").append(inferred.acceptor_site.getDinucleotideMotif(seq)).append(")");
                }
            } else
            {
                if (!inferred.hasIntron())
                {
                    // annotated site but inferred intron absence
                    sb.append("ss5(-),ss3(-)");
                } else
                {
                    int donor_shift = inferred.donor_site.getOffsetFromAnnotation(seq);
                    int acceptor_shift = inferred.acceptor_site.getOffsetFromAnnotation(seq);
                    sb.append("ss5(").append(donor_shift).append("/").append(inferred.donor_site.getDinucleotideMotif(seq)).append(")");
                    sb.append(",");
                    sb.append("ss3(").append(acceptor_shift).append("/").append(inferred.acceptor_site.getDinucleotideMotif(seq)).append(")");
                }
            }
            sb.append(" {");
            for (SpliceSite S: known_sites.keySet())
            {
                sb.append(S).append("/").append(S.getDinucleotideMotif(seq)).append(known_sites.get(S)?"+":"-");
                sb.append(" ");
            }
            sb.append("}");
        }
        return sb.toString();
    }
    

    private Set<IntronPlacement> getAllPlacements()
    {
        Set<IntronPlacement> all_states = new HashSet<>();

        all_states.add(NO_INTRON);
        all_states.add(MISSING_VALUE);
        
        for (DNAStrand seq: putative_sites.keySet())
        {
            Map<SpliceSite, Boolean> known_sites = putative_sites.get(seq);
            for (SpliceSite don: known_sites.keySet())
            {
                if (don.type == SiteType.DONOR)
                {
                    for (SpliceSite acc: known_sites.keySet())
                        if (acc.type==SiteType.ACCEPTOR && acc.phase==don.phase)
                        {
                            all_states.add(new IntronPlacement(don, acc));
                        }
                    }
                }
        }
        return all_states;
    }
    
    private static final int MIN_INTRON_LENGTH = 24;
    
    
    public void reportHistory(java.io.PrintStream out, TreeNode root, Map<TreeNode, IntronPlacement> node_labeling)
    {
        TreeNode[] nodes = root.getTraversal().getDFT();
        for (TreeNode N: nodes)
        {
            IntronPlacement node_state = node_labeling.get(N);
            IntronPlacement parent_state = N.isRoot()?node_state:node_labeling.get(N.getParent());
            
            String ss_motifs="";
            if (N.isLeaf() && node_state.hasIntron())
            {
                String donor_site_motif = node_state.donor_site.getDinucleotideMotif(gene_strands.get(N.getName()));
                String acceptor_site_motif = node_state.acceptor_site.getDinucleotideMotif(gene_strands.get(N.getName()));
                ss_motifs = "\tss:"+donor_site_motif+"-"+acceptor_site_motif;
            }
            
            out.println("#HISTORY "
                    +(N.isLeaf()?"leaf ":"node ")
                    +N.newickName()
                    +"\t"+(N.isLeaf()?node_state.mapToOrganism(gene_strands.get(N.getName())):node_state)
                    +"\t"+node_state.getEventFrom(parent_state)
                    +"\t"+(node_state.equals(parent_state)?"==":parent_state.toString())
                    +(N.isRoot()?"":"\t".concat(N.getParent().newickName()))
                    +(N.isLeaf()?"\t"+getSiteDescription(node_state, gene_strands.get(N.getName())):"")
                    +ss_motifs
            );
        }
    }
    
    
    
    /**
     * Tab-delimited description of events along the tree.
     * @param root phylogeny root
     * @param node_labeling if null, a header line is constructed; otherwise return value of {@link #ancestralParsimony(splice.model.TreeNode, int, int, int, int) } is expected here
     * @return 
     */
    public String getAncestralDescription(TreeNode root, Map<TreeNode, IntronPlacement> node_labeling)
    {
        TreeNode[] nodes = root.getTraversal().getDFT();
        StringBuilder description = new StringBuilder();
        if (node_labeling==null)
        {
            for (TreeNode N: nodes)
            {
                if (description.length()>0)
                    description.append("\t");
                description.append(N.newickName());
            }
        } else
        {
            int terminal_intron_count = 0;
            for (TreeNode N: nodes)
            {
                IntronPlacement node_state = node_labeling.get(N);
                IntronPlacement parent_state = N.isRoot()?node_state:node_labeling.get(N.getParent());
                if (description.length()>0)
                    description.append("\t");
                if (node_state.hasIntron())
                    description.append(node_state.donor_site.phase).append(':');
                description.append(node_state.getEventFrom(parent_state));

                if (N.isLeaf() && node_state.hasIntron())
                    terminal_intron_count++;
            }
            description.append("\t").append(terminal_intron_count);
        }
        return description.toString();        
    }
    
    /**
     * Ancestral reconstruction on a tree with given penalties.
     * 
     * @param root root of the evolutionary tree; leaves have {@link TreeNode#getName() } equal to origanism names
     * @param loss_pty penalty for complete intron loss
     * @param gain_pty penalty for intron gain
     * @param shift_pty penalty for shifted splice sites
     * @param reannotate_pty penalty for reannotating splice sites
     * @return map from treenodes to optimal-scoring parsimony reconstruction
     */
    public Map<TreeNode, IntronPlacement> ancestralParsimony(TreeNode root, int loss_pty, int gain_pty, int shift_pty, int reannotate_pty)
    {
        TreeNode[] nodes = root.getTraversal().getDFT();

        Map<TreeNode, Map<IntronPlacement, Integer>> subtree_scores = new HashMap<>();
        Map<TreeNode, Map<IntronPlacement, IntronPlacement>> score_backtracking = new HashMap<>();
        for (TreeNode N: nodes) 
            if (!N.isRoot())
                score_backtracking.put(N, new HashMap<>());
        
        Set<IntronPlacement> all_states = getAllPlacements();
        
        int very_large_pty = nodes.length*(loss_pty+gain_pty+shift_pty+reannotate_pty);
        
        for (int node_idx=0; node_idx<nodes.length; node_idx++)
        {
            TreeNode N = nodes[node_idx];
            Map<IntronPlacement, Integer> ancestral_states = new HashMap<>();
            if (N.isLeaf())
            {
                String org = N.getName();
                DNAStrand seq = gene_strands.get(org);
                
                if (putative_sites.containsKey(seq))
                {
                    Map<SpliceSite, Boolean> known_sites = putative_sites.get(seq);

                    for (SpliceSite don: known_sites.keySet())
                    {
                        if (don.type==SiteType.NOT_A_SITE)
                        {
                            int pty = (known_sites.get(don)?0:reannotate_pty);
                            ancestral_states.put(NO_INTRON, pty);
                        } else if (don.type == SiteType.DONOR)
                        {
                            int dpty = known_sites.get(don)?0:reannotate_pty;
                            for (SpliceSite acc: known_sites.keySet())
                                if (acc.type==SiteType.ACCEPTOR && acc.phase==don.phase)
                                {
                                    int apty = known_sites.get(acc)?0:reannotate_pty;
                                    IntronPlacement ichabod = new IntronPlacement(don, acc);
                                    // check if length is ok
                                    IntronPlacement mapped_here = ichabod.mapToOrganism(seq);

                                    if ((known_sites.get(don) && known_sites.get(acc)) || seq.getDownstreamStep()*(mapped_here.acceptor_site.pos-mapped_here.donor_site.pos)>=MIN_INTRON_LENGTH-1)
                                        ancestral_states.put(ichabod, dpty+apty);
                                }
                        }
                    }
                } else // unresolved sequence
                {
                    ancestral_states.put(MISSING_VALUE, 0);
                }
                        
//                System.out.print("#*SSH.aP leaf "+N+"\torg "+org+"\tseq "+seq+"\t"+getSiteDescription(seq));
//                for (IntronPlacement ichabod: ancestral_states.keySet())
//                {
//                    System.out.print("\t"+ichabod+"/"+ancestral_states.get(ichabod));
//                }
//                System.out.println();
            } else
            {
                // inner node
                for (IntronPlacement parent_state: all_states)
                {
                    int pty = 0;
                    for (int cidx=0; cidx<N.getNumChildren(); cidx++)
                    {
                        TreeNode C = N.getChild(cidx);
                        Map<IntronPlacement, Integer> child_scores = subtree_scores.get(C);
                        Map<IntronPlacement, IntronPlacement> child_backtrack = score_backtracking.get(C);

                        IntronPlacement best_child_placement =  null;
                        int min_child_score = -1;
                        for (IntronPlacement child_state: child_scores.keySet())
                        {
                            int subtree_pty = child_scores.get(child_state);
                            int edge_pty;
                            
                            if (parent_state == MISSING_VALUE)
                                edge_pty = (child_state==MISSING_VALUE?0:very_large_pty);
                            else
                            {
                                if (child_state == NO_INTRON)
                                {
                                    edge_pty = (parent_state==NO_INTRON?0:loss_pty);
                                } else if (child_state == MISSING_VALUE)
                                {
                                    edge_pty = 0;
                                } else
                                {
                                    if (parent_state==NO_INTRON)
                                    {
                                        edge_pty = gain_pty;
                                    } else
                                    {
                                        int dpty = (child_state.donor_site.pos == parent_state.donor_site.pos?0:shift_pty);
                                        int apty = (child_state.acceptor_site.pos == parent_state.acceptor_site.pos?0:shift_pty);
                                        edge_pty = apty+dpty;
                                    }
                                }
                            }
                            int sc = subtree_pty + edge_pty;
                            if (best_child_placement == null 
                                    || sc<min_child_score
                                    || (sc==min_child_score 
                                        && child_state.equals(parent_state))) // implements DELTRAN: state transitions are preferred farther from the root
                            {
                                best_child_placement = child_state;
                                min_child_score = sc;
                            }
                        } // for all child states
                        child_backtrack.put(parent_state, best_child_placement);
                        pty += min_child_score;
                    } // for each child
                    ancestral_states.put(parent_state, pty);
                    
                } // for all parental states
//                System.out.print("#*SSH.aP inner "+N);
//                for (IntronPlacement ichabod: ancestral_states.keySet())
//                {
//                    System.out.print("\t"+ichabod+"/"+ancestral_states.get(ichabod));
//                }
//                System.out.println();
            }
            subtree_scores.put(N, ancestral_states);
        } // for all nodes
        
        Map<TreeNode, IntronPlacement> node_labeling = new HashMap<>();
        for (int node_idx=nodes.length; node_idx>0; )
        {
            --node_idx;
            TreeNode N = nodes[node_idx];   
            if (N.isRoot())
            {
                Map<IntronPlacement, Integer> root_scores = subtree_scores.get(N);
                IntronPlacement best_root_label = null;
                int min_root_score = -1;
                for (IntronPlacement root_label: root_scores.keySet())
                {
                    int pty = root_scores.get(root_label);
                    if (best_root_label==null || min_root_score>pty)
                    {
                        best_root_label = root_label;
                        min_root_score = pty;
                    }
                }
                node_labeling.put(N, best_root_label);
            } else
            {
                TreeNode P = N.getParent();
                IntronPlacement parent_label = node_labeling.get(P);
                IntronPlacement child_label = score_backtracking.get(N).get(parent_label);
                node_labeling.put(N, child_label);
            }
        }
        
        return node_labeling;

    }
    
    private enum SiteType {DONOR, ACCEPTOR, NOT_A_SITE};
    public class SpliceSite
    {
        private final SiteType type;
        private final int pos;
        private final int phase;
        
        private SpliceSite(SiteType type, int pos, int phase)
        {
            this.pos = pos;
            this.type = type;
            this.phase = phase;
        }
        
        @Override
        public boolean equals(Object o)
        {
            if (o instanceof SpliceSite)
            {
                SpliceSite s = (SpliceSite)o;
                return this.type==s.type && this.pos==s.pos;
            } else
            {
                return super.equals(o);
            }
        }

        @Override
        public int hashCode() 
        {
            int hash = Objects.hashCode(this.type)*3+this.pos; 
            return hash;
        }
        
        @Override
        public String toString()
        {
            StringBuilder sb = new StringBuilder();
            sb.append(type.toString().charAt(0)).append(phase<0?"":Integer.toString(phase)).append(".").append(Integer.toString(pos));
            return sb.toString();
        }
        
        public int getPosition()
        {
            return pos;
        }
        
        public int getPhase()
        {
            return phase;
        }
        
        public SpliceSite mapToSequence(DNAStrand seq)
        {
            if (type==SiteType.DONOR)
            {
                int dpos = donor_homologies.projectAbsolutePosition(seq, pos);
                SpliceSite donor = new SpliceSite(type, dpos, phase);
                return donor;
            } else if (type==SiteType.ACCEPTOR)
            {
                int apos = acceptor_homologies.projectAbsolutePosition(seq, pos);
                SpliceSite acceptor = new SpliceSite(type, apos, phase);
                return acceptor;
            } else
                return this;
        }
        
        public String getDinucleotideMotif(DNAStrand seq)
        {
            if (type==SiteType.NOT_A_SITE)
                return "";
            SpliceSite mapped_site = mapToSequence(seq);
            if (mapped_site.pos<0 || mapped_site.pos>=seq.getDNA().getLength())
            {
                return "..";
            }
            byte nuc1 = seq.getNucleotideAt(mapped_site.pos);
            byte nuc2;
            if (type == SiteType.DONOR)
            {
                nuc2 = seq.getNucleotideAt(seq.add(mapped_site.pos, 1));
            } else
            {
                nuc2 = nuc1;
                nuc1 = seq.getNucleotideAt(seq.add(mapped_site.pos, -1));
            }
            char[] motif = new char[2];
            motif[0] = DNASequence.toChar(nuc1);
            motif[1] = DNASequence.toChar(nuc2);
            return new String(motif);
        }
        
        public int getOffsetFromAnnotation(DNAStrand seq)
        {
            SpliceSite annotated_site = getAnnotatedSite(seq, type);
            if (equals(annotated_site))
                return 0;
            
            if (type == SiteType.NOT_A_SITE)
            {
                return Integer.MIN_VALUE; // there was some annotation for this seq
            } else if (annotated_site==null)
            {
                return Integer.MAX_VALUE;    
            } else
            {
                SpliceSite mapped_this = mapToSequence(seq);
                SpliceSite mapped_that = annotated_site.mapToSequence(seq);
                int offset = seq.sub(mapped_this.pos, mapped_that.pos);
                return offset;
            }
        }
    }    
    
    private enum EventType
    {
        INTRON_GAIN("Gn"),
        INTRON_LOSS("Ls"),
        DONOR_INTRONIZATION("I5"),
        DONOR_EXONIZATION("E5"),
        ACCEPTOR_INTRONIZATION("I3"),
        ACCEPTOR_EXONIZATION("E3"),
        DONOR_CONSERVATION("C5"),
        ACCEPTOR_CONSERVATION("C3"),
        NO_INTRON("NO"),
        MISSING_DATA("na");
        
        EventType(String short_name)
        {
            this.short_name = short_name;
        }
        
        String getShortName(){ return short_name;}
        private final String short_name;
    }
    
    private static class Event
    {
        Event(EventType type)
        {
            this(type, 0, type, 0);
        }
        Event(EventType type, int amount)
        {
            this(type, amount, type, amount);
        }
        Event(EventType donor_event, int donor_amount, EventType acceptor_event, int acceptor_amount)
        {
            this.donor_event = donor_event;
            this.donor_amount = donor_amount;
            this.acceptor_event = acceptor_event;
            this.acceptor_amount= acceptor_amount;
        }

        private final int donor_amount;
        private final EventType donor_event;
        private final int acceptor_amount;
        private final EventType acceptor_event;
        
        @Override
        public String toString()
        {
            StringBuilder sb = new StringBuilder(donor_event.getShortName());
            if (donor_amount!=0)
                sb.append("(").append(donor_amount).append(")");
            if (donor_event != acceptor_event)
            {
                sb.append(",").append(acceptor_event.getShortName());
                if (acceptor_amount!=0)
                    sb.append("(").append(acceptor_amount).append(")");
            }
            return sb.toString();
        }
    }
    
            
    
    
//    private String getSiteMotif(IntronPlacement P, DNAStrand seq, boolean is_donor)
//    {
//        assert (P.donor_site.type != SiteType.NOT_A_SITE);
//        IntronPlacement true_pos = mapToSequence(P, seq);
//        int site_pos = (is_donor?true_pos.donor_site.pos:true_pos.acceptor_site.pos);
//        byte nuc1 = seq.getNucleotideAt(site_pos);
//        byte nuc2 = seq.isReverseStrand()==is_donor?seq.getNucleotideAt(site_pos-1):seq.getNucleotideAt(site_pos+1);
//        char[] motif = new char[2];
//        if (is_donor)
//        {
//            motif[0]=DNASequence.toChar(nuc1);
//            motif[1]=DNASequence.toChar(nuc2);
//        } else
//        {
//            motif[0]=DNASequence.toChar(nuc2);
//            motif[1]=DNASequence.toChar(nuc1);
//        }
//        return new String(motif);
//    }
    
    
    public class IntronPlacement
    {
        private final SpliceSite donor_site;
        private final SpliceSite acceptor_site;
        
        private IntronPlacement(SpliceSite donor, SpliceSite acceptor)
        {
            this.donor_site = donor;
            this.acceptor_site = acceptor;
            assert (donor.phase==acceptor.phase);
        }
        
        @Override
        public String toString()
        {
            return "I["+donor_site+", "+acceptor_site+"]";
        }
        
        private int getAcceptorOffset(){ return acceptor_site.pos-donor_site.pos;}
        public int getLength(){ return Math.abs(getAcceptorOffset())+1;}
        
        public SpliceSite getDonorSite(){ return donor_site;}
        public SpliceSite getAcceptorSite(){ return acceptor_site;}
        
        public boolean hasIntron(){ return donor_site.type != SiteType.NOT_A_SITE;}
        
        @Override
        public int hashCode()
        {
            return donor_site.hashCode() * 17 + acceptor_site.hashCode();
        }
        
        @Override
        public boolean equals(Object o)
        {
            if (o instanceof IntronPlacement)
            {
                IntronPlacement P = (IntronPlacement) o;
                return donor_site.equals(P.donor_site) && acceptor_site.equals(P.acceptor_site);
            } else
                return super.equals(o);
        }
        
        private Event getEventFrom(IntronPlacement parent)
        {
            if (!hasIntron())
            {
                if (donor_site.pos == MISSING_VALUE.donor_site.pos)
                {
                    return new Event(EventType.MISSING_DATA);
                } else if (parent==null || !parent.hasIntron())
                {
                    return new Event(EventType.NO_INTRON);
                } else
                {
                    return new Event(EventType.INTRON_LOSS);//, parent.getLength());
                }
            } else if (parent == null || !parent.hasIntron())
            {
                return new Event(EventType.INTRON_GAIN);//, getLength());
            } else
            {
                EventType donor_event;
                int donor_shift=this.donor_site.pos-parent.donor_site.pos;
                if (donor_shift==0)
                {
                    donor_event = EventType.DONOR_CONSERVATION;
                } else if (donor_shift*getAcceptorOffset()>0) // (Integer.signum(getAcceptorOffset())==Integer.signum(donor_shift)) 
                {
                    // donor site shifted towards acceptor site
                    donor_event = EventType.DONOR_EXONIZATION;
                } else
                {
                    donor_event = EventType.DONOR_INTRONIZATION;
                }
                EventType acceptor_event;
                int acceptor_shift=this.acceptor_site.pos-parent.acceptor_site.pos;
                if (acceptor_shift==0)
                {
                    acceptor_event = EventType.ACCEPTOR_CONSERVATION;
                } else if (acceptor_shift*getAcceptorOffset()>0) 
                {
                    // acceptor site shifted downstream
                    acceptor_event = EventType.ACCEPTOR_INTRONIZATION;
                } else
                {
                    acceptor_event = EventType.ACCEPTOR_EXONIZATION;
                }
                return new Event(donor_event, Math.abs(donor_shift), acceptor_event, Math.abs(acceptor_shift));
            }
        }
        
        /**
         * Maps the donor and acceptor sites to a sequence. 
         * 
         * @param seq
         * @return 
         */
        public IntronPlacement mapToOrganism(DNAStrand seq)
        {
            if (!hasIntron()) 
                return this;
            else
            {
                int dpos = donor_homologies.projectAbsolutePosition(seq, donor_site.pos);
                int apos = acceptor_homologies.projectAbsolutePosition(seq, acceptor_site.pos);
                SpliceSite donor = new SpliceSite(donor_site.type, dpos, donor_site.phase);
                SpliceSite acceptor = new SpliceSite(acceptor_site.type, apos, acceptor_site.phase);
                return new IntronPlacement(donor, acceptor);
            }
        }
    }
}

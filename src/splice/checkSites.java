package splice;

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

import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import splice.annotation.AnnotatedGenomes;
import splice.annotation.GenePred;
import splice.empirical.Aggregator;
import splice.model.Parser;
import splice.model.SpliceSiteHistory;
import splice.model.SpliceSiteHistory.IntronPlacement;
import splice.model.SpliceSiteHomology;
import splice.model.TreeNode;
import splice.sequence.DNAStrand;
import splice.sequence.ProteinSequence;
import splice.sequence.TranslationTable;

/**
 *
 * Computes ancestral reconstructions.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class checkSites 
{
    
    private static final int DEFAULT_LOSS_PTY = 5;
    private static final int DEFAULT_GAIN_PTY = 12;
    private static final int DEFAULT_REANNOTATE_PTY = 1;
    private static final int DEFAULT_SHIFT_PTY = 5; 
    
    private static final String PROGRAM_NAME = "ssss";
    
    private checkSites(){}
    
    private int pty_loss = DEFAULT_LOSS_PTY;
    private int pty_gain = DEFAULT_GAIN_PTY;
    private int pty_reannotate=DEFAULT_REANNOTATE_PTY;
    private int pty_shift = DEFAULT_SHIFT_PTY;
    private boolean show_alignment=false;

    private static final String REPORT_PREFIX = "#|";
    public void reportLaunch(PrintStream out, String[] args)
    {
        int max_args_reported = 24;
        
        java.util.Date now = java.util.Calendar.getInstance().getTime();
        java.util.Properties Props=System.getProperties();
        String system = Props.getProperty("os.name", "[unknown OS]")+" "+Props.getProperty("os.version","[Unknown version]")+" "+Props.getProperty("os.arch","[Unknown architecture]");
        String java = Props.getProperty("java.vm.name","[Unknown VM]")+" "+Props.getProperty("java.vm.version","[Unknown version]")+" ("+Props.getProperty("java.runtime.version","[Unknwon runtime]")+") "+Props.getProperty("java.vm.info","")+", "+Props.getProperty("java.vm.vendor","[Uknown vendor]");
        out.println(REPORT_PREFIX+getClass().getName()+" "+now);
        out.println(REPORT_PREFIX+"System: "+system);
        out.println(REPORT_PREFIX+"Java: "+java);
        out.println(REPORT_PREFIX+"Current directory: "+Props.getProperty("user.dir","[unknown]"));
        out.print(REPORT_PREFIX+"Arguments:");
        for (int i=0; i<Math.min(max_args_reported,args.length); i++) out.print(" "+args[i]);
        if (args.length>max_args_reported) out.print(" ...["+(args.length-max_args_reported)+" more arguments]");
        out.println();
    }

    private static final String OPTION_LOSS = "-loss";
    private static final String OPTION_GAIN = "-gain";
    private static final String OPTION_REANNOTATE = "-annot";
    private static final String OPTION_SHIFT = "-shift";
    private static final String OPTION_SHOW_ALI = "-showali";
    private static final String OPTION_OUTPUT = "-o";
    private static final String OPTION_SAVE = "-save";
    
    /**
     * File extension used in output Fasta files for protein sequences of modified gene structures
     */
    private static final String EXTENSION_FASTA = ".faa";
    /**
     * File extension used in output GFF3 files for modified gene structures
     */
    private static final String EXTENSION_GFF = ".gff";
    
    private void bumbumbum(String[] args) throws Exception
    {
        if (args.length==0 || args[0].equals("-h")){
            System.err.println("Call as "+getClass().getName()+" [options] org1=g1.fa,org2=g2.fa annot.txt agg.data tree.newick alignment1.fa alignment2.fa ...");
            System.err.println("Checks splice site homologies and reconstructs ancestral status.");
            System.err.println("\torg=g.fa: genomic sequence for organism encoded as 'org' is in Fasta file g.fa");
            System.err.println("\tagg.data: output statistics from "+collectStatistics.class.getName());
            System.err.println("\ttree: Newick format phylogenetic tree for the genomes");
            System.err.println("\talignment.fa: alignment of orthologous protein-sequences in Fasta format, terminal node names match org1, org2, ...");
            System.err.println("\t--- all input files may be URL or local path to a file, and also compressed by gzip (reognized by the .gz extension)");
            System.err.println("Options:");
            System.err.println("\t"+OPTION_LOSS+" x\tintron loss penalty");
            System.err.println("\t"+OPTION_GAIN+" x\tintron gain penalty");
            System.err.println("\t"+OPTION_SHIFT+" x\tsplice site shift penalty");
            System.err.println("\t"+OPTION_REANNOTATE+" x\treannotation penalty");
            System.err.println("\t"+OPTION_SHOW_ALI+" yes|no|true|false\tshow aligned codons and intron sequences");
            System.err.println("\t"+OPTION_OUTPUT+" file\toutput (text messages) redirect");
            System.err.println("\t"+OPTION_SAVE+" directory\twrites GFF and Fasta files for the reannotated gene structures into the given directory");
            
            System.exit(9);
        }
        
        int arg_idx=0;
        
        PrintStream out = System.out;
        
        File save_dir = null;
        while (arg_idx<args.length && args[arg_idx].startsWith("-"))
        {
            String option = args[arg_idx++];
            String value = args[arg_idx++];
            if (OPTION_LOSS.equals(option))
                pty_loss = Integer.parseInt(value);
            else if (OPTION_GAIN.equals(option))
                pty_gain = Integer.parseInt(value);
            else if (OPTION_SHIFT.equals(option))
                pty_shift = Integer.parseInt(value);
            else if (OPTION_REANNOTATE.equals(option))
                pty_reannotate = Integer.parseInt(value);
            else if (OPTION_OUTPUT.equals(option))
                out = "-".equals(value)?System.out:new PrintStream(value);
            else if (OPTION_SHOW_ALI.equals(option))
                show_alignment = value.equalsIgnoreCase("true") || value.equalsIgnoreCase("yes");
            else if (OPTION_SAVE.equals(option))
                save_dir = new File(value);
            else
                throw new IllegalArgumentException("Unknown option "+option);
        }
        
        reportLaunch(out, args);
        
        out.println(REPORT_PREFIX+"Intron loss penalty:  \t"+pty_loss);
        out.println(REPORT_PREFIX+"Intron gain penalty:  \t"+pty_gain);
        out.println(REPORT_PREFIX+"Site shift penalty:   \t"+pty_shift);
        out.println(REPORT_PREFIX+"Site reannotation pty:\t"+pty_reannotate);
        if (save_dir != null)
            out.println(REPORT_PREFIX+"Updated annotations saved in directory "+save_dir+" ("+(save_dir.exists()?"already exists)":"will be created)"));
        
        String genome_fa_list = args[arg_idx++];
        String annotation_file = args[arg_idx++];
        String scoring_file = args[arg_idx++];
        String tree_file = args[arg_idx++];
        
        AnnotatedGenomes annotations = new AnnotatedGenomes();
        List<String> wanted_organisms = annotations.readMultipleGenomes(genome_fa_list);
        annotations.readAnnotations(annotation_file); 
        
        Aggregator all_scoring = Aggregator.readData(new GeneralizedFileReader(scoring_file));

        TreeNode root =Parser.readNewick(new GeneralizedFileReader(tree_file));
        
        Map<String, PrintStream> gene_annotation_files = new HashMap<>();
        
        // make sure we have a directorz to save to if asked
        if (save_dir != null)
        {
            if (!save_dir.exists())
            {
                if (!save_dir.mkdir())
                    throw new java.io.IOException("Cannot create directory "+save_dir);


            }
            
            for (String org: wanted_organisms)
            {
                File org_file = new File(save_dir, org+EXTENSION_GFF);
                PrintStream org_gff = new PrintStream(org_file);
                org_gff.println("##gff-version 3");
                gene_annotation_files.put(org, org_gff);
                
                out.println("#FILE gene structures for "+org+" written into GFF3-format file "+org_file);
            }
        }

        
        int num_sites=0;
        while (arg_idx<args.length)
        {
            String alignment_file = args[arg_idx++];
            int slash_idx = alignment_file.lastIndexOf('/');
            int dot_idx = alignment_file.indexOf('.', slash_idx+1);
            dot_idx = dot_idx==-1?alignment_file.length():dot_idx;
            String alignment_name = alignment_file.substring(slash_idx+1, dot_idx);
            ProteinSequence[] alignment= ProteinSequence.readMultiFasta(alignment_file);
            SpliceSiteHomology H = new SpliceSiteHomology(all_scoring);
            if (H.computeAlignmentColumns(annotations, alignment))
            {
                if (show_alignment)
                {
                    H.reportGenes(out);                
                    H.printAlignment(out);
                }
                
                Map<String, GenePred> gene_structures = new HashMap<>();
                Map<String, String> edit_log = new HashMap<>();
                
                for (int gene_idx=0; gene_idx<alignment.length; gene_idx++)
                {
                    GenePred gene = H.getGene(gene_idx);
                    String org = gene.getChrom().getOrganism();
                    gene_structures.put(org, gene);
                }
                
                for (SpliceSiteHomology.IntronContext context: H.calculateIntronContexts())
                {
                    SpliceSiteHistory history = context.collectSites();
                    if (history != null)
                    { 
                        if (num_sites==0)
                            out.println("#ANCESTRAL\t"+history.getAncestralDescription(root, null));
                        //System.out.println("#*cS.z context "+context);
                        Map<TreeNode, IntronPlacement> reconstruction = history.ancestralParsimony(root, pty_loss, pty_gain, pty_shift, pty_reannotate);
                        out.print("#ANCESTRAL\t"+history.getAncestralDescription(root, reconstruction));
                        out.print("\t");
                        out.print(alignment_name);
                        out.print("\t");
                        out.println(context);
                        
                        history.reportHistory(out, root, reconstruction);
                        
                        for (TreeNode N: root.getTraversal().getLeaves())
                        {
                            IntronPlacement node_state = reconstruction.get(N);
                            String org = N.getName();
                            
                            DNAStrand seq = history.getSequence(org);
                            IntronPlacement annotated_state = history.getAnnotation(seq);
                            IntronPlacement mapped_state = node_state.mapToOrganism(seq);
                            
                            
                            if (node_state.hasIntron())
                            {
                                if (annotated_state == null || !annotated_state.hasIntron())
                                {
                                    int dpos = mapped_state.getDonorSite().getPosition();
                                    int apos = mapped_state.getAcceptorSite().getPosition();
                                    out.println("#REANNOTATE\t"+org+"\tintronify\t"+mapped_state+"\t"+node_state.getDonorSite().getDinucleotideMotif(seq)+">"+node_state.getAcceptorSite().getDinucleotideMotif(seq));
                                    GenePred modified = gene_structures.get(org).editIntronify(dpos, apos);
                                    gene_structures.put(org, modified);
                                    String edit_msg = "intronify("+dpos+".."+apos+")";
                                    if (edit_log.containsKey(org))
                                    {
                                        edit_msg = edit_log.get(org) + ","+edit_msg;
                                    }
                                    edit_log.put(org, edit_msg);
                                } else
                                {
                                    int donor_update = node_state.getDonorSite().getOffsetFromAnnotation(seq);
                                    
                                    if (donor_update!=0)
                                    {
                                        out.println("#REANNOTATE\t"+org+"\tdonor\t"+mapped_state.getDonorSite().getPosition()+"\t"+donor_update+"\t"+node_state.getDonorSite().getDinucleotideMotif(seq)+">"+node_state.getAcceptorSite().getDinucleotideMotif(seq));
                                        GenePred modified = gene_structures.get(org).editDonorShift(mapped_state.getDonorSite().getPosition());
                                        gene_structures.put(N.getName(), modified);
                                        
                                        String edit_msg = "5'shift("+donor_update+")";
                                        if (edit_log.containsKey(org))
                                            edit_msg = edit_log.get(org) + "," + edit_msg;
                                        edit_log.put(org, edit_msg);
                                    }
                                    
                                    int acceptor_update = node_state.getAcceptorSite().getOffsetFromAnnotation(seq);
                                    if (acceptor_update!=0)
                                    {
                                        out.println("#REANNOTATE\t"+N.getName()+"\tacceptor\t"+mapped_state.getAcceptorSite().getPosition()+"\t"+acceptor_update+"\t"+node_state.getAcceptorSite().getDinucleotideMotif(seq)+"<"+node_state.getDonorSite().getDinucleotideMotif(seq));
                                        GenePred modified = gene_structures.get(N.getName()).editAcceptorShift(mapped_state.getAcceptorSite().getPosition());
                                        gene_structures.put(N.getName(), modified);
                                        String edit_msg = "3'shift("+acceptor_update+")";
                                        if (edit_log.containsKey(org))
                                            edit_msg = edit_log.get(org) + "," + edit_msg;
                                        edit_log.put(org, edit_msg);
                                    }
                                }
                            } else
                            {
                                if (annotated_state != null && annotated_state.hasIntron())
                                {
                                    IntronPlacement original_annotation = annotated_state.mapToOrganism(seq);
                                    int dpos = original_annotation.getDonorSite().getPosition();
                                    int apos = original_annotation.getAcceptorSite().getPosition();
                                    out.println("#REANNOTATE\t"+org+"\texonify\t"+original_annotation);
                                    GenePred modified = gene_structures.get(org).editExonify(dpos, apos);
                                    gene_structures.put(org, modified);
                                    String edit_msg = "exonify("+dpos+".."+apos+")";
                                    if (edit_log.containsKey(org))
                                        edit_msg = edit_log.get(org) + "," +edit_msg;
                                    edit_log.put(org, edit_msg);
                                }
                            }
                        }
                        num_sites++;
                    }
                } // for all sites
                
                if (!edit_log.isEmpty())
                {
                    for (String org: edit_log.keySet())
                    {
                        GenePred gene = gene_structures.get(org);
                        
                        String fasta_name = org + "-" + alignment_name+EXTENSION_FASTA;
                        File org_fasta = new File(save_dir, fasta_name);
                        PrintStream fasta_out = new PrintStream(org_fasta);
                        
                        String fasta_defline = gene.getName()+" /reannotate="+edit_log.get(org);
                        String fasta_seq = gene.getTranslation(fasta_defline, TranslationTable.STANDARD_TABLE).toFasta();
                        fasta_out.println(fasta_seq);
                        fasta_out.close();
                        out.println("#FILE aa sequence for modified gene "+gene.getName()+" (alignment "+alignment_name+", genome "+org+") written into "+org_fasta);
                        
                        PrintStream gff_out = gene_annotation_files.get(org);
                        for (String line: gene.toGFF3(PROGRAM_NAME))
                        {
                            gff_out.print(line);
                            gff_out.print(";reannotate=");
                            gff_out.print(edit_log.get(org));
                            gff_out.println();
                        }
                    }
                }
            } else
            {
                // skip sequence
            }
            
            
        } // for all alignment files
        
        if (save_dir != null)
            for (PrintStream gff: gene_annotation_files.values())
                gff.close();
        
        if (out != System.out)
            out.close();
    }
    
    public static void main(String[] args) throws Exception
    {
        (new checkSites()).bumbumbum(args);
    }
    
}

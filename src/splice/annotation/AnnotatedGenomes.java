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
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import splice.GeneralizedFileReader;
import splice.empirical.Counter;
import splice.sequence.DNASequence;
import splice.sequence.DNASequence.Genome;
import splice.sequence.ProteinSequence;

/**
 * Set of genomes with accompanying annotations.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public final class AnnotatedGenomes 
{
    /**
     * Instantiation with empty sets.
     */
    public AnnotatedGenomes()
    {
        this(new HashMap<>(), new HashMap<>());
    }
    
    /**
     * Instantiation with genomes and annotations.
     * 
     * @param genomes maps organisms names to genomes
     * @param annotations maps organisms names to (gene name, gene description) pairs
     */
    private AnnotatedGenomes(Map<String,Genome> genomes, Map<String, Map<String, GenePred>> annotations)
    {
        this.annotations = annotations;
        this.genomes = genomes;
    }
    private final Map<String, Map<String, GenePred>> annotations;
    private final Map<String, Genome> genomes;
    
    public GenePred getAnnotationForGene(ProteinSequence seq)
    {
        String seq_id = seq.getIdent();
        String organism = seq.getOrganism();
        Map<String, GenePred> genome_annotations = annotations.get(organism);
        if (genome_annotations==null)
        {
            for (String org: annotations.keySet())
            {
                if (annotations.get(org).containsKey(seq_id))
                    return annotations.get(org).get(seq_id);
            }
            return null;
        }
        else
            return genome_annotations.get(seq_id);
    }
    
    /**
     * Adds a set of chromosomes to the known genome set.
     * 
     * @param organism organism identifier
     * @param refseq array of chromosome sequences (here we also set organism identifiers for {@link DNASequence#getOrganism() })
     */
    public void addGenome(String organism, DNASequence[] refseq)
    {
        Genome g;
        if (genomes.containsKey(organism))
        {
            g = genomes.get(organism);
            g.addAll(refseq);
        } else
        {
            g = new Genome(refseq);
            genomes.put(organism, g);
        }
        for (DNASequence s:refseq)
        {
            s.setOrganism(organism);
        }
    }
    
    public void addAnnotation(GenePred gene)
    {
        if (gene != null)
        {
            DNASequence chrom = gene.getChrom();
            String organism = chrom.getOrganism();
            Map<String,GenePred> genome_annotations;
            if (annotations.containsKey(organism))
            {
                genome_annotations = annotations.get(organism);
            } else
            {
                genome_annotations = new HashMap<>();
                annotations.put(organism, genome_annotations);
            }
            
            String gene_name = gene.getName();
            if (gene_name != null)
                genome_annotations.put(gene_name, gene);
        }
    }
    
    public Map<String, Genome> getGenomes(){ return genomes;}  
    
    public Genome getGenome(String organism)
    {
        return genomes.get(organism);
    }
    
    public Collection<GenePred> getAllAnnotations(String organism)
    {
        Map<String, GenePred> org_annot = annotations.get(organism);
        return org_annot.values();
    }
    
    /**
     * Reads multiple genomes (one Fasta file per genome)
     * 
     * @param genome_fa_list comma separated list of <var>organism</var><code>:</code><var>file</var> entries
     * @return list of organism codes
     * @throws java.io.IOException if file reading fails
     */
    public List<String> readMultipleGenomes(String genome_fa_list) throws java.io.IOException
    {
        List<String> listed_organisms = new ArrayList<>();
        for (String f: genome_fa_list.split(","))
        {
            int split_at = f.indexOf('=');
            String organism = (split_at==-1)?"":f.substring(0,split_at);
            String fasta_file = f.substring((split_at==-1)?0:split_at+1,f.length());
            //System.out.println("#*AG.rMG "+organism+" genome.fa:"+fasta_file);
            DNASequence[] refseq = DNASequence.readMultiFasta(fasta_file);
            addGenome(organism, refseq);
            listed_organisms.add(organism);
        }
        return listed_organisms;
    }

    private void readAnnotationsSingle(String annotation_file) throws java.io.IOException
    {
        BufferedReader R=new GeneralizedFileReader(annotation_file);
        String line ;
        do
        {
            line = R.readLine();
            if (line != null && !line.startsWith("#"))
            {
                line = line.trim();
                GenePred gene = GenePred.parseTabbed(line, getGenomes());
                addAnnotation(gene);
            }
        } while (line != null);
        R.close();   
    }    
    
    private void readAnnotationsSingle(String annotation_file, Genome genome) throws java.io.IOException
    {
        BufferedReader R=new GeneralizedFileReader(annotation_file);
        String line ;
        do
        {
            line = R.readLine();
            if (line != null && !line.startsWith("#"))
            {
                line = line.trim();
                GenePred gene = GenePred.parseTabbed(line, genome);
                addAnnotation(gene);
            }
        } while (line != null);
        R.close();   
    }
    
    /**
     * Reads annotations from a file or a list of files. 
     * List elements are separated by comma. A genome may be specified for each annotation file 
     * in 'org=annot.txt' syntax. 
     * 
     * @param annotation_file_or_list list of annotation files 
     * @throws java.io.IOException if reading fails
     */
    public void readAnnotations(String annotation_file_or_list) throws java.io.IOException
    {
        for (String f: annotation_file_or_list.split(","))
        {
            int split_at = f.indexOf('=');
            if (split_at == -1)
            {
                readAnnotationsSingle(f);
            } else
            {
                String organism = f.substring(0,split_at);
                String annotation_file = f.substring(split_at+1);
                Genome genome = getGenome(organism);
                readAnnotationsSingle(annotation_file, genome);
            }
        }
    }
    
    
    private static final int GFF_SEQID  =0;
    private static final int GFF_SOURCE =1;
    private static final int GFF_TYPE   =2;
    private static final int GFF_START  =3;
    private static final int GFF_END    =4;
    private static final int GFF_SCORE  =5;
    private static final int GFF_STRAND =6;
    private static final int GFF_PHASE  =7;
    private static final int GFF_ATTR   =8;
    
    public void readAnnotationsGFF3(String annotation_file,  String gff_feature_type, String prefix_identifier_with) 
            throws java.io.IOException
    {
        readGFFLike(annotation_file, gff_feature_type,
                "ID", '=', prefix_identifier_with);
    }
    
    public void readAnnotationsGTF(String annotation_file,String gtf_feature_type, String prefix_identifier_with) 
            throws java.io.IOException
    {
        readGFFLike(annotation_file, gtf_feature_type, 
                "transcript_id", ' ', prefix_identifier_with);
    }
    
    public void readAnnotationsGFF2(String annotation_file, String gff_feature_type, String protein_id_attribute, String prefix_identifier_with)
            throws java.io.IOException
    {
        readGFFLike(annotation_file, gff_feature_type,
                protein_id_attribute, ' ', prefix_identifier_with);
    }
    
    /**
     * Reading from multiple tab-delimited GFF2/GTF/GFF3 files.
     * 
     * @param annotation_file_or_list comma-separated list of org=annotation.txt entries
     * @param gene_annotation_type format 
     * @param gene_identifier_tag tag used to identify gene id
     * @param key_value_separator space or = separating key/value pairs in attribute column
     * @param prefix_identifier_with prepended value (if '|', then organism code is added)
     * @throws java.io.IOException 
     */
    public void readGFFLike(String annotation_file_or_list, String gene_annotation_type, 
            String gene_identifier_tag, 
            char key_value_separator,
            String prefix_identifier_with) throws java.io.IOException
    {
        for (String f: annotation_file_or_list.split(","))
        {
            int split_at = f.indexOf('=');
            if (split_at == -1)
            {
                throw new IllegalArgumentException("Must specify organism for annotation files in \"org=file.txt\" syntax: got "+f);
            } else
            {
                String organism = f.substring(0,split_at);
                String annotation_file = f.substring(split_at+1);
                Genome genome = getGenome(organism);
                String org_prefix = "|".equals(prefix_identifier_with)?organism+"|":prefix_identifier_with;
                readGFFLike(annotation_file, genome, gene_annotation_type, gene_identifier_tag, key_value_separator, org_prefix);
            }
        }
    }
    
    private void readGFFLike(String annotation_file, Genome genome, 
            String gene_annotation_type, 
            String gene_identifier_tag, 
            char key_value_separator,
            String prefix_identifier_with) throws java.io.IOException
    {
        BufferedReader R=new GeneralizedFileReader(annotation_file);
        Map<String,List<String>> annotated_intervals = new HashMap<>();
        
        String line ;
        do
        {
            line = R.readLine();
            if (line != null && !line.startsWith("#"))
            {
                line = line.trim();
                String[] fields = line.split("\\t");
                if (gene_annotation_type.equals(fields[GFF_TYPE]))
                {
                    Map<String, String> attributes = this.parseGFFAttributes(fields[GFF_ATTR], key_value_separator);
                    String gene_ident = attributes.getOrDefault(gene_identifier_tag, "unknown");

                    if ("unknown".equals(gene_ident))
                    {
                        System.err.println("#* SKIPPING --- NO IDENTIFIER (looking for "+gene_identifier_tag+": "+line);
                    } else if (annotated_intervals.containsKey(gene_ident))
                    {
                        List<String> exons = annotated_intervals.get(gene_ident);
                        exons.add(line);
                    } else
                    {
                        List<String> exons = new ArrayList<>();
                        exons.add(line);
                        annotated_intervals.put(gene_ident, exons);
                    }
                }
            }
        } while (line != null);
        R.close();   
        
        for (String gene_ident: annotated_intervals.keySet())
        {
            List<String> exons = annotated_intervals.get(gene_ident);
            int num_exons = exons.size();
            int[] exon_starts = new int[num_exons];
            int[] exon_ends = new int[num_exons];
            
            char strand = '.';
            DNASequence chromo=null;
            
            for (int eidx=0; eidx<num_exons; eidx++)
            {
                String[] gff_fields = exons.get(eidx).split("\\t");
                
                if (strand=='.')
                {
                    strand = gff_fields[GFF_STRAND].charAt(0);
                    chromo = genome.getChromosome(gff_fields[GFF_SEQID]);
                } else
                {
                    if (!(gff_fields[GFF_STRAND].charAt(0) == '.' || strand == gff_fields[GFF_STRAND].charAt(0)))
                    {
                        for (int j=0; j<num_exons; j++)
                        {
                            System.err.println("#*BAD STRAND "+gene_ident+" (exon "+j+(j==eidx?"**":"")+") "+exons.get(j)+"\t// strand "+strand+"\there "+gff_fields[GFF_STRAND]);
                        }
                    }
                    assert (gff_fields[GFF_STRAND].charAt(0) == '.' || strand == gff_fields[GFF_STRAND].charAt(0));
                    assert (chromo == genome.getChromosome(gff_fields[GFF_SEQID]));
                }
                int s = exon_starts[eidx] = Integer.parseInt(gff_fields[GFF_START])-1; // 1-based inclusive -> 0-based inclusive start
                int e = exon_ends[eidx]= Integer.parseInt(gff_fields[GFF_END]); // 1-based inclusive -> 0-based exclusive
            }
            
            // sort exons
            if (num_exons>1)
            {
                Integer[] sorted_exon_indexes = new Integer[num_exons];
                for (int eidx=0; eidx<num_exons; eidx++)
                {
                    sorted_exon_indexes[eidx] = eidx;
                }
                Arrays.sort(sorted_exon_indexes, new java.util.Comparator<Integer>()
                {

                    @Override
                    public int compare(Integer o1, Integer o2) 
                    {
                        return Integer.compare(exon_starts[o1], exon_starts[o2]);
                    }
                });
                
                int[] tmp_val = new int[num_exons];
                System.arraycopy(exon_starts, 0, tmp_val, 0, num_exons);
                for (int eidx=0; eidx<num_exons; eidx++)
                    exon_starts[eidx] = tmp_val[sorted_exon_indexes[eidx]];
                System.arraycopy(exon_ends, 0, tmp_val, 0, num_exons);
                for (int eidx=0; eidx<num_exons; eidx++)
                    exon_ends[eidx] = tmp_val[sorted_exon_indexes[eidx]];
            }
            int cds_start = exon_starts[0];
            int cds_end = exon_ends[num_exons-1];

            // cut class:ident style 
            String protein_identifier = gene_ident;
            if (gene_ident.startsWith("\"")) // strip quotes
            {
                int last_quote = gene_ident.lastIndexOf("\"");
                protein_identifier = gene_ident.substring(1, last_quote);
            }
            
            int colon_pos = protein_identifier.indexOf(':');
            if (colon_pos >=0 )
            {
                protein_identifier = protein_identifier.substring(colon_pos+1);
            }
            if (prefix_identifier_with != null)
                protein_identifier = prefix_identifier_with+protein_identifier;
            GenePred gene = new GenePred(protein_identifier, chromo, strand, exon_starts, exon_ends, cds_start, cds_end);
            addAnnotation(gene);
        }
    }
    
    private Map<String,String> parseGFFAttributes(String field_text, char key_value_separator) //throws java.net.URISyntaxException
    {
        Map<String, String> key_value_pairs = new HashMap<>();
        
        String[] entries = field_text.split("\\;\\s*");
        for (String pair : entries)
        {
            int sep_idx = pair.indexOf(key_value_separator);
            if (sep_idx==-1)
                throw new IllegalArgumentException("Attribute field with no key/value separator: "+pair+" (all: "+entries+")");
            String key = pair.substring(0,sep_idx);
            String value = pair.substring(sep_idx+1);
//            if (key_value_separator == '=')
//            {
//                value = (new java.net.URI(value)).toString();
//            }
            key_value_pairs.put(key, value);
        }
        return key_value_pairs;
    }
    
    
    
//    /**
//     * Filters protein sequences with known gene annotations.
//     * 
//     * @param sequences array of protein sequences
//     * @return list of known sequences
//     */
//    public List<ProteinSequence> keepAnnotatedSequences(ProteinSequence[] sequences)
//    {
//        List<ProteinSequence> kept_sequences = new ArrayList<>();
//        for(ProteinSequence pseq: sequences)
//        {
//            if (getAnnotationForGene(pseq)!=null)
//                kept_sequences.add(pseq);
//        }
//        return kept_sequences;
//    }
    
//    private static void incrementCount(Map<Integer, Integer> counts, int val)
//    {
//        if (counts.containsKey(val))
//        {
//            int n = counts.get(val);
//            counts.put(val, n+1);
//        } else
//        {
//            counts.put(val, 1);
//        }
//    }
    
    private void reportStatistics(String org, java.io.PrintStream out)
    {
        Counter<Integer> exon_counts=new Counter<>();
        Counter<Integer> exon_lengths = new Counter<>();
        Counter<Integer> intron_lengths = new Counter<>();
        for (GenePred gene: getAllAnnotations(org))
        {
            int num_exons = gene.getExonCount();
            exon_counts.increment(num_exons);

            {
                int e_idx=0;
                int e_end = Integer.MAX_VALUE; // will not be used
                do 
                {
                    int e_start = gene.getExonStart(e_idx);
                    if (e_idx>0)
                    {
                        intron_lengths.increment(e_start-e_end);
                    }
                    e_end = gene.getExonEnd(e_idx);
                    exon_lengths.increment(e_end-e_start);
                    e_idx++;
                } while (e_idx<num_exons);
            }
        }
        
        out.println("# Exon count distribution for "+org);
        Integer[] sorted_ecounts = exon_counts.keySet().toArray(new Integer[0]);
        Arrays.sort(sorted_ecounts);
        for (Integer num_e: sorted_ecounts)
        {
            out.printf("%d\t%d\n", num_e, exon_counts.getCount(num_e));
        }
        
        out.println("# Exon length distribution for "+org);
        Integer[] sorted_elengths = exon_lengths.keySet().toArray(new  Integer[0]);
        Arrays.sort(sorted_elengths);
        for (Integer len_e: sorted_elengths)
        {
            out.printf("%d\t%d\n", len_e, exon_lengths.getCount(len_e));
        }
        
        out.println("# Intron length distribution for "+org);
        Integer[] sorted_ilengths = intron_lengths.keySet().toArray(new Integer[0]);
        Arrays.sort(sorted_ilengths);
        for (Integer len_i: sorted_ilengths)
        {
            out.printf("%d\t%d\n", len_i, intron_lengths.getCount(len_i));
        }
    }
    
    /** 
     * Test code: reads annotations, and calculates statistics on exon count, intron and exon lengths.
     * 
     * @param args Command-line arguments
     * @throws Exception if something went wrong
     */
    public static void main(String[] args) throws Exception 
    {
        if (args.length==0 || args[0].equals("-h"))
        {
            System.err.println("Tests: call as $0 org1=g1.fa,org2=g2.fa annot.txt ");
            System.err.println("Reads gene annotations from file and reports statistics (exon counts, exon and intron lengths).");
            System.exit(9);
        }
        int arg_idx=0;
        String genome_fa_list = args[arg_idx++];
        String annotation_file = args[arg_idx++];
        
        AnnotatedGenomes A = new AnnotatedGenomes();
        List<String> wanted_organisms = A.readMultipleGenomes(genome_fa_list);
        A.readAnnotations(annotation_file);
        
        for (String org: wanted_organisms)
        {
            A.reportStatistics(org, System.out);
        }
        
    }
}
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

import java.io.PrintStream;
import java.util.List;
import splice.annotation.AnnotatedGenomes;
import splice.annotation.GenePred;

/**
 * Converts gene annotation files into a format we use.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class extractAnnotations 
{
    private extractAnnotations(){}

    private static final String OPTION_OUTPUT = "-o";
    private static final String OPTION_TYPE = "-feature";
    private static final String OPTION_SEPARATOR = "-key-value_separator";
    private static final String OPTION_PREFIX = "-prefix";
    private static final String OPTION_IDENTIFIER = "-ident";

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
    
    
    public void zimezum(String[] args) throws Exception
    {
        if (args.length==0 || args[0].equals("-h"))
        {
            System.err.println("Call as "+getClass().getName()+" [options] type org1=g1.fa,org2=g2.fa annotations ");
            System.err.println("Writes annotation file from GFF3, GTF and JGI/GFF formats, for one or more genomes.");
            System.err.println("Types (case-insensitive):");
            System.err.println("\tgff3\tGFF3 format input file(s), protein identifiers taken from ID attribute in column 9");
            System.err.println("\tgtf\tGTF format input file(s), protein identifiers taken from transcript_id attribute in column 9");
            System.err.println("\tjgi\tJGI's GFF2 format, protein identifiers taken from ProteinId attribute in column 9");
            System.err.println("\tgff\t\"Generic\" GFF format, specifiy attribute for protein identifiers");
            System.err.println("Genome files: comma separated list of organism abbreviation - Fasta file pairs");
            System.err.println("\torg=g.fa: genomic sequence for organism encoded as 'org' is in Fasta file g.fa");
            System.err.println("Annotations: comma-separated list of organism abbreviations - annotation file pairs");
            System.err.println("\torg=annot.txt: gene structure annotations for organism encoded as 'org' is in annot.txt");
            System.err.println("Options:");
            System.err.println("\t"+OPTION_TYPE+" feature_type\tfeature type in GFF formats (default CDS)");
            System.err.println("\t"+OPTION_IDENTIFIER+" attribute key for gene identifier");
            System.err.println("\t"+OPTION_PREFIX+" prefix for constructing gene identifier from the value associated with "+OPTION_IDENTIFIER+"; use '|' in order to prepend organism code and vertical bar separator");
            System.err.println("\t"+OPTION_SEPARATOR+" one-character separator between kety and value pairs ('=' in GFF3, ' ' in GTF and GFF2)");
            System.err.println("\t"+OPTION_OUTPUT+" file\toutput redirect");
            
            System.exit(9);
        }
        
        int arg_idx=0;
        
        PrintStream out = System.out;
        String feature_type = "CDS";
        String protein_id_attribute = null;
        char key_value_separator = ' ';
        String prefix_id_with = null;
        
        // set columns in anticipation of genepred format
        GenePred.setAccnoColumn(0);
        GenePred.setChromoColumn(1);
        GenePred.setExonColumns(8,9);
        GenePred.setStrandColumn(2);
        
        while (arg_idx<args.length && args[arg_idx].startsWith("-"))
        {
            int column_indices_parsed = GenePred.parseCommandLineOptions(args, arg_idx);
            if (column_indices_parsed == arg_idx)
            {
                String option = args[arg_idx++];
                String value = (arg_idx==args.length?null:args[arg_idx++]);
                if (OPTION_OUTPUT.equals(option))
                    out = "-".equals(value)?System.out:new PrintStream(value);    
                else if (OPTION_IDENTIFIER.equals(option))
                    protein_id_attribute = value;
                else if (OPTION_TYPE.equals(option))
                    feature_type = value;
                else if (OPTION_SEPARATOR.equals(option))
                    key_value_separator = value.charAt(0);
                else if (OPTION_PREFIX.equals(option))
                    prefix_id_with = value;
                else
                    throw new IllegalArgumentException("Unknown option "+option);
            } else
                arg_idx = column_indices_parsed;
        }
        
        reportLaunch(out,args);
        out.println(REPORT_PREFIX+"Feature type                :\t"+feature_type);
        out.println(REPORT_PREFIX+"Protein identifier attribute:\t"+(protein_id_attribute == null?"(default)":protein_id_attribute));
        
        String annotation_type = args[arg_idx++].toLowerCase();
        
        String genome_fa_list = args[arg_idx++];
        String annotation_file = args[arg_idx++];
        
        AnnotatedGenomes annotations = new AnnotatedGenomes();
        List<String> wanted_organisms = annotations.readMultipleGenomes(genome_fa_list);  
        
        if ("gff3".equals(annotation_type))
        {
            annotations.readAnnotationsGFF3(annotation_file, feature_type, prefix_id_with);
        } else if ("gtf".equals(annotation_type))
        {
            annotations.readAnnotationsGTF(annotation_file, feature_type, prefix_id_with);
        } else if ("jgi".equals(annotation_type))
        {
            annotations.readAnnotationsGFF2(annotation_file, feature_type, "proteinId", prefix_id_with);
        } else if ("gff".equals(annotation_type))
        {
            if (protein_id_attribute == null)
                throw new IllegalArgumentException("Specifiy "+OPTION_IDENTIFIER+" with generic GFF");
            annotations.readGFFLike(annotation_file, feature_type, protein_id_attribute, key_value_separator, prefix_id_with);
        } else if ("genepred".equals(annotation_type)) 
        {
            annotations.readAnnotations(annotation_file);
        }
        
        out.println("# seq\tgene\taccno\tCDS\tstrand\tlength\tipos\tilen\tgi\torganism\tleftmost\tGP-exonstart\tGP-exonend");
        for (String org:wanted_organisms)
        {
            for (GenePred gene: annotations.getAllAnnotations(org))
            {
                out.println(gene.toAnnotationString());
            }
        }
    }
    
    
    public static void main(String[] args) throws Exception
    {
        (new extractAnnotations()).zimezum(args);
    }
    
}

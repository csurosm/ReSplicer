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
import splice.empirical.Aggregator;
import splice.sequence.ProteinSequence;

/**
 *
 * Calculates statistics for the probabilistic models. Output 
 * is used in {@link checkSites}.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class collectStatistics 
{
    private collectStatistics(){}
    
    private static final String OPTION_OUTPUT = "-o";
    
    public void recefice(String[] args) throws Exception
    {
        if (args.length==0 || args[0].equals("-h"))
        {
            System.err.println("Call as "+getClass().getName()+" [options] org1=g1.fa,org2=g2.fa annot.txt ali1.fa ali2.fa ... ");
            System.err.println("Collects statstics about sequnce and structure in ortholog genes.");
            System.err.println("\torg=g.fa: genomic sequence for organism encoded as 'org' is in Fasta file g.fa");
            System.err.println("\talignment.fa: alignment of orthologous protein-sequences in Fasta format, terminal node names match org1, org2, ...");
            System.err.println("\t--- all input files may be URL or local path to a file, and also compressed by gzip (reognized by the .gz extension)");
            System.err.println("Options:");
            System.err.println("\t"+OPTION_OUTPUT+" file\toutput redirect");
            
            System.exit(9);
        }
        
        int arg_idx=0;
        
        PrintStream out = System.out;
        while (arg_idx<args.length && args[arg_idx].startsWith("-"))
        {
            int column_indices_parsed = GenePred.parseCommandLineOptions(args, arg_idx);
            if (column_indices_parsed == arg_idx)
            {
                String option = args[arg_idx++];
                String value = (arg_idx==args.length?null:args[arg_idx++]);
                if (OPTION_OUTPUT.equals(option))
                    out = "-".equals(value)?System.out:new PrintStream(value);               
                else
                    throw new IllegalArgumentException("Unknown option "+option);
            } else
                arg_idx = column_indices_parsed;
        }
        
        String genome_fa_list = args[arg_idx++];
        String annotation_file = args[arg_idx++];
        
        AnnotatedGenomes annotations = new AnnotatedGenomes();
        List<String> wanted_organisms = annotations.readMultipleGenomes(genome_fa_list);
        annotations.readAnnotations(annotation_file);  

        Aggregator agg = new Aggregator(wanted_organisms, annotations);

        while (arg_idx<args.length)
        {
            String alignment_file = args[arg_idx++];
            agg.countAlignment(annotations, ProteinSequence.readMultiFasta(alignment_file));
        }
        agg.writeData(out);

        if (out != System.out)
            out.close();

    }
    
    public static void main(String[] args) throws Exception
    {
        (new collectStatistics()).recefice(args);
    }
    
}

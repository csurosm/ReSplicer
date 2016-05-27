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

package splice.sequence;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import splice.GeneralizedFileReader;

/**
 *
 * Class for protein sequence representation, with gaps allowed.
 * An instance stores the residue sequence and the accompanying defline for the Fasta file
 * (not the &gt;). 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public final class ProteinSequence extends FastaSequence implements Iterable<AminoAcid>
{
    private final List<AminoAcid> residues;
//    private final String defline;
    
    /**
     * Instantiation with empty sequence. 
     * 
     * @param defline definition line (without starting &gt;)
     */
    public ProteinSequence(String defline)
    {
        super(defline);
        residues = new ArrayList<>();
    }

    /**
     * Adds the contained standard IUBMB residues to the end of the sequence.
     * 
     * @param seq sequence to append at end
     */
    @Override
    public void append(String seq)
    {
        for (int pos=0;pos<seq.length();pos++)
        {
            char c = seq.charAt(pos);
            AminoAcid r = AminoAcid.getResidue(c);
            residues.add(r);
        }
    }

    /**
     * Adds a single residue at the end.
     * 
     * @param c residue to add
     */
    public void append(char c)
    {
        AminoAcid r = AminoAcid.getResidue(c);
        residues.add(r);
    }
    
    /**
     * Adds a single residue at the end of the sequence. 
     * @param aa residue to add
     */
    public void append(AminoAcid aa)
    {
        residues.add(aa);
    }
    
//    /**
//     * Definition line from the input Fasta file.
//     * 
//     * @return 
//     */
//    public String getDefLine()
//    {
//        return defline;
//    }
//    
//    /**
//     * The sequence identifier comes from the defline up to the first space character. 
//     * @return 
//     */
//    public String getIdent()
//    {
//        int first_white_space = defline.indexOf(' ');
//        if (first_white_space<0)
//            return defline;
//        else 
//            return defline.substring(0,first_white_space);
//    }

    /**
     * Amino acid at the specified position
     * 
     * @param pos 0-based position in the sequence
     * @return amino acid residue in this position
     */
    public AminoAcid getResidueAt(int pos)
    {
        return residues.get(pos);
    }
    
    /**
     * Number of residues in the sequence.
     * 
     * @return sequence length
     */
    public int getLength() 
    {
        return residues.size();
    }

    
    /**
     * Reads multiple protein sequences contained in a Fasta file. 
     * 
     * @param file_name local or remote (using <code>http:</code> or <code>ftp:</code> file name, possible gzipped (<code>.gz</code> extension) 
     * @return amino acid sequences in the file
     * @throws IOException if reading fails
     */
    public static ProteinSequence[] readMultiFasta(String file_name) throws IOException
    {
        List<ProteinSequence> sequences = new ArrayList<>();
        // ensuring a very general range of input file specifications: URLs, compression.
        BufferedReader fileReader    = new GeneralizedFileReader(file_name);
        String         line;
        String         id;
        //StringBuffer   contentBuffer = new StringBuffer();
        // current sequence is collected in a DNASequence object
        ProteinSequence current_seq = null;


        // do not catch expections here: propagate to caller and handle it there if necessary
        // (e.g., by informing the user with a well-formatted message...)
        //
        // try {

        do
        {
            line = fileReader.readLine();
            if (line != null)
            {
                line = line.trim();
                if (line.isEmpty())
                {
                    continue;
                }
                char firstChar = line.charAt(0);

                if (firstChar == '>')
                {
                    //addToSequenceList(id, contentBuffer);

                    // now can get the new id > ..
                    id = line.substring(1).trim();
                    current_seq = new ProteinSequence(id);
                    // add to our list immediately
                    sequences.add(current_seq);
                } else if (firstChar == ';')
                {
                    // comment line, skip it
                } else if (current_seq != null)
                {
                    // carry on reading sequence content
                    current_seq.append(line); // .trim()); // already trimmed
                }
            }
        } while (line != null);
        fileReader.close();

        return sequences.toArray(new ProteinSequence[0]);
    }

    /**
     * Iterator over the residues in the sequence.
     * 
     * @return iterator over residues
     */
    @Override
    public Iterator<AminoAcid> iterator() 
    {
        return residues.iterator();
    }

    /**
     * Converts the sequence to a string that can be directly written into a Fasta file.
     * 
     * @param chars_per_line how many characters are allowed per sequence line 
     * @return formatted Fasta-conform representation (defline with &gt; prefix and sequence)
     */
    public String toFasta(int chars_per_line)
    {
        StringBuilder sb = new StringBuilder(">");
        sb.append(getDefLine());
        sb.append("\n");
        int pos = 0;
        int len = getLength();
        while (pos<len)
        {
            AminoAcid r = getResidueAt(pos++);
            sb.append(r.getCode());
            if (pos<len && (pos % chars_per_line == 0))
            {
                sb.append("\n");
            }
        }
        return sb.toString();
    }
    
    /**
     * Fasta representation with default number of characters per line (50). 
     * 
     * @return Fasta-encoded string representation for the sequence
     */
    public String toFasta(){ return toFasta(50);}
}

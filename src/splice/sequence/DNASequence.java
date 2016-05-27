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
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import splice.GeneralizedFileReader;

/**
 * A class for storing DNA sequences. Every sequence has a name in addition
 * to the sequence of nucleotides. The actual sequence is stored by using one
 * byte per nucleotide.
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 *
 */
public class DNASequence extends FastaSequence implements NucleotideEncoding, Comparable<DNASequence>
{
    /**
     * Static array for char -&gt; index conversions
     */
    private static final int[]  char2index;
    /**
     * Static array for index -&gt; char conversions
     */
    private static final char[] index2char;

    /**
     * Static array for char -&gt; byte conversions
     */
    private static final byte[] char2byte;
    /**
     * Static array for byte -&gt; char conversions
     */
    private static final char[] byte2char;

    private static final byte[] complements;
    
    /**
     * Initialization of the static arrays used in conversions
     */
    static
    {
        // char <-> index conversion
        char2index=new int[128];
        java.util.Arrays.fill(char2index,-1);
        char2index[CHAR_A]=IDX_A;
        char2index[CHAR_C]=IDX_C;
        char2index[CHAR_G]=IDX_G;
        char2index[CHAR_T]=IDX_T;
        char2index[CHAR_U]=IDX_T;

        index2char=new char[5];
        index2char[IDX_A]=CHAR_A;
        index2char[IDX_C]=CHAR_C;
        index2char[IDX_G]=CHAR_G;
        index2char[IDX_T]=CHAR_T;
        index2char[IDX_GAP]=CHAR_GAP;

        // byte<-> char conversions
        byte2char = new char[16];
        byte2char[BYTE_A]=CHAR_A;
        byte2char[BYTE_C]=CHAR_C;
        byte2char[BYTE_G]=CHAR_G;
        byte2char[BYTE_T]=CHAR_T;
        byte2char[BYTE_N]=CHAR_N;
        byte2char[BYTE_R]=CHAR_R;
        byte2char[BYTE_K]=CHAR_K;
        byte2char[BYTE_M]=CHAR_M;
        byte2char[BYTE_Y]=CHAR_Y;
        byte2char[BYTE_W]=CHAR_W;
        byte2char[BYTE_S]=CHAR_S;
        byte2char[BYTE_B]=CHAR_B;
        byte2char[BYTE_D]=CHAR_D;
        byte2char[BYTE_H]=CHAR_H;
        byte2char[BYTE_V]=CHAR_V;
        byte2char[BYTE_ERROR]=' ';// also used as spacer character
        char2byte = new byte[128];
        java.util.Arrays.fill(char2byte,BYTE_ERROR);
        for (byte b=0; b<16; b++)
            char2byte[byte2char[b]]=b;
        char2byte[CHAR_U]=BYTE_T;
        char2byte[CHAR_X]=BYTE_N;

        // byte complements
        complements = new byte[16];

        complements[BYTE_A]=BYTE_T;
        complements[BYTE_T]=BYTE_A;
        complements[BYTE_C]=BYTE_G;
        complements[BYTE_G]=BYTE_C;
        complements[BYTE_W]=BYTE_W;
        complements[BYTE_S]=BYTE_S;
        complements[BYTE_N]=BYTE_N;
        complements[BYTE_R]=BYTE_Y;
        complements[BYTE_Y]=BYTE_R;
        complements[BYTE_K]=BYTE_M;
        complements[BYTE_M]=BYTE_K;
        complements[BYTE_B]=BYTE_V;
        complements[BYTE_D]=BYTE_H;
        complements[BYTE_H]=BYTE_D;
        complements[BYTE_V]=BYTE_B;
        complements[BYTE_ERROR]=BYTE_ERROR;
    }

    /**
     * Converts a byte-encoding into 0..3 indexes.
     *
     * It is assumed that the nucleotide is not ambiguous,
     * otherwise the conversion is bogus (for the sake of speed, it is not verified whether the nuc
     * is ambiguous or not).
     * @param nuc byte-encoded nucleotide
     * @return integer (0-3) for ACGT
     */
    public static int toIndex(byte nuc)
    {
        assert (!isAmbiguous(nuc));

        return nuc & 3;
    }

    /**
     * Converts a nucleotide(char) into a nucleotide(index).
     *
     * @param nuc (upper-case) A,C,G, or T
     * @return index 0..3 for real nucleotide, otherwise -1
     */
    public static int toIndex(char nuc)
    {
        return char2index[nuc];
    }

    public static final char toChar(int idx)
    {
        return index2char[idx];
    }

    /**
     * Converts from character to byte encoding.
     * 
     * @param nuc upper-case nucleotide code 
     * @return character decoding of the nucleotide encoded in the lower 4 bits
     */
    public static final char toChar(byte nuc){
        return byte2char[nuc];
    }
    
    public static boolean isAmbiguous(char c)
    {
        return isAmbiguous(char2byte[c]);
    }

    public static boolean isAmbiguous(byte nuc)
    {
        boolean a=(nuc & 12) != 4;
        //System.out.println("#**DNAS.iA "+nuc+"/"+toChar(nuc)+" "+a);
        return a;
    }

    /**
     * @param nuc one of CHAR_A, CHAR_G, CHAR_C, CHAR_T, CHAR_U, CHAR_R, CHAR_Y
     * @return CHAR_R (purine) or CHAR_Y (pyrimidine); or CHAR_N if cannot be determined
     */
    public static final char toRY(char nuc)
    {
        if (nuc==CHAR_A || nuc == CHAR_G || nuc==CHAR_R)
            return CHAR_R;
        else if (nuc==CHAR_T || nuc==CHAR_U || nuc==CHAR_C || nuc==CHAR_Y)
            return CHAR_Y;
        else
            return CHAR_N;
    }


    /**
     * Creates a new DNASequence object with no actual sequence in it, only a name.
     *
     * @param name the sequence name
     */
    public DNASequence(String name)
    {
        super(name);
    }

    /**
     * Creates a new DNASequence object with no actual sequence in it, and empty String name
     * (not null).
     */
    public DNASequence()
    {
        this("");
    }

    /**
     * Instantiation with the same name and content as another
     * sequence. The data structure for the sequence is shared
     * (adding nucleotides to one will add nucleotides to the other too)
     *
     * @param orig original DNA sequence object
     */
    protected DNASequence(DNASequence orig)
    {
        super(orig.getDefLine());
//        setName(orig.name);
        this.sequence = orig.sequence;

    }

            
//    /**
//     * Sequence name
//     */
//    private String name;
    /**
     * Data structure storing the actual sequence
     */
    private BigByteArray sequence;

//    /**
//     * @return The name of this DNA sequence (defline in original Fasta file)
//     */
//    public String getName(){return name;}
//    /**
//     * Sets the name to the given String.
//     * @param name the new name
//     */
//    public final void setName(String name){this.name=name;}

//    /**
//     * Identifier constructed from the first characters of the name up to the first space
//     * 
//     * @return 
//     */
//    public final String getIdent()
//    {
//        int space = name.indexOf(' ');
//        if (space >=0)
//            return name.substring(0,space);
//        else 
//            return name;
//    }
    
//    /**
//     * Equality by name only. 
//     * 
//     * @param o
//     * @return 
//     */
//    @Override
//    public boolean equals(Object o)
//    {
//        if (o instanceof DNASequence)
//        {
//            DNASequence s = (DNASequence)o;
//            return name.equals(s.name);
//        } else
//        {
//            return super.equals(o);
//        }
//    }
//    
//    @Override
//    public int hashCode()
//    {
//        return name.hashCode();
//    }
//    
    /**
     * Access to a position in the sequence.
     *
     * @param pos 0-based position along the sequence
     * @return byte-encoded nucleotide at the given position (index boundaries are not checked)
     */
    public byte getNucleotideAt(int pos)
    {
        return sequence.get(pos);
    }

    /**
     * Number of nucleotides in the sequence.
     *
     * @return Sequence length
     */
    public int getLength()
    {
        return sequence.size();
    }


    /**
     * Appends a DNA sequence at the end. If the
     * sequence was not allocated yet, it is allocated here.
     * Otherwise, the proper capacity is ensured, and
     * the byte-equivalents of the character-encoded nucleotides
     * are added at the end one by one.
     * 
     * @param seq character-encoded DNA sequence; automatically converted to upper case
     */
    public void append(String seq)
    {
        int seq_len = seq.length();
        if (sequence == null) sequence = new BigByteArray(seq_len);
        if (sequence.ensureCapacity(sequence.size()+seq_len))
        {
            //System.out.println("#**DNAS.app "+name+"\t"+sequence.capacity+"\t("+sequence.getLength()+"+"+seq_len+")");
        }
            

        for (int i=0; i<seq_len; i++)
            sequence.addNoTest(char2byte[Character.toUpperCase(seq.charAt(i))]);
        //System.out.println("#*DNA.a "+seq_len+"\ttot "+sequence.getLength()+"\t"+seq);
    }
    
    /**
     * Counts the occurrences of a nucleotide
     * @param nuc nucleotide in byte=-encoding
     * @param start_pos 0-based start pos, inclusive
     * @param end_pos 0-based end pos, exclusive
     * @return how many times that nucleotide occurs 
     */
    public int countNucleotide(byte nuc, int start_pos, int end_pos)
    {
        return sequence.countOccurrences(nuc, start_pos, end_pos);
    }
    
//    public DNASequence getIntervalView(int start, int end)
//    {
//        return new IntervalSequence(this, start, end);
//    }
//    
//    private static class IntervalSequence extends DNASequence
//    {
//        private IntervalSequence(DNASequence parent, int start, int end)
//        {
//            super(parent);
//            this.start = start;
//            this.end = end;
//            this.parent = parent;
//        }
//        private final int start;
//        private final int end;
//        private final DNASequence parent;
//        
//        @Override 
//        public byte getNucleotideAt(int pos)
//        {
//            return parent.getNucleotideAt(pos+start);
//        }
//
//        @Override 
//        public int getLength()
//        {
//            return end-start;
//        }        
//        
//    }
        
    /**
     * Computes a string representation for the nucleotide sequence.
     * If start is larger than end, then reverse complement is calculated.
     *
     * @param start 0-based start position [inclusive]
     * @param end 0-based end position [exclusive]
     * @return String for the denoted interval. In Genbank flatfile 1-based closed-interval 
     *        coordinates, it is either the string start+1..end or complement(end+1..start)
     */
    public String getSubstring(int start, int end)
    {
        byte[] seq_slice;
        if (sequence == null)
            return "%%null%%";
        if (start<end)
            seq_slice = sequence.slice(start,end-start);
        else{
            seq_slice = sequence.slice(end,start-end);
            int i=0, j=seq_slice.length-1;
            while (i<=j)
            {
                byte ni = seq_slice[i];
                byte nj = seq_slice[j];
                seq_slice[j]=complements[ni];
                seq_slice[i]=complements[nj];
                ++i; --j;
            }
        }
        for (int i=0; i<seq_slice.length; i++)
        {
            byte b=seq_slice[i];
            char c=byte2char[b];
            //System.out.println("#**DNAS.gS i "+i+" "+b+"->"+c);
            seq_slice[i]=(byte) c; // new String(byte[]) for "platfrom-dependent" instantiation
        }
        return new String(seq_slice);
    }    

    /**
     * Complements a byte-encoded nucleotide
     * 
     * @return the complemented nucleotide
     */
    public static final byte complement(byte b){
        return complements[b];
    }

    /**
     * Complements a char-encoded nucleotide
     * 
     * @return the complemented nucleotide
     */
    public static final char complement(char c){
        return byte2char[complements[char2byte[c]]];
    }

    /**
     * Watson-Crick base pairing.
     *
     * @param idx 0..3 index
     * @return the complemented nucleotide's index
     */
    public static final int complement(int idx)
    {
        return 3-idx;
    }

    /**
     * Computes the reverse complement of a DNA String. 
     * E.g., reverseComplement("AWATC") produces "GATWT".
     * 
     * @param seq DNA sequence to be complemented
     * @return reverse complement
     */ 
    public static String reverseComplement(String seq)
    {
        int seq_len=seq.length();
        char[] zz=new char[seq_len];
        for (int i=0; i<seq_len; i++){
            char c=seq.charAt(seq_len-i-1);
            //System.out.println("#**DNA.rC "+i+"/"+seq_len+" c="+c+"  // seq='"+seq+"'");
            zz[i]=complement(c);
        }
        return new String(zz);
    }
    
    public static final boolean isPurine(byte nuc)
    {
        return (nuc==BYTE_A || nuc==BYTE_G || nuc==BYTE_R);
    }
    
    public static final boolean isPurine(int nuc_idx)
    {
        return (nuc_idx==IDX_A || nuc_idx==IDX_G );
    }
    
    
    /**
     * Reads in a Fasta file with multiple DNA sequences.
     * Sequences names include the Fasta comment sign &gt;
     *
     * @param file_name path to the Fasta file: may be GZIP compressed, and URLs starting with ftp: or http: are recognized
     * @return array of sequences listed in the file
     * @throws IOException if the read fails
     *
     * @author Eric Barake &amp;
     */
    public static DNASequence[] readMultiFasta(String file_name) throws IOException
    {
        return readMultiFasta(file_name, null);
    }

    /**
     * Reads in a Fasta file with multiple DNA sequences.
     * Sequences names include the Fasta comment sign &gt;
     *
     * @param file_name path to the Fasta file: may be GZIP compressed, and URLs starting with ftp: or http: are recognized
     * @param wanted_chromosomes set of requested file names, null means all of them
     * @return array of sequences listed in the file
     * @throws IOException if the read fails
     *
     * @author Eric Barake &amp;
     */
    public static DNASequence[] readMultiFasta(String file_name, Set<String> wanted_chromosomes) throws IOException
    {
        ArrayList<DNASequence> sequences = new ArrayList<>();
        // ensuring a very general range of input file specifications: URLs, compression.
        BufferedReader fileReader    = new GeneralizedFileReader(file_name);
        String         line;
        String         id;
        //StringBuffer   contentBuffer = new StringBuffer();
        // current sequence is collected in a DNASequence object
        DNASequence current_seq = null;


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
                    if (wanted_chromosomes==null || wanted_chromosomes.contains(id))
                    {

                        // contentBuffer = new StringBuffer();
                        // start a new DNA sequence
                        current_seq = new DNASequence(id);
                        // add to our list immediately
                        sequences.add(current_seq);
                    } else
                    {
                        current_seq = null;
                    }
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

//        for (int sidx=0; sidx<Math.min(100,sequences.size()); sidx++)
//        {
//            DNASequence S = sequences.get(sidx);
//            System.out.println("#*DNA.rMF read ref"+sidx+"\t"+S);
//            for (int pos=0; pos<Math.min(12,S.getLength()); pos++)
//            {
//                if ((pos) % 50==0)
//                {
//                    System.out.print("\n#*DNA.rMF "+sidx+"\t"+(1+pos)+"\t|");
//                }
//                System.out.print(toChar(S.getNucleotideAt(pos)));
//            }
//            System.out.println();
//        }
//        if (sequences.getLength()>=100)
//        {
//            System.out.println("#*DNA.rmF 100.."+Integer.toString(sequences.getLength()-1)+"\t[...]");
//        }

        return sequences.toArray(new DNASequence[0]);
    }

    /**
     * Sorting by identifiers (by {@link #getIdent() }). 
     * @param o other DNA sequence
     * @return negative value if this identifier is lexicographically before the other identifier
     */
    @Override
    public int compareTo(DNASequence o) 
    {
        return this.getIdent().compareTo(o.getIdent());
    }


//    /**
//     * Sets an input reader by parsing the file path: if it looks like an URL,
//     * a URL connection is initiated; if it ends with <tt>gz</tt>, then
//     * it is uncompressed on the fly.
//     *
//     * @param file_name URL or file path name
//     * @return reader for the (possibly uncompressed file content)
//     * @throws IOException if URL access fails, or the file cannot be opened for reading
//     */
//    public static Reader guessReaderForInput(String file_name) throws IOException
//    {
//        java.io.InputStream base ;
//        if (file_name.startsWith("ftp:") || file_name.startsWith("http:"))
//        {
//            URL url = new URL(file_name);
//
//            base = url.openStream();
//        } else
//        {
//            base=new java.io.FileInputStream(file_name);
//        }
//
//        if (file_name.endsWith(".gz"))
//            base = new GZIPInputStream(base);
//        else if (file_name.endsWith("zip"))
//            base = new ZipInputStream(base);
//        return new java.io.InputStreamReader(base);
//    }


    /**
     *
     * A class implementing the functionality of the {@link java.util.Vector} and {@link java.util.ArrayList} classes,
     * for byte-type elements.
     * 
     * The data structure relies on an array of arrays.
     * Memory is allocated by blocks of 8Mb (or whatever is set by {@link #BLOCK_BITS},
     *   and the array of blocks expands automatically
     *   as it becomes necessary. For very short sequences (less than the block getLength),
     *   the requested allocation getLength is respected, but when more bytes or added, or
     *   larger initial capacity is requested, the block-based allocation policy is used.
     *   If the capacity needs to be increased, it will be increased by a multiple of the block getLength.
     *   The implementation has the advantage that there is almost no extra working memory used when the
     *   space needs to be expanded.
     *
     * @author Mikl&oacute;s Cs&#369;r&ouml;s
     *
     */
    private static class BigByteArray
    {
        /**
         * Constant determining the block getLength in a {@link BigByteArray}.
         */
        private static final int BLOCK_BITS=20; // 20=1Mb 21=2Mb 22=4Mb 23=8Mb 24=16Mb 25=32Mb 26=64Mb 27=128Mb
        /**
         * Integer value of block getLength (2<sup>{@link #BLOCK_BITS}</sup>).
         */
        private static final int BLOCK_SIZE=1<<BLOCK_BITS;
        /**
         * Bitmask selecting lower bits for offset within a block.
         */
        private static final int BLOCK_OFFSET=BLOCK_SIZE-1;

        /** Creates a new instance of BigByteArray */
        BigByteArray()
        {
            this(BLOCK_SIZE);
        }

        BigByteArray(int initialCapacity)
        {
            init(initialCapacity);
        }

        /**
         * The actual data storage: an array of blocks.
         */
        private byte[][] elementData;

        /**
         * Maximum number of elements before expansion is needed.
         */
        private int capacity;
        /**
         * Actual count of elements.
         */
        private int elementCount;

        /**
         * @return number of elements addded so far
         */
        int size(){return elementCount;}

        /**
         * Access to individual bytes in the data structure. For the sake of
 speed, it is not verified if index is within the allowable range
 (less than getLength()): the return value is a 0 byte
 in that case, or a {@link IndexOutOfBoundsException} exception is thrown.
         *
         * @return the element at the given position
         */
        byte get(int index)
        {
            int hi = (int)(index>>BLOCK_BITS);
            int lo = (int)(index & BLOCK_OFFSET);

//            if (hi>=elementData.length || lo>=elementData[hi].length)
//            {
//                throw new IllegalArgumentException("needs idx= "+index+"\thi "+hi+"\tlo "+lo+"\te[] "+elementData.length
//                        +(hi<elementData.length?"\te[][] "+elementData[hi].length:""));
//            }
//
            return elementData[hi][lo];
        }

        /**
         * Adds a new element at the end.
         * @param element value to be added
         */
        void add(byte element)
        {
            ensureCapacity(elementCount+1);
            addNoTest(element);
        }


         /* Access to indivudal bytes in the data structure. For the sake of
         * speed, it is not verified it is not verified if index is within the allowable range
         * (less than getLength()): the return value is a 0 byte
         * in that case, or a {@link IndexOutOfBoundsException} exception is thrown.
         *
         * @return the old element at the given position
         */
        void set(int index, byte element)
        {
            int hi = (int)(index>>BLOCK_BITS);
            int lo = (int)(index & BLOCK_OFFSET);
            elementData[hi][lo]=element;
        }

        /**
         * Adds a ne wbyte at the end, without checking for capacity
         *
         * @param element
         */
        private void addNoTest(byte element)
        {
            set(elementCount++,element);
        }


        /**
         * Allocates space in {@link #elementData};
         * initializes {@link #capacity} and {@link #elementCount} variables.
         */
        private void init(int initialCapacity)
        {
            if (initialCapacity < BLOCK_SIZE) // tiny sequence
            {
                elementData=new byte[1][(int)initialCapacity];
                capacity = initialCapacity;
            } else 
            {
                // round up to complete blocks
                int nblocks = (initialCapacity+BLOCK_SIZE-1) >>> BLOCK_BITS; // rounding up to BLOCK_SIZE
                assert (nblocks == (int)nblocks); // better not ask for longer than 2^51 ...

                elementData = new byte[(int)nblocks][BLOCK_SIZE];
                capacity = nblocks<<BLOCK_BITS;
            }
            elementCount=0;
        }
        
        /**
         * Computes a slice of the vector and returns it as an array. Does not check whether
         * the indexes surpass the element count (for speed).
         *
         * @param pos starting index for the slice
         * @param length length of the slice
         * @return a byte[] array of given length for the elements pos..pos+length-1
         */
        public byte[] slice(int pos, int length)
        {
            return slice(pos, length, new byte[length]);
        }

        /**
         * Computes a slice of the vector and returns it as an array. Does not check whether
         * the indexes surpass the element count (for speed). The elements of the slice are put 
         * in the array given as a parameter: possible elements beyond the length are not changed.
         *
         * @param pos starting index for the slice
         * @param length length of the slice
         * @param destination where the elements are copied (with indexes 0 to length-1)
         * @return destination
         */
        public byte[] slice(int pos, int length, byte[] destination)
        {
            int offset = pos & BLOCK_OFFSET;
            int hi = pos >> BLOCK_BITS;
            for (int ncopied=0; ncopied<length; )
            {
                int in_this_block = Math.min(BLOCK_SIZE-offset,length-ncopied);
                System.arraycopy(elementData[hi], offset, destination, ncopied, in_this_block);
                offset=0;
                hi++;
                ncopied+=in_this_block;
            }
            
            return destination;
        }
        
        int countOccurrences(byte b, int start_pos, int end_pos)
        {
            int nocc=0;
            int length = end_pos-start_pos;
            int offset = start_pos & BLOCK_OFFSET;
            int hi = start_pos >> BLOCK_BITS;
            for (int nseen=0; nseen<length; )
            {
                int in_this_block = Math.min(BLOCK_SIZE-offset,length-nseen);
                final byte[] block = elementData[hi];
                for (int j=0; j<in_this_block; j++)
                {
                    if (block[offset]==b) 
                    {
                        nocc++;
                    }
                    offset++;
                }
                offset=0;
                hi++;
                nseen+=in_this_block;
            }
            return nocc;
        }
        

        /**
         * Enlarges the underlying data structure when necessary.
         *
         * @param minCapacity minimum capacity after enlargement (it is rounded up to the nearest multiple of {@link #BLOCK_SIZE} if expansion is performed).
         */
        private boolean ensureCapacity(int minCapacity)
        {
            boolean retval=false;
            
            //System.out.println("#**DNAS.BBA.eC "+minCapacity+"\tcap "+capacity+"\tblock "+BLOCK_SIZE);
            
            if (minCapacity > capacity && capacity < BLOCK_SIZE)
            { // first align the first array entry up to a complete block
                int first_block_size = capacity;
                while (first_block_size<minCapacity)
                {
                    first_block_size += (first_block_size % 3==0?first_block_size/3:first_block_size/2);
                }
                first_block_size = first_block_size>BLOCK_SIZE?BLOCK_SIZE:first_block_size;
                byte[] newData=new byte[first_block_size];
                System.arraycopy(elementData[0],0,newData,0,(int)elementCount);
                assert (elementData.length==1); // either we have complete blocks, or just one truncated block
                elementData[0]=null;
                elementData[0]=newData;
                capacity = first_block_size;
                //System.out.println("#**DNAS.BBA.eC/1 "+minCapacity+"\tcap "+capacity);
                retval = true;
            }

            if (minCapacity>capacity)
            {
                // at this point, capacity is a multiple of BLOCK_SIZE 
                int nblocks = (minCapacity+BLOCK_SIZE-1) >>> BLOCK_BITS; // rounding up to BLOCK_SIZE: this is how many blocks are needed
                assert (nblocks == (int)nblocks); // better not ask for longer than 2^51 ...
                assert (nblocks*BLOCK_SIZE >= minCapacity); // arithmetics was ok

                int block_capacity = elementData.length;
//                System.out.println("#*DNA.BBA.eC req "+minCapacity+"\tcap "+capacity+";\tnblocks "+nblocks+"\tbc "+block_capacity);

                if (nblocks > block_capacity)
                {
                    // figure out how large it should be
                    do block_capacity += block_capacity; while (nblocks>block_capacity);
                    byte[][] newData = new byte[(int)block_capacity][];

                    // existing blocks are kept
                    System.arraycopy(elementData, 0, newData, 0, elementData.length);
                    elementData = newData;
                }
                // allocate new blocks as necessary: since requested capacity is larger than current capacity, at least one new block needs to be allocated
                int block_idx=(int)(capacity >>> BLOCK_BITS);
                do
                {
                    assert (elementData[block_idx]==null);
                    elementData[block_idx] =  new byte[BLOCK_SIZE];
                    capacity += BLOCK_SIZE;
                    block_idx++;
                } while (block_idx<nblocks);
                retval = true;
                //System.out.println("#**DNAS.BBA.eC/n "+minCapacity+"\tcap "+capacity);
            }
            return retval;
        }
    }

    /**
     * @return a String description of this Object: for short sequences, the complete sequence
     * is included in the description, otherwise only a reasonably sized prefix.
     */
    @Override
    public String toString()
    {
        StringBuilder sb=new StringBuilder();
        sb.append("DNA sequence of length ");
        sb.append(getLength());
        sb.append(" '");
        sb.append(getDefLine());
        sb.append("'");
        return sb.toString();
    }
    
//    /**
//     * Generic pretty print for fasta files.
//     * 
//     * @param chars_per_line characters per line 
//     * @param defline definition line, prefixed with <code>&gt;/code>
//     * @param sequence sequence to be printed, any characters are allowed
//     * @return formatted string  
//     */
//    public static String toFasta(int chars_per_line, String defline, String sequence)
//    {
//        StringBuilder sb = new StringBuilder(">");
//        sb.append(defline);
//        sb.append("\n");
//        int pos = 0;
//        int len = sequence.length();
//        while (pos<len)
//        {
//            
//            int chunk_length = Math.min(chars_per_line, len-pos);
//            String ss = sequence.substring(pos, pos+chunk_length);
//            sb.append(ss);
//            pos += chunk_length;
//            if (pos<len)
//                sb.append("\n");
//        }
//        return sb.toString();
//    }

    /**
     * Class for representing a set of chromosomes, searchable by chromosome identifiers
     * (using {@link DNASequence#getIdent() }).
     */
    public static class Genome implements Iterable<DNASequence>
    {
        /**
         * Instantiation with a set of chromosomes.
         * 
         * @param chromosomes array of chromosomes
         */
        public Genome(DNASequence[] chromosomes)
        {
            chromo_names = new HashMap<>();
            for (DNASequence seq: chromosomes)
            {
                //String name = seq.getDefLine();
                String id = seq.getIdent();
                //System.out.println("#**G.() '"+id+"' -> "+seq+"\t// "+name);
                chromo_names.put(id, seq);
            }
        }

        public void addAll(DNASequence[] chromosomes)
        {
            for (DNASequence seq: chromosomes)
            {
                //String name = seq.getName();
                String id = seq.getIdent();
                //System.out.println("#**G.() '"+id+"' -> "+seq+"\t// "+name);
                chromo_names.put(id, seq);
            }
        }

        private final Map<String, DNASequence> chromo_names;

        /**
         * Finds the chromosome with the given identifier
         * 
         * @param chromo_id chromosome identifier (matching {@link #getIdent() }
         * @return null if not there
         */
        public DNASequence getChromosome(String chromo_id)
        {
            return chromo_names.get(chromo_id); // null if not there
        }

        /**
         * Chromosomes in arbitrary order.
         * 
         * @return iterator over the chromosomes
         */
        @Override
        public Iterator<DNASequence> iterator() 
        {
            return chromo_names.values().iterator();
        }    
    }


    /**
     * Test code.
     *
     * @param args command-line arguments: fasta-file
     * @throws IOException if reading fails
     */
    public static void main(String[] args)  throws IOException
    {
        if (args.length==0) throw new IllegalArgumentException("Call with Fasta file name in argument");

        String fastaFile = args[0];
        DNASequence[] refSeq = readMultiFasta(fastaFile);
        long tot_len = 0L;
        long tot_cap = 0L;
        for (int i=0; i<refSeq.length; i++)
        {
            DNASequence s = refSeq[i];
            System.out.println(i+"\t"+s.getDefLine() + "---" + s.getLength()+" nucleotides ("+s.sequence.capacity+" bytes)");
            tot_len += s.getLength();
            tot_cap += s.sequence.capacity;
        }
        
        System.out.println("Loaded " + refSeq.length+ " sequence"+(refSeq.length>1?"s":"")+".");
        System.out.println("Total length "+tot_len+" nucleotides\t(using "+tot_cap+" bytes)");
    }
    
    
}

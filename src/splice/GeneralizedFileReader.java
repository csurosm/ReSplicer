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

import java.io.IOException;
import java.io.Reader;
import java.net.URL;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

/**
 *
 * Common interface to reading files and URLs, possibly compressed.
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class GeneralizedFileReader extends java.io.BufferedReader
{
    
    public GeneralizedFileReader(String path_name) throws IOException
    {
        super(guessReaderForInput(path_name));
    }
    /**
     * Sets an input reader by parsing the file path: if it looks like an URL,
     * a URL connection is initiated; if it ends with <tt>gz</tt>, then
     * it is uncompressed on the fly.
     *
     * @param file_name URL or file path name
     * @return reader for the (possibly uncompressed file content)
     * @throws IOException if URL access fails, or the file cannot be opened for reading
     */
    public static Reader guessReaderForInput(String file_name) throws IOException
    {
        java.io.InputStream base ;
        if (file_name.startsWith("ftp:") || file_name.startsWith("http:"))
        {
            URL url = new URL(file_name);

            base = url.openStream();
        } else
        {
            base=new java.io.FileInputStream(file_name);
        }

        if (file_name.endsWith(".gz"))
            base = new GZIPInputStream(base);
        else if (file_name.endsWith("zip"))
            base = new ZipInputStream(base);
        return new java.io.InputStreamReader(base);
    }
}